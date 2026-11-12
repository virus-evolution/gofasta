/*
Package pipeline processes a stream of records concurrently and hands the results on in input order.

Records are read one at a time, processed on several goroutines, and passed to a sink, which is called
from the goroutine that called Run, one record at a time, in the order the records were read. Every
record travels as a Result, which carries any error from reading or processing it, and every Result
reaches the sink (unless the pipeline is stopped), so the sink is the one place that decides what to do
with failed records: stop at the first (FirstError, or a strict Handler), or warn about them, skip them
and carry on (Handler).

Errors from reading (next) arrive as a *ReadError and end the input: a Handler always treats them as
fatal. Errors from processing a record (work) belong to that record; wrapping them with Named gives the
record's name to the warning.

The number of records between reading and the sink is bounded, so a record that is slow to process does
not let the others pile up without limit.
*/
package pipeline

import (
	"context"
	"errors"
	"fmt"
	"io"
	"runtime/debug"
	"sync"
	"sync/atomic"
)

// Result is one record on its way through a pipeline.
type Result[T any] struct {
	Idx int   // 0-based position of the record in the input
	Val T     // the processed record (the zero value if Err is set)
	Err error // why the record could not be read or processed, if it could not
}

// ReadError is the error of a Result whose record could not be read. Reading stops after it.
type ReadError struct {
	Err error
}

func (e *ReadError) Error() string { return e.Err.Error() }
func (e *ReadError) Unwrap() error { return e.Err }

// RecordError is the error of one record that could not be processed, with the record's name.
type RecordError struct {
	Name string
	Err  error
}

func (e *RecordError) Error() string { return e.Err.Error() }
func (e *RecordError) Unwrap() error { return e.Err }

// Named attaches a record's name to the error from processing it (nil stays nil).
func Named(name string, err error) error {
	if err == nil {
		return nil
	}
	return &RecordError{Name: name, Err: err}
}

// PanicError is the value Run panics with when next or work panicked. It holds the original panic value
// and the stack of the goroutine that panicked.
type PanicError struct {
	Value any
	Stack []byte
}

func (p *PanicError) Error() string {
	return fmt.Sprintf("%v\n\ngoroutine stack at the original panic:\n%s", p.Value, p.Stack)
}

// Unwrap returns the original panic value if it was an error.
func (p *PanicError) Unwrap() error {
	if err, ok := p.Value.(error); ok {
		return err
	}
	return nil
}

// inFlightPerThread bounds how many records may be between reading and the sink at once, per worker.
const inFlightPerThread = 8

// Run calls next until it returns io.EOF, passes each record to work on threads goroutines (at least
// one), and passes the results to sink in input order, from the calling goroutine.
//
// An error from next becomes a Result with that error (as a *ReadError), at the position of the record
// that could not be read, and reading stops there. An error from work becomes the Err of that record's Result.
//
// Run stops early if sink returns an error, which Run then returns, or if ctx is cancelled before the end
// of the input has been read and every record has reached the sink, in which case it returns ctx.Err()
// (the sink is not called again once ctx is cancelled). Either way it waits for its workers to finish before returning; next may still
// be in the middle of one call, whose result is discarded. Functions that should notice cancellation
// themselves can close over ctx.
//
// If next or work panics, Run stops, waits for its workers, and panics in the calling goroutine with a
// *PanicError. (A panic in a call of next that is still running when Run returns for another reason is
// discarded.)
func Run[In, Out any](ctx context.Context, threads int, next func() (In, error), work func(In) (Out, error), sink func(Result[Out]) error) error {
	if threads < 1 {
		threads = 1
	}
	inner, cancel := context.WithCancel(ctx)
	defer cancel()

	// the first panic in next or work; recordPanic must be deferred directly by the goroutine that may panic
	var panicked atomic.Pointer[PanicError]
	recordPanic := func() {
		if v := recover(); v != nil {
			panicked.CompareAndSwap(nil, &PanicError{Value: v, Stack: debug.Stack()})
			cancel()
		}
	}

	type job struct {
		idx int
		in  In
		err error
	}
	jobs := make(chan job)
	results := make(chan Result[Out], threads)
	tokens := make(chan struct{}, threads*inFlightPerThread)

	// The reader sets exhausted once it has handed on its last record (at io.EOF or a read error), and
	// sent counts the records it handed on, so that Run can tell whether everything reached the sink.
	var (
		exhausted atomic.Bool
		sent      atomic.Int64
	)

	// reader
	go func() {
		defer close(jobs)
		defer recordPanic()
		for i := 0; ; i++ {
			select {
			case tokens <- struct{}{}:
			case <-inner.Done():
				return
			}
			in, err := next()
			if err == io.EOF {
				exhausted.Store(true)
				return
			}
			if err != nil {
				err = &ReadError{Err: err}
			}
			select {
			case jobs <- job{idx: i, in: in, err: err}:
				sent.Add(1)
			case <-inner.Done():
				return
			}
			if err != nil {
				exhausted.Store(true)
				return
			}
		}
	}()

	// workers
	var wg sync.WaitGroup
	wg.Add(threads)
	for w := 0; w < threads; w++ {
		go func() {
			defer wg.Done()
			defer recordPanic()
			for {
				var j job
				var ok bool
				select {
				case j, ok = <-jobs:
					if !ok {
						return
					}
				case <-inner.Done():
					return
				}
				r := Result[Out]{Idx: j.idx, Err: j.err}
				if j.err == nil {
					r.Val, r.Err = work(j.in)
				}
				select {
				case results <- r:
				case <-inner.Done():
					return
				}
			}
		}()
	}
	go func() {
		wg.Wait()
		close(results)
	}()

	// in-order delivery
	pending := make(map[int]Result[Out])
	nextIdx := 0
	var sinkErr error
	for r := range results {
		if sinkErr != nil || inner.Err() != nil {
			continue
		}
		pending[r.Idx] = r
		for inner.Err() == nil {
			p, ok := pending[nextIdx]
			if !ok {
				break
			}
			delete(pending, nextIdx)
			nextIdx++
			<-tokens
			if sinkErr = sink(p); sinkErr != nil {
				cancel()
				break
			}
		}
	}

	if p := panicked.Load(); p != nil {
		panic(p)
	}
	if sinkErr != nil {
		return sinkErr
	}
	if !exhausted.Load() || int64(nextIdx) != sent.Load() {
		return ctx.Err()
	}
	return nil
}

// FirstError is a sink adapter: it returns a Result's error if it has one, and otherwise passes the
// value to f. It gives a pipeline the policy of stopping at the first record that failed.
func FirstError[T any](f func(T) error) func(Result[T]) error {
	return func(r Result[T]) error {
		if r.Err != nil {
			return r.Err
		}
		return f(r.Val)
	}
}

// Noun names the records of a pipeline in a Handler's warnings, e.g. Noun{"record", "records"}.
type Noun struct {
	One, Many string
}

// Handler is a sink that applies gofasta's policy for records that fail. If strict, the first failure
// stops the pipeline with that record's error. Otherwise each failed record is reported on warn, as
//
//	warning: skipping NAME (NOUN N): ERROR
//
// and skipped, and the others are passed to f. After Run, Finish reports how many were skipped. Read
// errors (*ReadError) always stop the pipeline.
type Handler[T any] struct {
	strict  bool
	warn    io.Writer
	noun    Noun
	f       func(T) error
	total   int
	skipped int
}

// NewHandler returns a Handler that passes the records that did not fail to f.
func NewHandler[T any](strict bool, warn io.Writer, noun Noun, f func(T) error) *Handler[T] {
	return &Handler[T]{strict: strict, warn: warn, noun: noun, f: f}
}

// Sink is the sink to pass to Run.
func (h *Handler[T]) Sink(r Result[T]) error {
	if r.Err == nil {
		h.total++
		return h.f(r.Val)
	}
	var readErr *ReadError
	if h.strict || errors.As(r.Err, &readErr) {
		return r.Err
	}
	h.total++
	h.skipped++
	name := ""
	var recErr *RecordError
	if errors.As(r.Err, &recErr) && recErr.Name != "" {
		name = recErr.Name + " "
	}
	_, err := fmt.Fprintf(h.warn, "warning: skipping %s(%s %d): %v\n", name, h.noun.One, r.Idx+1, r.Err)
	return err
}

// Total is the number of records that reached the sink (not counting a read error).
func (h *Handler[T]) Total() int { return h.total }

// Skipped is the number of records that failed and were skipped.
func (h *Handler[T]) Skipped() int { return h.skipped }

// Finish writes a summary if any records were skipped, and returns an error if there were records and
// every one of them failed.
func (h *Handler[T]) Finish() error {
	if h.skipped == 0 {
		return nil
	}
	if _, err := fmt.Fprintf(h.warn, "warning: skipped %d of %d %s\n", h.skipped, h.total, h.noun.Many); err != nil {
		return err
	}
	if h.skipped == h.total {
		return fmt.Errorf("none of the %d %s could be processed", h.total, h.noun.Many)
	}
	return nil
}
