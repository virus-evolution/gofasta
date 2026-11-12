/*
Package pipeline processes a stream of records concurrently and hands the results on in input order.

Records are read one at a time, processed on several goroutines, and passed to a sink, which is called
from the goroutine that called Run, one record at a time, in the order the records were read. Every
record travels as a Result, which carries any error from reading or processing it, and every Result
reaches the sink (unless the pipeline is stopped), so the sink is the one place that decides what to do
with failed records: stop at the first, as FirstError does, or skip them and carry on.

The number of records between reading and the sink is bounded, so a record that is slow to process does
not let the others pile up without limit.
*/
package pipeline

import (
	"context"
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
// An error from next becomes a Result with that error, at the position of the record that could not be
// read, and reading stops there. An error from work becomes the Err of that record's Result.
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
