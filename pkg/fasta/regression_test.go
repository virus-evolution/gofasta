package fasta

// Regression tests.

import (
	"errors"
	"io"
	"strings"
	"testing"
)

// a header on the last line with no trailing newline is returned like one that ends in a newline:
// as a record with an empty sequence.
func TestFinalHeaderWithoutNewlineIsReturned(t *testing.T) {
	for _, input := range []string{">a\nAC\n>b", ">a\nAC\n>b\n"} {
		r := NewReader(strings.NewReader(input))

		rec, err := r.Read()
		if err != nil || rec.ID != "a" || rec.Seq != "AC" {
			t.Fatalf("%q: first record should be read normally, got %+v, %v", input, rec, err)
		}

		rec, err = r.Read()
		if err != nil || rec.ID != "b" || rec.Seq != "" {
			t.Errorf("%q: expected record b with an empty sequence, got %+v, %v", input, rec, err)
		}

		_, err = r.Read()
		if err != io.EOF {
			t.Errorf("%q: expected io.EOF after the last record, got %v", input, err)
		}
	}
}

func TestEmptyHeaderDoesNotPanic(t *testing.T) {
	for _, in := range []string{">", ">\n", ">\nAC\n", "> \nAC\n>x\nGT\n", ">\r\nAC\r\n"} {
		r := NewReader(strings.NewReader(in))
		rec, err := r.Read()
		if err != nil || rec.ID != "" {
			t.Errorf("%q: got ID %q, err %v", in, rec.ID, err)
		}
	}
}

// a read error part way through a record is returned as an error
type failingReader struct{ n int }

func (f *failingReader) Read(p []byte) (int, error) {
	if f.n == 0 {
		return 0, errors.New("boom")
	}
	f.n--
	return copy(p, ">a\n"), nil
}

func TestReadErrorMidRecord(t *testing.T) {
	r := NewReader(&failingReader{n: 1})
	if _, err := r.Read(); err == nil {
		t.Errorf("expected the read error")
	}
}
