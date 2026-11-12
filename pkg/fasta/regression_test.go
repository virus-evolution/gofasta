package fasta

// Regression tests.

import (
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
