package snps

// Regression tests.

import (
	"bytes"
	"strings"
	"testing"
)

// with --hard-gaps, a gap in both the reference and the query is not a SNP.
func TestSNPsHardGapsGapVsGap(t *testing.T) {
	ref := strings.NewReader(">ref\nAC-T\n")
	query := strings.NewReader(">q1\nAC-T\n>q2\nAA-G\n")

	out := new(bytes.Buffer)
	if err := SNPs(ref, query, true, false, 0.0, out, true); err != nil {
		t.Fatal(err)
	}

	// q1 is identical to the reference. q2 differs at positions 2 and 4 only
	// (the shared gap at position 3 is not a difference).
	want := "query,SNPs\nq1,\nq2,C2A|T4G\n"
	if out.String() != want {
		t.Errorf("want:\n%s\ngot:\n%s", want, out.String())
	}
}
