package sam

import (
	"io"
	"strings"
	"testing"
)

// Malformed sam input gives an error (not a hang, and not a panic in a worker goroutine).
func TestMalformedSamGivesError(t *testing.T) {
	hdr := "@SQ\tSN:r\tLN:45\n"
	rec := func(pos, cigar, seq string) string {
		return "q1\t0\tr\t" + pos + "\t60\t" + cigar + "\t*\t0\t0\t" + seq + "\t*\n"
	}
	for name, in := range map[string]string{
		"empty input":               "",
		"not sam":                   "this is not a sam file\n",
		"no @SQ":                    "@HD\tVN:1.6\n" + rec("1", "10M", "ATGAAACCCG"),
		"malformed header":          "@SQ\t\tr",
		"CIGAR ends in a number":    hdr + rec("1", "5M5", "ATGAAACCCG"),
		"SEQ * with a CIGAR":        hdr + rec("1", "10M", "*"),
		"CIGAR shorter than SEQ":    hdr + rec("1", "5M", "ATGAAACCCG"),
		"CIGAR longer than SEQ":     hdr + rec("1", "20M", "ATGAAACCCG"),
		"past the end of reference": hdr + rec("40", "10M", "ATGAAACCCG"),
		"deletion past the end":     hdr + rec("40", "5M10D", "ATGAA"),
		"second record is bad":      hdr + rec("1", "10M", "ATGAAACCCG") + "q2\t0\tr\t1\t60\t10M\t*\t0\t0\tA\t*\n",
	} {
		quietStderr(t)
		check := func(fn string, err error) {
			if err == nil {
				t.Errorf("%s (%s): expected an error", name, fn)
			}
		}
		runWithTimeout(t, func() {
			check("toMultiAlign", ToMultiAlign(strings.NewReader(in), io.Discard, 60, -1, -1, false, 2, true))
		})
		runWithTimeout(t, func() {
			discardStdout(t)
			check("toPairAlign", ToPairAlign(strings.NewReader(in), strings.NewReader(">r\n"+fuzzRef+"\n"), "stdout", 60, -1, -1, false, false, 2, true))
		})
		runWithTimeout(t, func() {
			check("variants", Variants(strings.NewReader(in), nil, false, strings.NewReader(fuzzGenbank), "gb", io.Discard, -1, -1, false, 0, false, 2, true))
		})
	}
}

func TestReferenceLengthMustMatchSamHeader(t *testing.T) {
	in := "@SQ\tSN:r\tLN:99\nq1\t0\tr\t1\t60\t10M\t*\t0\t0\tATGAAACCCG\t*\n"
	discardStdout(t)
	runWithTimeout(t, func() {
		if err := ToPairAlign(strings.NewReader(in), strings.NewReader(">r\n"+fuzzRef+"\n"), "stdout", 60, -1, -1, false, false, 2, true); err == nil {
			t.Errorf("toPairAlign: expected an error for a reference of the wrong length")
		}
		if err := Variants(strings.NewReader(in), nil, false, strings.NewReader(fuzzGenbank), "gb", io.Discard, -1, -1, false, 0, false, 2, true); err == nil {
			t.Errorf("variants: expected an error for a reference of the wrong length")
		}
	})
}

// records mapped to a reference other than the first @SQ are skipped; the other records are processed as usual
func TestRecordsOnOtherReferencesAreSkipped(t *testing.T) {
	in := "@SQ\tSN:r\tLN:45\n@SQ\tSN:s\tLN:90\n" +
		"q1\t0\ts\t50\t60\t10M\t*\t0\t0\tATGAAACCCG\t*\n" +
		"q2\t0\tr\t1\t60\t45M\t*\t0\t0\t" + fuzzRef + "\t*\n" +
		"q2\t2048\ts\t1\t60\t10M\t*\t0\t0\tATGAAACCCG\t*\n"
	quietStderr(t)
	out := new(strings.Builder)
	runWithTimeout(t, func() {
		if err := ToMultiAlign(strings.NewReader(in), out, -1, -1, -1, false, 2, true); err != nil {
			t.Errorf("unexpected error: %v", err)
		}
	})
	if want := ">q2\n" + fuzzRef + "\n"; out.String() != want {
		t.Errorf("got %q, want %q", out.String(), want)
	}
}
