package closest

// Regression tests.

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

// a target with no comparable sites has a NaN distance, which must never be reported as the closest,
// even when it is listed first.
func TestClosestNaNTargetIsNeverReplaced(t *testing.T) {
	targetData := ">allN\nNNNNNNNNNN\n>near\nACGTACGTAA\n"
	queryData := ">q\nACGTACGTAC\n"

	for _, measure := range []string{"raw", "tn93"} {
		out := new(bytes.Buffer)
		err := Closest(strings.NewReader(queryData), strings.NewReader(targetData), measure, out, 1, true)
		if err != nil {
			t.Fatal(err)
		}

		lines := strings.Split(strings.TrimSpace(out.String()), "\n")
		if len(lines) != 2 {
			t.Fatalf("%s: unexpected output:\n%s", measure, out.String())
		}
		fields := strings.Split(lines[1], ",")
		if fields[1] != "near" {
			t.Errorf("%s: closest should be 'near' (finite distance), got %q (full line: %s)", measure, fields[1], lines[1])
		}
		if strings.Contains(lines[1], "NaN") {
			t.Errorf("%s: output should never contain NaN: %s", measure, lines[1])
		}
	}
}

// a NaN distance must not pass the --max-dist filter in ClosestN.
func TestClosestNNaNBypassesMaxDist(t *testing.T) {
	targetData := ">allN\nNNNNNNNNNN\n>near\nACGTACGTAA\n"
	queryData := ">q\nACGTACGTAC\n"

	out := new(bytes.Buffer)
	err := ClosestN(0, 0.5, strings.NewReader(queryData), strings.NewReader(targetData), "raw", out, false, 1, true)
	if err != nil {
		t.Fatal(err)
	}

	if strings.Contains(out.String(), "allN") {
		t.Errorf("a target with an undefined distance should not pass the --max-dist filter:\n%s", out.String())
	}
}

func TestClosestNNegativeNumberIsAnError(t *testing.T) {
	err := ClosestN(-3, 5, strings.NewReader(">a\nACGT\n"), strings.NewReader(">b\nACGT\n"), "raw", io.Discard, false, 1, true)
	if err == nil {
		t.Errorf("expected an error")
	}
}
