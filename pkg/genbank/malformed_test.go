package genbank

import (
	"strings"
	"testing"
)

// Malformed input gives an error (or is tolerated), never a panic or a huge allocation.
func TestMalformedLocationsDoNotPanic(t *testing.T) {
	for _, s := range []string{
		"", " ", "comp0", "(join())", "complement(join(1,5..9.))", "join(1..3,5)", "join()", "join(,)", "complement(",
		"join(1..3,", "1..122222222222222222", "join(1..99999999999,2..3)", "9..1", "0..3", "-1..5", "1..2..3",
		"complement()", "join(complement(),1..3)", strings.Repeat("join(", 1000) + "1..3" + strings.Repeat(")", 1000),
	} {
		l := Location{Representation: s}
		if _, err := l.GetPositions(); err == nil {
			t.Errorf("%q: expected an error", s)
		}
		if _, err := l.IsReverse(); err == nil {
			t.Errorf("%q: expected an error from IsReverse", s)
		}
	}
}

func TestNestedLocations(t *testing.T) {
	for s, want := range map[string][]int{
		"join(complement(1..3),5..6)":             {3, 2, 1, 5, 6},
		"complement(join(1..3,5..6))":             {6, 5, 3, 2, 1},
		"join(1..2,join(4..5,7..8))":              {1, 2, 4, 5, 7, 8},
		"complement(join(1..2,complement(4..5)))": {4, 5, 2, 1},
	} {
		got, err := Location{Representation: s}.GetPositions()
		if err != nil || len(got) != len(want) {
			t.Errorf("%q: got %v, %v", s, got, err)
			continue
		}
		for i := range want {
			if got[i] != want[i] {
				t.Errorf("%q: got %v want %v", s, got, want)
			}
		}
	}
}

func TestMalformedFeaturesDoNotPanic(t *testing.T) {
	for name, in := range map[string]string{
		"whitespace-only line":      "LOCUS X\nFEATURES             Location/Qualifiers\n     CDS             1..9\n   \n                     /gene=\"a\"\nORIGIN\n        1 atgaaataa\n//\n",
		"qualifier before feature":  "LOCUS X\nFEATURES\n                     /gene=\"a\"\n     CDS             1..9\nORIGIN\n",
		"no features, only a blank": "LOCUS X\nFEATURES\n  \nORIGIN\n",
	} {
		_, err := ReadGenBank(strings.NewReader(in))
		if name == "qualifier before feature" && err == nil {
			t.Errorf("%s: expected an error", name)
		}
		if name != "qualifier before feature" && err != nil {
			t.Errorf("%s: unexpected error %v", name, err)
		}
	}
}
