package variants

// Regression tests.

import (
	"bytes"
	"reflect"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/fasta"
)

func encodeForTest(t *testing.T, s string) []byte {
	t.Helper()
	enc, err := fasta.Record{Seq: s}.Encode()
	if err != nil {
		t.Fatal(err)
	}
	return enc.Seq
}

// insertion positions are in reference coordinates, also after earlier insertions in the reference
// (which are gap columns in the alignment).
func TestInsertionPositionAfterEarlierInsertion(t *testing.T) {
	// reference has two 2-base gaps: one after reference base 10 and one after base 20
	ref := encodeForTest(t, "ACGTACGTAC--GTACGTACGT--ACGTACGTAC")
	que := encodeForTest(t, "ACGTACGTACAAGTACGTACGTCCACGTACGTAC")

	indels := getIndelsPair(ref, que)

	want := []Variant{
		{Position: 10, Changetype: "ins", Length: 2},
		{Position: 20, Changetype: "ins", Length: 2},
	}

	if !reflect.DeepEqual(want, indels) {
		t.Errorf("insertion positions should be in reference coordinates\nwant: %+v\ngot:  %+v", want, indels)
	}
}

// an insertion at the very end of the alignment follows the same convention as one in the middle:
// its position is the number of reference bases before it.
func TestInsertionPositionAtEndOfAlignment(t *testing.T) {
	ref := encodeForTest(t, "ACGTACGTAC--")
	que := encodeForTest(t, "ACGTACGTACAA")

	indels := getIndelsPair(ref, que)

	want := []Variant{{Position: 10, Changetype: "ins", Length: 2}}

	if !reflect.DeepEqual(want, indels) {
		t.Errorf("end-of-alignment insertion should follow the same convention as mid-sequence insertions\nwant: %+v\ngot:  %+v", want, indels)
	}
}

// The order of --aggregate output must be deterministic: mutations that tie on position, type and
// alternative allele (e.g. deletions of different lengths at the same position) are ordered by their
// remaining fields.
func TestAggregateOutputOrderIsDeterministic(t *testing.T) {
	vs := []Variant{
		{Changetype: "del", Position: 100, Length: 12},
		{Changetype: "del", Position: 100, Length: 9},
		{Changetype: "del", Position: 100, Length: 3},
		{Changetype: "ins", Position: 100, Length: 6},
		{Changetype: "ins", Position: 100, Length: 3},
		{Changetype: "aa", Position: 100, Feature: "S", RefAl: "A", QueAl: "G", Residue: 5},
		{Changetype: "aa", Position: 100, Feature: "N", RefAl: "A", QueAl: "G", Residue: 9},
		{Changetype: "aa", Position: 100, Feature: "N", RefAl: "A", QueAl: "G", Residue: 3},
		{Changetype: "aa", Position: 100, Feature: "N", RefAl: "C", QueAl: "G", Residue: 3},
	}

	want := "mutation,frequency\n" +
		"aa:N:A3G,1.000000000\n" +
		"aa:N:C3G,1.000000000\n" +
		"aa:N:A9G,1.000000000\n" +
		"aa:S:A5G,1.000000000\n" +
		"del:100:3,1.000000000\n" +
		"del:100:9,1.000000000\n" +
		"del:100:12,1.000000000\n" +
		"ins:100:3,1.000000000\n" +
		"ins:100:6,1.000000000\n"

	for i := 0; i < 50; i++ {
		out := new(bytes.Buffer)
		agg, err := NewAggregator(out, -1, -1, false, "ref")
		if err != nil {
			t.Fatal(err)
		}
		if err := agg.Add(AnnoStructs{Queryname: "q", Vs: vs}); err != nil {
			t.Fatal(err)
		}
		if err := agg.Write(0.0); err != nil {
			t.Fatal(err)
		}

		if out.String() != want {
			t.Fatalf("run %d: output order is not deterministic/as expected\nwant:\n%s\ngot:\n%s", i, want, out.String())
		}
	}
}

// deletion positions are in reference coordinates, also after insertions in the reference (gap columns),
// and a deletion at the very start of the alignment is not reported
func TestDeletionPositionAfterEarlierInsertion(t *testing.T) {
	ref := encodeForTest(t, "ACGTACGTAC--GTACGTACGTACGTACGTAC")
	que := encodeForTest(t, "--GTACGTACAAGTACG---CGTACG-ACGTAC")[:32]
	indels := getIndelsPair(ref, que)

	want := []Variant{
		{Position: 10, Changetype: "ins", Length: 2},
		{Position: 16, Changetype: "del", Length: 3},
		{Position: 25, Changetype: "del", Length: 1},
	}
	if !reflect.DeepEqual(want, indels) {
		t.Errorf("want: %+v\ngot:  %+v", want, indels)
	}
}
