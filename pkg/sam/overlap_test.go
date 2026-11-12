package sam

import (
	"fmt"
	"math/rand"
	"strconv"
	"strings"
	"testing"

	biogosam "github.com/biogo/hts/sam"
)

// Where the records of one query overlap in reference coordinates, the alignment of the FIRST record is
// used for the overlapping stretch and the later records only contribute outside it.

func overlapTestRef() string {
	rng := rand.New(rand.NewSource(7))
	b := make([]byte, 60)
	for i := range b {
		b[i] = "ACGT"[rng.Intn(4)]
	}
	return string(b)
}

func samLines(ref string, lines ...string) string {
	return "@SQ\tSN:ref\tLN:" + strconv.Itoa(len(ref)) + "\n" + strings.Join(lines, "\n") + "\n"
}

func rec(name string, flag int, pos int, cigar string, seq string) string {
	return fmt.Sprintf("%s\t%d\tref\t%d\t60\t%s\t*\t0\t0\t%s\t*", name, flag, pos, cigar, seq)
}

// The same insertion carried by both of two overlapping records is ONE insertion, not two.
func TestOverlappingRecordsWithTheSameInsertion(t *testing.T) {
	ref := overlapTestRef()
	ins := "CCC"

	single := pairAlignFromSam(t, samLines(ref,
		rec("q", 0, 1, "20M3I40M", ref[:20]+ins+ref[20:])), ref)

	// piece 1 covers reference 1-40, piece 2 (supplementary) covers 11-60; both carry the insertion after 20
	two := pairAlignFromSam(t, samLines(ref,
		rec("q", 0, 1, "20M3I20M", ref[:20]+ins+ref[20:40]),
		rec("q", 2048, 11, "10M3I40M", ref[10:20]+ins+ref[20:])), ref)

	if string(two.ref) != string(single.ref) || string(two.query) != string(single.query) {
		t.Errorf("two overlapping records should give the alignment of a single record\nwant ref %s\n     que %s\ngot  ref %s\n     que %s",
			single.ref, single.query, two.ref, two.query)
	}
}

// Two equally valid placements of an insertion in a repeat (here a tandem duplication that can be placed
// after reference base 17 or after base 20) are the same event, and the first record's placement is used.
func TestOverlappingRecordsWithEquivalentInsertionPlacements(t *testing.T) {
	ref := overlapTestRef()
	dup := ref[17:20]

	// the query is the reference with an extra copy of ref[17:20]
	single := pairAlignFromSam(t, samLines(ref,
		rec("q", 0, 1, "17M3I43M", ref[:17]+dup+ref[17:])), ref)

	two := pairAlignFromSam(t, samLines(ref,
		rec("q", 0, 1, "17M3I23M", ref[:17]+dup+ref[17:40]),      // placed after base 17
		rec("q", 2048, 11, "10M3I40M", ref[10:20]+dup+ref[20:])), // placed after base 20
		ref)

	if string(two.ref) != string(single.ref) || string(two.query) != string(single.query) {
		t.Errorf("want the first record's placement\nwant ref %s\n     que %s\ngot  ref %s\n     que %s",
			single.ref, single.query, two.ref, two.query)
	}
}

// An insertion at the very end of the earlier record's alignment is not inside the earlier record's
// coverage, so it is kept.
func TestOverlappingRecordsKeepInsertionAtTheEndOfTheEarlierRecord(t *testing.T) {
	ref := overlapTestRef()
	ins := "GGG"

	two := pairAlignFromSam(t, samLines(ref,
		rec("q", 0, 1, "20M", ref[:20]),
		rec("q", 2048, 11, "10M3I40M", ref[10:20]+ins+ref[20:])), ref)

	wantRef := ref[:20] + "---" + ref[20:]
	wantQue := ref[:20] + ins + ref[20:]
	if string(two.ref) != wantRef || string(two.query) != wantQue {
		t.Errorf("want ref %s\n     que %s\ngot  ref %s\n     que %s", wantRef, wantQue, two.ref, two.query)
	}
}

// A record that lies entirely inside the earlier record's coverage adds nothing, and a later record that
// spans an earlier one still contributes outside it.
func TestOverlappingRecordsContainment(t *testing.T) {
	ref := overlapTestRef()

	// the later record is entirely inside the first: its (different) bases are ignored
	other := strings.Repeat("A", 20)
	inside := pairAlignFromSam(t, samLines(ref,
		rec("q", 0, 1, "60M", ref),
		rec("q", 2048, 21, "20M", other)), ref)
	if string(inside.query) != ref {
		t.Errorf("a record inside the first one's coverage should be ignored\nwant %s\ngot  %s", ref, inside.query)
	}

	// the first record is inside the later one: the later one fills in around it, the first one's bases win
	mid := ref[20:40]
	flipped := strings.Repeat("T", 20)
	around := pairAlignFromSam(t, samLines(ref,
		rec("q", 0, 21, "20M", mid),
		rec("q", 2048, 1, "60M", ref[:20]+flipped+ref[40:])), ref)
	want := ref[:20] + mid + ref[40:]
	if string(around.query) != want {
		t.Errorf("want\n%s\ngot\n%s", want, around.query)
	}
}

// resolveOverlaps keeps every record valid (CIGAR and SEQ lengths agree), and afterwards no reference
// position is aligned by more than one record, every position that was covered still is, and each position
// still carries the base from the first record that covered it.
func TestResolveOverlapsRandomised(t *testing.T) {
	rng := rand.New(rand.NewSource(3))
	const L = 200

	for iter := 0; iter < 300; iter++ {
		ref := make([]byte, L)
		for i := range ref {
			ref[i] = "ACGT"[rng.Intn(4)]
		}
		var recs []biogosam.Record
		for k := 0; k < 1+rng.Intn(4); k++ {
			pos := rng.Intn(L - 20)
			var cigar []string
			var seq []byte
			rpos := pos
			for rpos < L && rng.Intn(8) != 0 {
				n := 1 + rng.Intn(15)
				switch rng.Intn(4) {
				case 0: // insertion, not at the start
					if len(cigar) > 0 {
						cigar = append(cigar, fmt.Sprintf("%dI", n))
						for j := 0; j < n; j++ {
							seq = append(seq, "ACGT"[rng.Intn(4)])
						}
					}
				case 1: // deletion, not at the start
					if len(cigar) > 0 && rpos+n < L {
						cigar = append(cigar, fmt.Sprintf("%dD", n))
						rpos += n
					}
				default:
					if rpos+n > L {
						n = L - rpos
					}
					cigar = append(cigar, fmt.Sprintf("%dM", n))
					for j := 0; j < n; j++ {
						seq = append(seq, "ACGT"[rng.Intn(4)]) // the records disagree about the bases
					}
					rpos += n
				}
			}
			if len(cigar) == 0 || strings.HasSuffix(cigar[len(cigar)-1], "I") || strings.HasSuffix(cigar[len(cigar)-1], "D") {
				continue
			}
			text := "@SQ\tSN:ref\tLN:" + strconv.Itoa(L) + "\n" +
				fmt.Sprintf("q\t%d\tref\t%d\t60\t%s\t*\t0\t0\t%s\t*\n", 0, pos+1, strings.Join(cigar, ""), seq)
			r, err := biogosam.NewReader(strings.NewReader(text))
			if err != nil {
				t.Fatal(err)
			}
			rr, err := r.Read()
			if err != nil {
				t.Fatalf("%v\n%s", err, text)
			}
			recs = append(recs, *rr)
		}
		if len(recs) == 0 {
			continue
		}

		out := resolveOverlaps(append([]biogosam.Record(nil), recs...))

		// the owner of every reference position, from the original records
		owner := make([]byte, L)
		for _, r := range recs {
			start := r.Pos
			for i, b := range referenceView(r) {
				if owner[start+i] == 0 && b != 0 {
					owner[start+i] = b
				}
			}
		}
		got := make([]byte, L)
		for _, r := range out {
			_, qlen := r.Cigar.Lengths()
			hard := 0
			for _, op := range r.Cigar {
				if op.Type() == biogosam.CigarHardClipped {
					hard += op.Len()
				}
			}
			if qlen-hard != len(r.Seq.Expand()) {
				t.Fatalf("iteration %d: CIGAR %v and SEQ (%d) disagree", iter, r.Cigar, len(r.Seq.Expand()))
			}
			for i, b := range referenceView(r) {
				if b == 0 {
					continue
				}
				if got[r.Pos+i] != 0 {
					t.Fatalf("iteration %d: reference position %d is aligned by more than one record", iter, r.Pos+i)
				}
				got[r.Pos+i] = b
			}
		}
		if string(got) != string(owner) {
			t.Fatalf("iteration %d: the alignment of the first record at each position was not kept\nwant %q\ngot  %q", iter, owner, got)
		}
	}
}
