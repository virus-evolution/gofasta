package sam

// Regression tests.

import (
	"math/rand"
	"strconv"
	"strings"
	"testing"

	biogosam "github.com/biogo/hts/sam"
)

// pairAlignFromSam runs one query's SAM records through the same steps as toPairAlign (grouping the records
// of the query, which resolves overlaps between them, then blockToPairwiseAlignment with insertions retained)
// and returns the resulting pairwise alignment.
func pairAlignFromSam(t *testing.T, samText string, ref string) alignPair {
	t.Helper()
	return pairAlignFromSamWith(t, samText, ref, true)
}

// pairAlignFromSamUnresolved is as pairAlignFromSam but without resolving overlaps between the records, so
// that blockToSeqPair is given overlapping rows. resolveOverlaps normally prevents that, but blockToSeqPair
// must be correct for any block.
func pairAlignFromSamUnresolved(t *testing.T, samText string, ref string) alignPair {
	t.Helper()
	return pairAlignFromSamWith(t, samText, ref, false)
}

func pairAlignFromSamWith(t *testing.T, samText string, ref string, resolve bool) alignPair {
	t.Helper()

	r, err := biogosam.NewReader(strings.NewReader(samText))
	if err != nil {
		t.Fatal(err)
	}
	group := samRecords{}
	for {
		rec, err := r.Read()
		if err != nil {
			break
		}
		group.records = append(group.records, *rec)
	}
	if resolve {
		group.records = resolveOverlaps(group.records)
	}

	pair, err := blockToPairwiseAlignment(group, []byte(ref), false)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	return pair
}

// when a query has several SAM records (e.g. primary + supplementary) and one
// contains an insertion, the bases after the insertion point are kept in every row.
func TestBlockToSeqPairSupplementaryAfterInsertion(t *testing.T) {
	ref := "ACGTACGTACGTACGTACGT"

	// primary: 5M3I15M. supplementary (flag 2048): 20M, which disagrees with the primary
	// at two of the columns they share after the insertion, so those must come out as N.
	samText := "@SQ\tSN:ref\tLN:20\n" +
		"q\t0\tref\t1\t60\t5M3I15M\t*\t0\t0\tACGTAGGGCGTACGTACGTACGT\t*\n" +
		"q\t2048\tref\t1\t60\t20M\t*\t0\t0\tACGTATTTACGTACGTACGT\t*\n"

	pair := pairAlignFromSamUnresolved(t, samText, ref)

	wantRef := "ACGTA---CGTACGTACGTACGT"
	wantQue := "ACGTAGGGNNTACGTACGTACGT"

	if string(pair.ref) != wantRef {
		t.Errorf("ref:   want %s got %s", wantRef, pair.ref)
	}
	if string(pair.query) != wantQue {
		t.Errorf("query: want %s got %s", wantQue, pair.query)
	}
}

// a second insertion in the same record keeps the gaps for the first insertion in the other rows.
func TestBlockToSeqPairTwoInsertionsInOneRecord(t *testing.T) {
	ref := "ACGTACGTACGTACGTACGT"

	// primary: 5M3I5M3I10M. supplementary: 20M identical to the reference up to base 10,
	// then the complement of the reference for bases 11-20 (so all of those columns conflict).
	samText := "@SQ\tSN:ref\tLN:20\n" +
		"q\t0\tref\t1\t60\t5M3I5M3I10M\t*\t0\t0\tACGTAGGGCGTACTTTGTACGTACGT\t*\n" +
		"q\t2048\tref\t1\t60\t20M\t*\t0\t0\tACGTACGTACCATGCATGCA\t*\n"

	pair := pairAlignFromSamUnresolved(t, samText, ref)

	wantRef := "ACGTA---CGTAC---GTACGTACGT"
	wantQue := "ACGTAGGGCGTACTTTNNNNNNNNNN"

	if string(pair.ref) != wantRef {
		t.Errorf("ref:   want %s got %s", wantRef, pair.ref)
	}
	if string(pair.query) != wantQue {
		t.Errorf("query: want %s got %s", wantQue, pair.query)
	}
}

// a row's own insertions are counted when placing other rows' insertions. Here the
// primary has an insertion at reference base 5 and the supplementary has one at base 12; the
// primary's gaps for the supplementary's insertion must go after the primary's own insertion.
func TestBlockToSeqPairInsertionsInDifferentRecords(t *testing.T) {
	ref := "ACGTACGTACGTACGTACGT"

	// supplementary: 12M2I8M with CC inserted after reference base 12
	samText := "@SQ\tSN:ref\tLN:20\n" +
		"q\t0\tref\t1\t60\t5M3I15M\t*\t0\t0\tACGTAGGGCGTACGTACGTACGT\t*\n" +
		"q\t2048\tref\t1\t60\t12M2I8M\t*\t0\t0\tACGTACGTACGTCCACGTACGT\t*\n"

	pair := pairAlignFromSamUnresolved(t, samText, ref)

	wantRef := "ACGTA---CGTACGT--ACGTACGT"
	wantQue := "ACGTAGGGCGTACGTCCACGTACGT"

	if string(pair.ref) != wantRef {
		t.Errorf("ref:   want %s got %s", wantRef, pair.ref)
	}
	if string(pair.query) != wantQue {
		t.Errorf("query: want %s got %s", wantQue, pair.query)
	}
}

// insertions from two records at the SAME reference position are laid out one after the
// other (in record order), so every row gets the same columns.
func TestBlockToSeqPairInsertionsAtSamePosition(t *testing.T) {
	ref := "ACGTACGTACGTACGTACGT"

	samText := "@SQ\tSN:ref\tLN:20\n" +
		"q\t0\tref\t1\t60\t5M3I15M\t*\t0\t0\tACGTAGGGCGTACGTACGTACGT\t*\n" +
		"q\t2048\tref\t1\t60\t5M2I15M\t*\t0\t0\tACGTACCCGTACGTACGTACGT\t*\n"

	pair := pairAlignFromSamUnresolved(t, samText, ref)

	wantRef := "ACGTA-----CGTACGTACGTACGT"
	wantQue := "ACGTAGGGCCCGTACGTACGTACGT"

	if string(pair.ref) != wantRef {
		t.Errorf("ref:   want %s got %s", wantRef, pair.ref)
	}
	if string(pair.query) != wantQue {
		t.Errorf("query: want %s got %s", wantQue, pair.query)
	}
}

// rows that end before, or start after, another row's insertion.
// The primary covers only reference bases 1-10 and the supplementary only 11-20 (with an
// insertion after base 15), so neither row can be given gaps for an insertion it does not span.
func TestBlockToSeqPairRowsEndingEarlyOrStartingLate(t *testing.T) {
	ref := "ACGTACGTACGTACGTACGT"

	samText := "@SQ\tSN:ref\tLN:20\n" +
		"q\t0\tref\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n" +
		"q\t2048\tref\t11\t60\t5M2I5M\t*\t0\t0\tGTACGCCTACGT\t*\n"

	pair := pairAlignFromSam(t, samText, ref)

	wantRef := "ACGTACGTACGTACG--TACGT"
	wantQue := "ACGTACGTACGTACGCCTACGT"

	if string(pair.ref) != wantRef {
		t.Errorf("ref:   want %s got %s", wantRef, pair.ref)
	}
	if string(pair.query) != wantQue {
		t.Errorf("query: want %s got %s", wantQue, pair.query)
	}
}

// Randomised consistency check for blockToSeqPair with several records per query, each with
// random insertions and deletions: whatever the records look like, the aligned reference must
// degap to exactly the reference, and the reference and query rows must be the same length.
// (Bases in M operations are copied from the reference and inserted bases are never shared
// between records, so there are no conflicting bases.)
func TestBlockToSeqPairRandomisedConsistency(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	refBases := "ACGT"

	for iter := 0; iter < 300; iter++ {
		refLen := 20 + rng.Intn(20)
		refB := make([]byte, refLen)
		for i := range refB {
			refB[i] = refBases[rng.Intn(4)]
		}
		ref := string(refB)

		var sb strings.Builder
		sb.WriteString("@SQ\tSN:ref\tLN:" + strconv.Itoa(refLen) + "\n")

		nRecords := 1 + rng.Intn(3)
		for r := 0; r < nRecords; r++ {
			pos := rng.Intn(refLen / 2) // 0-based
			rpos := pos
			var cigar, seq strings.Builder
			for rpos < refLen {
				l := 1 + rng.Intn(5)
				switch op := rng.Intn(4); {
				case op == 0: // insertion (not at the very start or end of the record)
					if cigar.Len() == 0 {
						continue
					}
					cigar.WriteString(strconv.Itoa(l) + "I")
					for k := 0; k < l; k++ {
						seq.WriteByte("CG"[rng.Intn(2)])
					}
				case op == 1: // deletion
					if cigar.Len() == 0 || rpos+l >= refLen {
						continue
					}
					cigar.WriteString(strconv.Itoa(l) + "D")
					rpos += l
				default: // match
					if rpos+l > refLen {
						l = refLen - rpos
					}
					cigar.WriteString(strconv.Itoa(l) + "M")
					seq.WriteString(ref[rpos : rpos+l])
					rpos += l
				}
				// sometimes stop early, on a match
				if rng.Intn(6) == 0 && strings.HasSuffix(cigar.String(), "M") {
					break
				}
			}
			flag := 0
			if r > 0 {
				flag = 2048
			}
			sb.WriteString("q\t" + strconv.Itoa(flag) + "\tref\t" + strconv.Itoa(pos+1) + "\t60\t" + cigar.String() + "\t*\t0\t0\t" + seq.String() + "\t*\n")
		}

		for _, resolve := range []bool{false, true} {
			pair := pairAlignFromSamWith(t, sb.String(), ref, resolve)

			if got := strings.ReplaceAll(string(pair.ref), "-", ""); got != ref {
				t.Fatalf("iteration %d (resolve=%v): degapped reference row is wrong\nsam:\n%swant %s\ngot  %s\nrow  %s", iter, resolve, sb.String(), ref, got, pair.ref)
			}
			if len(pair.ref) != len(pair.query) {
				t.Fatalf("iteration %d (resolve=%v): reference (%d) and query (%d) rows differ in length\nsam:\n%s", iter, resolve, len(pair.ref), len(pair.query), sb.String())
			}
		}
	}
}
