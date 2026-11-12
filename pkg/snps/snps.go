/*
Package snps implements functions to call nucleotide changes between each sequence
in a fasta format alignment and a reference sequence.
*/
package snps

import (
	"context"
	"errors"
	"io"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
)

// snpLine is a struct for one fasta record's SNPs
type snpLine struct {
	queryname string
	snps      []string
}

// getSNPs gets the SNPs between the reference sequence and one fasta record
func getSNPs(refSeq []byte, FR fasta.EncodedRecord, DA [256]string) (snpLine, error) {
	if len(FR.Seq) != len(refSeq) {
		rl := strconv.Itoa(len(refSeq))
		ql := strconv.Itoa(len(FR.Seq))
		return snpLine{}, errors.New("Reference sequence (" + rl + " bases) and " + FR.ID + " (" + ql + " bases) are different lengths")
	}
	SNPs := make([]string, 0)
	for i, nuc := range FR.Seq {
		// identical codes are never a difference (this matters for --hard-gaps, where a
		// gap is encoded as 4 and 4&4 < 16)
		if refSeq[i] != nuc && (refSeq[i]&nuc) < 16 {
			SNPs = append(SNPs, DA[refSeq[i]]+strconv.Itoa(i+1)+DA[nuc])
		}
	}
	return snpLine{queryname: FR.ID, snps: SNPs}, nil
}

// writeLine writes one record's SNPs
func writeLine(w io.Writer, SL snpLine) error {
	_, err := w.Write([]byte(SL.queryname + "," + strings.Join(SL.snps, "|") + "\n"))
	return err
}

// writeAggregate writes the frequency of each SNP among n records, for those at or above threshold,
// sorted by position and then by alternative base
func writeAggregate(w io.Writer, propMap map[string]float64, n float64, threshold float64) error {
	order := make([]string, 0, len(propMap))
	for k := range propMap {
		order = append(order, k)
	}
	// a SNP is reference base, position, alternative base; the bases are single characters
	position := func(snp string) int {
		pos, _ := strconv.Atoi(snp[1 : len(snp)-1])
		return pos
	}
	sort.SliceStable(order, func(i, j int) bool {
		pos_i, pos_j := position(order[i]), position(order[j])
		alt_i := order[i][len(order[i])-1]
		alt_j := order[j][len(order[j])-1]
		return pos_i < pos_j || (pos_i == pos_j && alt_i < alt_j)
	})
	for _, snp := range order {
		if propMap[snp]/n < threshold {
			continue
		}
		_, err := w.Write([]byte(snp + "," + strconv.FormatFloat(propMap[snp]/n, 'f', 9, 64) + "\n"))
		if err != nil {
			return err
		}
	}
	return nil
}

// SNPs annotates snps for each record in a fasta-format alignment with respect to a reference sequence.
// A record that can't be processed (an invalid character, or a different length from the reference) is
// reported on stderr and skipped, unless strict, when it stops the run with an error.
func SNPs(ref, alignment io.Reader, hardGaps bool, aggregate bool, threshold float64, w io.Writer, strict bool) error {

	refs, err := fasta.LoadEncodeAlignment(ref, hardGaps, false, false)
	if err != nil {
		return err
	}
	if len(refs) > 1 {
		return errors.New("more than one record in --reference")
	}
	refSeq := refs[0].Seq

	DA := encoding.MakeDecodingArray()
	next := fasta.NewRecordSource(alignment)
	work := func(record fasta.Record) (snpLine, error) {
		FR, err := fasta.EncodeRecord(record, hardGaps, false, false)
		if err != nil {
			return snpLine{}, pipeline.Named(record.ID, err)
		}
		SL, err := getSNPs(refSeq, FR, DA)
		return SL, pipeline.Named(record.ID, err)
	}

	var write func(snpLine) error
	propMap := make(map[string]float64)
	n := 0.0
	if aggregate {
		if _, err := w.Write([]byte("SNP,frequency\n")); err != nil {
			return err
		}
		write = func(SL snpLine) error {
			n++
			for _, snp := range SL.snps {
				propMap[snp]++
			}
			return nil
		}
	} else {
		if _, err := w.Write([]byte("query,SNPs\n")); err != nil {
			return err
		}
		write = func(SL snpLine) error { return writeLine(w, SL) }
	}

	h := pipeline.NewHandler(strict, os.Stderr, pipeline.Noun{One: "record", Many: "records"}, write)
	if err := pipeline.Run(context.Background(), runtime.NumCPU(), next, work, h.Sink); err != nil {
		return err
	}
	if err := h.Finish(); err != nil {
		return err
	}
	if aggregate {
		return writeAggregate(w, propMap, n, threshold)
	}
	return nil
}
