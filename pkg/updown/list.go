package updown

import (
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"runtime"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
)

// updownLine is a struct for one records snps relative to a reference sequence and ambiguity tracts.
// It is meant as a compressed compressed version of a genome to facilitate fast distance calculations
type updownLine struct {
	id         string
	idx        int
	snps       []string
	snpsSorted []string
	snpsPos    []int
	snpCount   int
	ambs       []int //  [1,265,11083,11083,...,...] each pair constitutes the 1-based inclusive start/end positions of a tract of ambiguities
	ambCount   int   // total number of sites that are not ATGC
}

// writeLine writes one record's SNPs and ambiguities as a line of csv
func writeLine(w io.Writer, udLine updownLine) error {
	ambstrings := make([]string, 0)
	for i := 0; i < len(udLine.ambs); i += 2 {
		if udLine.ambs[i] == udLine.ambs[i+1] {
			ambstrings = append(ambstrings, strconv.Itoa(udLine.ambs[i]))
		} else {
			ambstrings = append(ambstrings, strconv.Itoa(udLine.ambs[i])+"-"+strconv.Itoa(udLine.ambs[i+1]))
		}
	}
	_, err := w.Write([]byte(udLine.id + "," + strings.Join(udLine.snps, "|") + "," + strings.Join(ambstrings, "|") + "," + strconv.Itoa(udLine.snpCount) + "," + strconv.Itoa(udLine.ambCount) + "\n"))
	return err
}

// List gets a list of ATGC SNPs with respect to reference + ambiguous sites for each query sequence in a fasta-format
// alignment, and writes it to file
func List(reference, alignment io.Reader, out io.Writer) error {

	temp, err := fasta.LoadEncodeAlignment(reference, false, false, false)
	if err != nil {
		return err
	}
	if len(temp) > 1 {
		return errors.New("More than one record in --reference")
	}

	refSeq := temp[0].Seq

	// every site in the reference would ideally be ∈ {A,T,G,C}
	DA := encoding.MakeDecodingArray()
	for _, nuc := range refSeq {
		if !(nuc&8 == 8) {
			fmt.Fprintf(os.Stderr, "Warning: there is at least one ambiguous nucleotide in the --reference sequence: %s. This isn't recommended\n", DA[nuc])
			break
		}
	}

	if _, err := out.Write([]byte("query,SNPs,ambiguities,SNPcount,ambcount\n")); err != nil {
		return err
	}
	next := fasta.NewEncodedAlignmentSource(alignment, false, false, false)
	work := func(FR fasta.EncodedRecord) (updownLine, error) { return getLine(refSeq, FR, DA) }
	return pipeline.Run(context.Background(), runtime.NumCPU(), next, work, pipeline.FirstError(func(udL updownLine) error {
		return writeLine(out, udL)
	}))
}
