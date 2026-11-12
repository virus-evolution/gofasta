package updown

import (
	"context"
	"encoding/csv"
	"errors"
	"io"
	"runtime"
	"sort"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
)

// getAmbArr parses the ambiguities field from one line of the output of gofasta
// updown list to an array of 1-based inclusive start-stop pairs of integers,
// which represent tracts of ambiguities
func getAmbArr(s string) ([]int, error) {
	// if there are no ambiguities:
	if len(s) == 0 {
		return make([]int, 0), nil
	}
	ambs := strings.Split(s, "|")
	A := make([]int, 0)
	for _, a := range ambs {
		asplit := strings.Split(a, "-")
		if len(asplit) == 1 {
			a1, err := strconv.Atoi(asplit[0])
			if err != nil {
				return make([]int, 0), errors.New("error parsing ambiguity range from file")
			}
			a2, err := strconv.Atoi(asplit[0])
			if err != nil {
				return make([]int, 0), errors.New("error parsing ambiguity range from file")
			}
			A = append(A, a1)
			A = append(A, a2)
		} else if len(asplit) == 2 {
			a1, err := strconv.Atoi(asplit[0])
			if err != nil {
				return make([]int, 0), errors.New("error parsing ambiguity range from file")
			}
			a2, err := strconv.Atoi(asplit[1])
			if err != nil {
				return make([]int, 0), errors.New("error parsing ambiguity range from file")
			}
			A = append(A, a1)
			A = append(A, a2)
		} else {
			return make([]int, 0), errors.New("error parsing ambiguity range from file")
		}
	}

	return A, nil
}

// headerEqual is a utility function to check that two slices of strings are equal. It is used
// to check the format of the csv-file input to the updown routines
func headerEqual(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i, _ := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// newCSVSource returns a function that reads one record per call from a csv file in the format written
// by gofasta updown list, numbering the records in file order, and returns io.EOF after the last one.
// badHeader is the error to return if the header is not updown list's.
func newCSVSource(in io.Reader, badHeader error) func() (updownLine, error) {
	r := csv.NewReader(in)
	header := true
	counter := 0
	return func() (updownLine, error) {
		for {
			record, err := r.Read()
			if err != nil {
				return updownLine{}, err
			}
			if header {
				if !headerEqual(record, []string{"query", "SNPs", "ambiguities", "SNPcount", "ambcount"}) {
					return updownLine{}, badHeader
				}
				header = false
				continue
			}
			a, err := getAmbArr(record[2])
			if err != nil {
				return updownLine{}, err
			}

			var snps []string
			var snpPos []int
			if len(record[1]) > 0 {
				snps = strings.Split(record[1], "|")
				snpPos = make([]int, len(snps))
				for i, snp := range snps {
					snpPos[i], err = strconv.Atoi(snp[1 : len(snp)-1])
					if err != nil {
						return updownLine{}, err
					}
				}
			} else {
				snps = make([]string, 0)
				snpPos = make([]int, 0)
			}

			amb_count, err := strconv.Atoi(record[4])
			if err != nil {
				return updownLine{}, err
			}

			snpsSorted := make([]string, len(snps))
			copy(snpsSorted, snps)
			sort.Slice(snpsSorted, func(i, j int) bool {
				return snpsSorted[i] < snpsSorted[j]
			})

			udL := updownLine{id: record[0], idx: counter, snps: snps, snpsSorted: snpsSorted, snpsPos: snpPos, ambs: a, ambCount: amb_count}
			counter++
			return udL, nil
		}
	}
}

// readCSVToUDLList reads a csv file in the format produced by gofasta updown
// list to an array of updownLine structs, each of which contains the snps and
// ambiguous positions for one query or target sequence
func readCSVToUDLList(in io.Reader) ([]updownLine, error) {
	next := newCSVSource(in, errors.New("bad header when parsing --query csv: is this file the output of gofasta updown list?"))
	LudL := make([]updownLine, 0)
	for {
		udL, err := next()
		if err == io.EOF {
			return LudL, nil
		}
		if err != nil {
			return make([]updownLine, 0), err
		}
		LudL = append(LudL, udL)
	}
}

// fastaToUDLList converts a fasta format alignment to an array of updownLine
// structs, given a reference sequence, each of which contains the snps and
// ambiguous positions for one query or target sequence
func fastaToUDLList(in io.Reader, refSeq []byte) ([]updownLine, error) {
	DA := encoding.MakeDecodingArray()
	next := fasta.NewEncodedAlignmentSource(in, false, false, false)
	work := func(FR fasta.EncodedRecord) (updownLine, error) { return getLine(refSeq, FR, DA) }
	udla := make([]updownLine, 0)
	err := pipeline.Run(context.Background(), runtime.NumCPU(), next, work, pipeline.FirstError(func(udL updownLine) error {
		udla = append(udla, udL)
		return nil
	}))
	if err != nil {
		return make([]updownLine, 0), err
	}
	return udla, nil
}

// getLine gets the mutation and ambiguity lists between the reference and one fasta record
func getLine(refSeq []byte, FR fasta.EncodedRecord, DA [256]string) (updownLine, error) {

	if len(FR.Seq) != len(refSeq) {
		return updownLine{}, errors.New("alignment and reference are not the same width")
	}

	var (
		snps      = make([]string, 0)
		snpPos    = make([]int, 0)
		ambs      = make([]int, 0) // [1,265,11083,11083,...]: each pair is the 1-based inclusive start/end of a tract of ambiguities
		snpCount  int
		ambCount  int
		amb_start int
		amb_stop  int
		cont      bool // in a tract of ambiguities
	)

	for i, que_nuc := range FR.Seq {

		// if query nucleotide is a known base
		if que_nuc&8 == 8 {
			// if it is different from the reference:
			if (refSeq[i] & que_nuc) < 16 {
				snps = append(snps, DA[refSeq[i]]+strconv.Itoa(i+1)+DA[que_nuc])
				snpPos = append(snpPos, i+1)
				snpCount++
			}
			// also need to finalise any previous ambiguity tract
			if cont {
				ambs = append(ambs, amb_start+1, amb_stop+1)
				cont = false
			}
			// otherwise update the ambiguities
		} else {
			ambCount++
			if cont {
				amb_stop = i
			} else {
				amb_start = i
				amb_stop = i
				cont = true
			}
		}
	}

	// need to finalise any ambiguity tract that reaches the edge of the sequence
	if cont {
		ambs = append(ambs, amb_start+1, amb_stop+1)
	}

	snpsSorted := make([]string, len(snps))
	copy(snpsSorted, snps)
	sort.Slice(snpsSorted, func(i, j int) bool {
		return snpsSorted[i] < snpsSorted[j]
	})

	return updownLine{
		id:         FR.ID,
		idx:        FR.Idx,
		snps:       snps,
		snpCount:   snpCount,
		snpsPos:    snpPos,
		ambs:       ambs,
		ambCount:   ambCount,
		snpsSorted: snpsSorted,
	}, nil
}
