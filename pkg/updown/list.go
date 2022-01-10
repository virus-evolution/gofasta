package updown

import (
	"errors"
	"fmt"
	"io"
	"os"
	"runtime"
	"strconv"
	"strings"
	"sync"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fastaio"
)

type updownLine struct {
	id       string
	idx      int
	snps     []string
	snpsPos  []int
	snpCount int
	ambs     []int //  [1,265,11083,11083,...,...] each pair constitutes the 1-based inclusive start/end positions of a tract of ambiguities
	ambCount int   // total number of sites that are not ATGC
}

// getLine gets the mutation + ambiguity lists between the reference and each Fasta record at a time
func getLines(refSeq []byte, cFR chan fastaio.EncodedFastaRecord, cUDs chan updownLine, cErr chan error) {

	DA := encoding.MakeDecodingArray()

	var snp string
	var snps []string
	var snpCount int
	var snpPos []int
	var ambs []int // [1,265,11083,11083,...,...] len(ambs) %% 2 must equal 0: each pair constitutes the 1-based inclusive start/end positions of a tract of ambiguities
	var ambCount int
	var udLine updownLine
	var amb_start int
	var amb_stop int
	var cont bool

	for FR := range cFR {

		if len(FR.Seq) != len(refSeq) {
			cErr <- errors.New("alignment and reference are not the same width")
		}

		cont = false // for tracts of ambiguities

		udLine = updownLine{}
		udLine.id = FR.ID
		udLine.idx = FR.Idx
		snps = make([]string, 0)
		snpPos = make([]int, 0)
		ambs = make([]int, 0)
		snpCount = 0
		ambCount = 0

		for i, que_nuc := range FR.Seq {

			// if query nucleotide is a known base
			if que_nuc&8 == 8 {
				// if it is different from the reference:
				if (refSeq[i] & que_nuc) < 16 {
					snp = DA[refSeq[i]] + strconv.Itoa(i+1) + DA[que_nuc]
					snps = append(snps, snp)
					snpPos = append(snpPos, i+1)
					snpCount++
				}
				// also need to finalise any previous ambiguity tract
				if cont {
					ambs = append(ambs, amb_start+1)
					ambs = append(ambs, amb_stop+1)
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
			ambs = append(ambs, amb_start+1)
			ambs = append(ambs, amb_stop+1)
		}

		udLine.snps = snps
		udLine.snpCount = snpCount
		udLine.snpsPos = snpPos
		udLine.ambs = ambs
		udLine.ambCount = ambCount
		cUDs <- udLine
	}

	return
}

// writeOutput writes the output to stdout or a file as it arrives.
// It uses a map to write things in the same order as they are in the input file.
func writeOutput(w io.Writer, cudLs chan updownLine, cErr chan error, cWriteDone chan bool) {

	outputMap := make(map[int]updownLine)

	counter := 0

	var err error

	_, err = w.Write([]byte("query,SNPs,ambiguities,SNPcount,ambcount\n"))
	if err != nil {
		cErr <- err
		return
	}

	var ambstrings []string

	for udL := range cudLs {
		outputMap[udL.idx] = udL

		if udLine, ok := outputMap[counter]; ok {
			ambstrings = make([]string, 0)
			for i := 0; i < len(udLine.ambs); i += 2 {
				if udLine.ambs[i] == udLine.ambs[i+1] {
					ambstrings = append(ambstrings, strconv.Itoa(udLine.ambs[i]))
				} else {
					ambstrings = append(ambstrings, strconv.Itoa(udLine.ambs[i])+"-"+strconv.Itoa(udLine.ambs[i+1]))
				}
			}
			_, err := w.Write([]byte(udLine.id + "," + strings.Join(udLine.snps, "|") + "," + strings.Join(ambstrings, "|") + "," + strconv.Itoa(udLine.snpCount) + "," + strconv.Itoa(udLine.ambCount) + "\n"))
			if err != nil {
				cErr <- err
				return
			}
			delete(outputMap, counter)
			counter++
		} else {
			continue
		}
	}

	for n := 1; n > 0; {
		if len(outputMap) == 0 {
			n--
			break
		}
		udLine := outputMap[counter]
		ambstrings = make([]string, 0)
		for i := 0; i < len(udLine.ambs); i += 2 {
			if udLine.ambs[i] == udLine.ambs[i+1] {
				ambstrings = append(ambstrings, strconv.Itoa(udLine.ambs[i]))
			} else {
				ambstrings = append(ambstrings, strconv.Itoa(udLine.ambs[i])+"-"+strconv.Itoa(udLine.ambs[i+1]))
			}
		}
		_, err := w.Write([]byte(udLine.id + "," + strings.Join(udLine.snps, "|") + "," + strings.Join(ambstrings, "|") + "," + strconv.Itoa(udLine.snpCount) + "," + strconv.Itoa(udLine.ambCount) + "\n"))
		if err != nil {
			cErr <- err
			return
		}
		delete(outputMap, counter)
		counter++
	}

	cWriteDone <- true
}

// List gets a list of ATGC SNPs and ambiguous sites for each query
func List(reference, alignment io.Reader, out io.Writer) error {

	cErr := make(chan error)

	cFR := make(chan fastaio.EncodedFastaRecord)
	cFRDone := make(chan bool)

	cudLs := make(chan updownLine, runtime.NumCPU())
	cudLsDone := make(chan bool)

	cWriteDone := make(chan bool)

	temp, err := fastaio.ReadEncodeAlignmentToList(reference, false)
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

	go fastaio.ReadEncodeAlignment(alignment, false, cFR, cErr, cFRDone)

	go writeOutput(out, cudLs, cErr, cWriteDone)

	var wgudLs sync.WaitGroup
	wgudLs.Add(runtime.NumCPU())

	for n := 0; n < runtime.NumCPU(); n++ {
		go func() {
			getLines(refSeq, cFR, cudLs, cErr)
			wgudLs.Done()
		}()
	}

	go func() {
		wgudLs.Wait()
		cudLsDone <- true
	}()

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cFRDone:
			close(cFR)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cudLsDone:
			close(cudLs)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cWriteDone:
			n--
		}
	}

	return nil
}
