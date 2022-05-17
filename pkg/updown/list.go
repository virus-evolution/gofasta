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

// List gets a list of ATGC SNPs with respect to reference + ambiguous sites for each query sequence in a fasta-format
// alignment, and writes it to file
func List(reference, alignment io.Reader, out io.Writer) error {

	cErr := make(chan error)

	cFR := make(chan fastaio.EncodedFastaRecord, runtime.NumCPU()+50)
	cFRDone := make(chan bool)

	cudLs := make(chan updownLine, runtime.NumCPU()+50)
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
