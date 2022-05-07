package updown

import (
	"encoding/csv"
	"errors"
	"io"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fastaio"
)

// getAmbArr parses the ambiguities field from one line of the output of gofasta updown list to an array of
// 1-based inclusive start-stop pairs of integers which represent tracts of ambiguities
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
// to check the format of the csv-file input to updown routines
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

// readCSVToUDLChan reads a csv file in the format produced by gofasta updown list to a channel of
// updownLine structs, each of which contains the snps and ambiguous positions for one query or target sequence
func readCSVToUDLChan(in io.Reader, cudL chan updownLine, cErr chan error, cReadDone chan bool) {

	var snps []string
	var snpPos []int

	header := true
	r := csv.NewReader(in)

	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			cErr <- err
			return
		}
		if header {
			if !headerEqual(record, []string{"query", "SNPs", "ambiguities", "SNPcount", "ambcount"}) {
				cErr <- errors.New("bad header when parsing --target csv: is this the output of gofasta updown list?")
				return
			}
			header = false
			continue
		}
		a, err := getAmbArr(record[2])
		if err != nil {
			cErr <- err
			return
		}

		if len(record[1]) > 0 {
			snps = strings.Split(record[1], "|")
			snpPos = make([]int, len(snps))
			for i, snp := range snps {
				snpPos[i], err = strconv.Atoi(snp[1 : len(snp)-1])
				if err != nil {
					cErr <- err
					return
				}
			}
		} else {
			snps = make([]string, 0)
			snpPos = make([]int, 0)
		}

		amb_count, err := strconv.Atoi(record[4])
		if err != nil {
			cErr <- err
			return
		}

		snpsSorted := make([]string, len(snps))
		copy(snpsSorted, snps)
		sort.Slice(snpsSorted, func(i, j int) bool {
			return snpsSorted[i] < snpsSorted[j]
		})

		udL := updownLine{id: record[0], snps: snps, snpsSorted: snpsSorted, snpsPos: snpPos, ambs: a, ambCount: amb_count}
		cudL <- udL
	}

	cReadDone <- true
}

// readCSVToUDLList reads a csv file in the format produced by gofasta updown list to an array of
// updownLine structs, each of which contains the snps and ambiguous positions for one query or target sequence
func readCSVToUDLList(in io.Reader) ([]updownLine, error) {

	LudL := make([]updownLine, 0)

	var snps []string
	var snpPos []int

	header := true
	r := csv.NewReader(in)

	counter := 0
	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return make([]updownLine, 0), err
		}
		if header {
			if !headerEqual(record, []string{"query", "SNPs", "ambiguities", "SNPcount", "ambcount"}) {
				return make([]updownLine, 0), errors.New("bad header when parsing --query csv: is this file the output of gofasta updown list?")
			}
			header = false
			continue
		}
		a, err := getAmbArr(record[2])
		if err != nil {
			return make([]updownLine, 0), err
		}

		if len(record[1]) > 0 {
			snps = strings.Split(record[1], "|")
			snpPos = make([]int, len(snps))
			for i, snp := range snps {
				snpPos[i], err = strconv.Atoi(snp[1 : len(snp)-1])
				if err != nil {
					return make([]updownLine, 0), err
				}
			}
		} else {
			snps = make([]string, 0)
			snpPos = make([]int, 0)
		}

		amb_count, err := strconv.Atoi(record[4])
		if err != nil {
			return make([]updownLine, 0), err
		}

		snpsSorted := make([]string, len(snps))
		copy(snpsSorted, snps)
		sort.Slice(snpsSorted, func(i, j int) bool {
			return snpsSorted[i] < snpsSorted[j]
		})

		udL := updownLine{id: record[0], snps: snps, snpsSorted: snpsSorted, snpsPos: snpPos, ambs: a, ambCount: amb_count}

		LudL = append(LudL, udL)
		counter++
	}

	return LudL, nil
}

// fastaToUDLslice converts a fasta format alignment to an array of updownLine structs, given a reference sequence,
// each of which contains the snps and ambiguous positions for one query or target sequence
func fastaToUDLList(in io.Reader, refSeq []byte) ([]updownLine, error) {

	var udla []updownLine

	cInternalErr := make(chan error)

	cFR := make(chan fastaio.EncodedFastaRecord)
	cFRDone := make(chan bool)
	cudLs := make(chan updownLine, runtime.NumCPU())
	cudLsDone := make(chan bool)
	cArrayDone := make(chan bool)

	go fastaio.ReadEncodeAlignment(in, false, cFR, cInternalErr, cFRDone)

	var wgudLs sync.WaitGroup
	wgudLs.Add(1)

	go func() {
		getLines(refSeq, cFR, cudLs, cInternalErr)
		wgudLs.Done()
	}()

	go func() {
		wgudLs.Wait()
		cudLsDone <- true
	}()

	go func() {
		for udl := range cudLs {
			udla = append(udla, udl)
		}
		cArrayDone <- true
	}()

	for n := 1; n > 0; {
		select {
		case err := <-cInternalErr:
			return make([]updownLine, 0), err
		case <-cFRDone:
			close(cFR)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cInternalErr:
			return make([]updownLine, 0), err
		case <-cudLsDone:
			close(cudLs)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case <-cArrayDone:
			n--
		}
	}

	return udla, nil
}

// readFastaToUDLChan converts each record in a fasta format alignment to an updownLine struct, given a reference sequence,
// and passes it to a channel
func readFastaToUDLChan(target io.Reader, refSeq []byte, cudL chan updownLine, cErr chan error, cReadDone chan bool) {
	cInternalErr := make(chan error)

	cFR := make(chan fastaio.EncodedFastaRecord)
	cFRDone := make(chan bool)

	cReOrder := make(chan updownLine)
	cReOrderDone := make(chan bool)

	cudLsDone := make(chan bool)

	go fastaio.ReadEncodeAlignment(target, false, cFR, cInternalErr, cFRDone)

	var wgudLs sync.WaitGroup
	wgudLs.Add(runtime.NumCPU())

	for n := 0; n < runtime.NumCPU(); n++ {
		go func() {
			getLines(refSeq, cFR, cReOrder, cInternalErr)
			wgudLs.Done()
		}()
	}

	go reorderRecords(cReOrder, cudL, cReOrderDone)

	go func() {
		wgudLs.Wait()
		cudLsDone <- true
	}()

	for n := 1; n > 0; {
		select {
		case err := <-cInternalErr:
			cErr <- err
			return
		case <-cFRDone:
			close(cFR)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cInternalErr:
			cErr <- err
			return
		case <-cudLsDone:
			close(cReOrder)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cInternalErr:
			cErr <- err
			return
		case <-cReOrderDone:
			cReadDone <- true
			n--
		}
	}
}

// reorderRecords reorders the records in a channel of updownLine structs according to the order they were
// in the input. It does this to ensure that given the same dataset, the results of the updown routines will
// be the same whether the input target file is in csv or fasta format. This is necessary because when there are
// ties for genetic distance and ambiguity content, sort.SliceStable will preserve input order as precedence.
func reorderRecords(cIn, cOut chan updownLine, cReorderDone chan bool) {

	reorderMap := make(map[int]updownLine)

	counter := 0

	for input := range cIn {
		reorderMap[input.idx] = input
		if output, ok := reorderMap[counter]; ok {

			cOut <- output

			delete(reorderMap, counter)
			counter++
		} else {
			continue
		}
	}

	for n := 1; n > 0; {
		if len(reorderMap) == 0 {
			n--
			break
		}
		output := reorderMap[counter]

		cOut <- output

		delete(reorderMap, counter)
		counter++
	}

	cReorderDone <- true
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

		snpsSorted := make([]string, len(snps))
		copy(snpsSorted, snps)
		sort.Slice(snpsSorted, func(i, j int) bool {
			return snpsSorted[i] < snpsSorted[j]
		})

		udLine.snps = snps
		udLine.snpCount = snpCount
		udLine.snpsPos = snpPos
		udLine.ambs = ambs
		udLine.ambCount = ambCount
		udLine.snpsSorted = snpsSorted

		cUDs <- udLine
	}

	return
}
