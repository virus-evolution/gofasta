package updown

import (
	"io"
	"os"
	"sync"
	"errors"
	"strings"
	"strconv"
	"runtime"
	"encoding/csv"

	"github.com/cov-ert/gofasta/pkg/fastaio"
)

func getAmbArr(s string) ([]int, error) {
	// if there are no ambiguities:
	if len(s) == 0 {
		return make([]int, 0), nil
	}
	ambs := strings.Split(s, "|")
	A := make([]int, 0)
	for _, a := range(ambs) {
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

func headerEqual(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i, _ := range(a) {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func readCSVToChan(inFile string, cudL chan updownLine, cErr chan error, cReadDone chan bool) {
	f, err := os.Open(inFile)
	if err != nil {
		cErr<- err
	}
	defer f.Close()

	var snps []string
	var snpPos []int

	header := true
	r := csv.NewReader(f)

	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			cErr<- err
		}
		if header {
			if ! headerEqual(record, []string{"query","SNPs","ambiguities","SNPcount","ambcount"}) {
				cErr<- errors.New("bad header when parsing --target csv: is this the output of gofasta updown list?")
			}
			header = false
			continue
		}
		a, err := getAmbArr(record[2])
		if err != nil {
			cErr<- err
		}

		if len(record[1]) > 0 {
			snps = strings.Split(record[1], "|")
			snpPos = make([]int, len(snps))
			for i, snp := range(snps) {
				snpPos[i], err = strconv.Atoi(snp[1:len(snp) - 1])
				if err != nil {
					cErr<- err
				}
			}
		} else {
			snps = make([]string, 0)
			snpPos = make([]int, 0)
		}

		amb_count, err := strconv.Atoi(record[4])
		if err != nil {
			cErr<- err
		}

		udL := updownLine{id: record[0], snps: snps, snpsPos: snpPos, ambs: a, ambCount: amb_count}
		cudL<- udL
	}

	cReadDone<- true
}

func readCSVToList(inFile string) ([]updownLine, error) {
	f, err := os.Open(inFile)
	if err != nil {
		return make([]updownLine, 0), err
	}
	defer f.Close()

	LudL := make([]updownLine, 0)

	var snps []string
	var snpPos []int

	header := true
	r := csv.NewReader(f)

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
			if ! headerEqual(record,[]string{"query","SNPs","ambiguities","SNPcount","ambcount"}) {
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
			for i, snp := range(snps) {
				snpPos[i], err = strconv.Atoi(snp[1:len(snp) - 1])
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
		udL := updownLine{id: record[0], idx: counter, snps: snps, snpsPos: snpPos, ambs: a, ambCount: amb_count}
		LudL = append(LudL, udL)
		counter++
	}

	return LudL, nil
}

func fastaToUDLslice(infile string, refSeq []byte) ([]updownLine, error) {

	var udla []updownLine

	cInternalErr := make(chan error)

	cFR := make(chan fastaio.EncodedFastaRecord)
	cFRDone := make(chan bool)
	cudLs := make(chan updownLine, runtime.NumCPU())
	cudLsDone := make(chan bool)
	cArrayDone := make(chan bool)

	go fastaio.ReadEncodeAlignment(infile, cFR, cInternalErr, cFRDone)

	var wgudLs sync.WaitGroup
	wgudLs.Add(1)

	go func() {
		getLines(refSeq, cFR, cudLs, cInternalErr)
		wgudLs.Done()
	}()

	go func() {
		wgudLs.Wait()
		cudLsDone<- true
	}()

	go func() {
		for udl := range(cudLs) {
			udla = append(udla, udl)
		}
		cArrayDone<- true
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

func reorderRecords(cIn, cOut chan updownLine, cReorderDone chan bool) {

	reorderMap := make(map[int]updownLine)

	counter := 0

	for input := range(cIn) {
		reorderMap[input.idx] = input
		if output, ok := reorderMap[counter]; ok {

			cOut<- output

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

		cOut<- output

		delete(reorderMap, counter)
		counter++
	}

	cReorderDone<- true
}

func readFastaToChan(target string, refSeq []byte, cudL chan updownLine, cErr chan error, cReadDone chan bool) {
	cInternalErr := make(chan error)

	cFR := make(chan fastaio.EncodedFastaRecord)
	cFRDone := make(chan bool)

	cReOrder := make(chan updownLine)
	cReOrderDone := make(chan bool)

	cudLsDone := make(chan bool)

	go fastaio.ReadEncodeAlignment(target, cFR, cInternalErr, cFRDone)

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
		cudLsDone<- true
	}()

	for n := 1; n > 0; {
		select {
		case err := <-cInternalErr:
			cErr<- err
		case <-cFRDone:
			close(cFR)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cInternalErr:
			cErr<- err
		case <-cudLsDone:
			close(cReOrder)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cInternalErr:
			cErr<- err
		case <-cReOrderDone:
			cReadDone<- true
			n--
		}
	}
}
