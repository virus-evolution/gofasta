package sam

import (
	"errors"
	"fmt"
	"os"
	"sync"

	"github.com/cov-ert/gofasta/pkg/fastaio"

	// biogosam "github.com/biogo/hts/sam"
)

// getFastaRecord returns a FastaRecord struct with a sequence ID and a sequence
// that has been optionally trimmed and padded
func getFastaRecord(rawseq []byte, id string, idx int, trim bool, pad bool, trimstart int,
	trimend int) fastaio.FastaRecord {

	var seq []byte

	if pad {
		seq = swapInNs(rawseq)
	} else {
		seq = swapInGapsNs(rawseq)
	}

	if trim {
		if pad {
			for i, _ := range seq {
				if i < trimstart || i >= trimend {
					seq[i] = 'N'
				}
			}
		} else {
			seq = seq[trimstart:trimend]
		}
	}

	FR := fastaio.FastaRecord{ID: id, Description: id, Seq: string(seq), Idx: idx}

	return FR
}

// sanity checks the trimming and padding arguments (given the length of the ref seq)
func checkArgs(refLen int, trim bool, pad bool, trimstart int, trimend int) error {

	if trim {
		if trimstart > refLen-2 || trimstart < 1 {
			return errors.New("error parsing trimming coordinates: check or include --trimstart")
		}
		if trimend > refLen-1 || trimend < 1 {
			return errors.New("error parsing trimming coordinates: check or include --trimend")
		}
		if trimstart >= trimend {
			return errors.New("error parsing trimming coordinates: check trimstart and trimend")
		}
	}

	if pad && !trim {
		_, err := fmt.Fprintln(os.Stderr, "warning: using --pad without --trim has no effect")
		if err != nil {
			return err
		}
	}

	return nil
}

// worker function that takes items from a channel of sam block structs (with indices)
// and writes the corresponding fasta records to a channel
func blockToFastaRecord(ch_in chan samRecords, ch_out chan fastaio.FastaRecord, ch_err chan error,
	refLen int, trim bool, pad bool, trimstart int, trimend int, includeInsertions bool) {

	for group := range ch_in {

		id := group.records[0].Name
		rawseq, err := getSeqFromBlock(group.records, refLen, includeInsertions)
		if err != nil {
			ch_err <- err
		}
		ch_out <- getFastaRecord(rawseq, id, group.idx, trim, pad, trimstart, trimend)
	}
	return
}

// writeAlignmentOut reads fasta records from a channel and writes them to a single
// outfile, in the order in which they are present in the input file.
// It passes a true to a done channel when the channel of fasta records is empty
func writeAlignmentOut(ch chan fastaio.FastaRecord, outfile string, cdone chan bool, cerr chan error) {

	outputMap := make(map[int]fastaio.FastaRecord)

	counter := 0

	var f *os.File
	var err error

	if outfile != "stdout" {
		f, err = os.Create(outfile)
		if err != nil {
			cerr <- err
		}
	} else {
		f = os.Stdout
	}

	defer f.Close()

	for FR := range ch {

		outputMap[FR.Idx] = FR

		if fastarecord, ok := outputMap[counter]; ok {
			_, err = f.WriteString(">" + fastarecord.ID + "\n")
			if err != nil {
				cerr <- err
			}
			_, err = f.WriteString(fastarecord.Seq + "\n")
			if err != nil {
				cerr <- err
			}
			delete(outputMap, counter)
			counter++
		} else {
			continue
		}
	}

	for n := 1; n > 0; {
		if len(outputMap) == 0 {
			break
		}
		fastarecord := outputMap[counter]
		_, err = f.WriteString(">" + fastarecord.ID + "\n")
		if err != nil {
			cerr <- err
		}
		_, err = f.WriteString(fastarecord.Seq + "\n")
		if err != nil {
			cerr <- err
		}
		delete(outputMap, counter)
		counter++
	}

	cdone <- true
}

// ToMultiAlign converts a SAM file to a fasta-format alignment
// Insertions relative to the reference are discarded.
func ToMultiAlign(infile string, reffile string, outfile string, trim bool, pad bool, trimstart int,
	trimend int, threads int) error {

	// samHeader, err := getSamHeader(infile)
	// if err != nil {
	// 	return err
	// }

	refA, _, err := fastaio.PopulateByteArrayGetNames(reffile)
	if err != nil {
		return err
	}

	refLen := len(refA[0])
	// fmt.Println(refLen)

	err = checkArgs(refLen, trim, pad, trimstart, trimend)
	if err != nil {
		return err
	}

	cSR := make(chan samRecords, threads)
	cReadDone := make(chan bool)

	cFR := make(chan fastaio.FastaRecord)
	cWriteDone := make(chan bool)

	cErr := make(chan error)

	cWaitGroupDone := make(chan bool)

	go groupSamRecords(infile, cSR, cReadDone, cErr)

	go writeAlignmentOut(cFR, outfile, cWriteDone, cErr)

	var wg sync.WaitGroup
	wg.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			blockToFastaRecord(cSR, cFR, cErr, refLen, trim, pad, trimstart, trimend, false)
			wg.Done()
		}()
	}

	go func() {
		wg.Wait()
		cWaitGroupDone <- true
	}()

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cReadDone:
			close(cSR)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cWaitGroupDone:
			close(cFR)
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
