package sam

import (
	"errors"
	"io"
	"sync"

	"github.com/virus-evolution/gofasta/pkg/fastaio"

	biogosam "github.com/biogo/hts/sam"
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
				if i < trimstart-1 || i >= trimend {
					seq[i] = 'N'
				}
			}
		} else {
			seq = seq[trimstart-1 : trimend]
		}
	}

	FR := fastaio.FastaRecord{ID: id, Description: id, Seq: string(seq), Idx: idx}

	return FR
}

// checkArgs sanity checks the trimming and padding arguments, given the length of the reference sequence
func checkArgs(refLen int, trimstart int, trimend int) (int, int, bool, error) {

	trim := false

	if trimstart == -1 {
		trimstart = 1
	} else {
		trim = true
	}

	if trimend == -1 {
		trimend = refLen
	} else {
		trim = true
	}

	if trimstart > refLen || trimstart < 1 {
		return 0, 0, false, errors.New("error parsing --start coordinate. --start must be > 0 && <= length(reference)")
	}
	if trimend > refLen || trimend < 1 {
		return 0, 0, false, errors.New("error parsing --end coordinate. --end must be > 0 &! > length(reference)")
	}
	if trimstart > trimend {
		return 0, 0, false, errors.New("error parsing trimming coordinates: --start must be <= --end")
	}

	return trimstart, trimend, trim, nil
}

// blockToFastaRecord is a worker function that takes items from a channel of sam block structs (with indices)
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

// ToMultiAlign converts a SAM file to a fasta-format alignment.
// Insertions relative to the reference are discarded, so all the sequences are the same (=reference) length
func ToMultiAlign(samIn io.Reader, out io.Writer, trimstart int, trimend int, pad bool, threads int) error {

	cSR := make(chan samRecords, threads)
	cReadDone := make(chan bool)

	cSH := make(chan biogosam.Header)

	cFR := make(chan fastaio.FastaRecord)
	cWriteDone := make(chan bool)

	cErr := make(chan error)

	cWaitGroupDone := make(chan bool)

	go groupSamRecords(samIn, cSH, cSR, cReadDone, cErr)

	header := <-cSH
	refLen := header.Refs()[0].Len()

	trimstart, trimend, trim, err := checkArgs(refLen, trimstart, trimend)
	if err != nil {
		return err
	}

	go fastaio.WriteAlignment(cFR, out, cWriteDone, cErr)

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
			close(cSH)
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
