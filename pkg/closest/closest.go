/*
Package closest provides routines to find the closest sequences
to a set of query sequences, by genetic distance
*/
package closest

import (
	"errors"
	"fmt"
	"io"
	"os"
	"runtime"
	"strconv"
	"strings"
	"sync"

	"github.com/cov-ert/gofasta/pkg/encoding"
	"github.com/cov-ert/gofasta/pkg/fastaio"
)

// resultsStruct is a struct that contains information about a query sequence and its (current)
// single closest target by raw genetic distance, including the snps that distinguish them and
// the total snp-distance between them.
type resultsStruct struct {
	qname        string
	qidx         int
	tname        string
	completeness int64
	distance     float64
	snps         []string
}

// scoreEncodedAlignment scores a single fasta record by how complete its genome is.
func scoreEncodedAlignment(cIn chan fastaio.EncodedFastaRecord, cOut chan fastaio.EncodedFastaRecord) {
	scoring := encoding.MakeEncodedScoreArray()
	var score int64

	for EFR := range cIn {
		score = 0
		for _, nuc := range EFR.Seq {
			score += scoring[nuc]
		}
		EFR.Score = score
		cOut <- EFR
	}

	return
}

// findClosest finds the single closest sequence by genetic distance among a set of target sequences to a query sequence
func findClosest(query fastaio.EncodedFastaRecord, cIn chan fastaio.EncodedFastaRecord, cOut chan resultsStruct) {
	var closest resultsStruct
	var distance float64
	var snps []string

	var n int64
	var d int64

	first := true

	decoding := encoding.MakeDecodingArray()

	for target := range cIn {
		n = 0
		d = 0
		for i, tNuc := range target.Seq {
			if (query.Seq[i] & tNuc) < 16 {
				n += 1
				d += 1
			}
			if (query.Seq[i]&8 == 8) && query.Seq[i] == tNuc {
				d += 1
			}
		}
		distance = float64(n) / float64(d)

		if first {
			snps = make([]string, 0)
			for i, tNuc := range target.Seq {
				if (query.Seq[i] & tNuc) < 16 {
					snps = append(snps, strconv.Itoa(i+1)+decoding[query.Seq[i]]+decoding[tNuc])
				}
			}
			closest = resultsStruct{tname: target.ID, completeness: target.Score, distance: distance, snps: snps}
			first = false
			continue
		}

		if distance < closest.distance {
			snps = make([]string, 0)
			for i, tNuc := range target.Seq {
				if (query.Seq[i] & tNuc) < 16 {
					snps = append(snps, strconv.Itoa(i+1)+decoding[query.Seq[i]]+decoding[tNuc])
				}
			}
			closest = resultsStruct{tname: target.ID, completeness: target.Score, distance: distance, snps: snps}

		} else if distance == closest.distance {
			if target.Score > closest.completeness {
				snps = make([]string, 0)
				for i, tNuc := range target.Seq {
					if (query.Seq[i] & tNuc) < 16 {
						snps = append(snps, strconv.Itoa(i+1)+decoding[query.Seq[i]]+decoding[tNuc])
					}
				}
				closest = resultsStruct{tname: target.ID, completeness: target.Score, distance: distance, snps: snps}
			}
		}
	}

	closest.qname = query.ID
	closest.qidx = query.Idx

	cOut <- closest
}

// splitInput fans out target sequences over an array of query sequences, so that each target is passed over each query.
func splitInput(queries []fastaio.EncodedFastaRecord, cIn chan fastaio.EncodedFastaRecord, cOut chan resultsStruct, cErr chan error, cSplitDone chan bool) {

	nQ := len(queries)

	// make an array of channels, one for each query
	QChanArray := make([]chan fastaio.EncodedFastaRecord, nQ)
	for i := 0; i < nQ; i++ {
		QChanArray[i] = make(chan fastaio.EncodedFastaRecord)
	}

	for i, q := range queries {
		go findClosest(q, QChanArray[i], cOut)
	}

	targetCounter := 0
	for EFR := range cIn {
		if targetCounter == 0 {
			if len(EFR.Seq) != len(queries[0].Seq) {
				cErr <- errors.New("query and target alignments are not the same width")
			}
		}
		targetCounter++

		for i, _ := range QChanArray {
			QChanArray[i] <- EFR
		}
	}

	fmt.Fprintf(os.Stderr, "number of sequences in target alignment: %d\n", targetCounter)

	for i, _ := range QChanArray {
		close(QChanArray[i])
	}

	cSplitDone <- true
}

// writeClosest parses an array of resultsStructs in order to write them, usually to stdout or file
func writeClosest(results []resultsStruct, w io.Writer) error {

	var err error

	_, err = w.Write([]byte("query,closest,SNPdistance,SNPs\n"))
	if err != nil {
		return err
	}

	for _, result := range results {
		w.Write([]byte(result.qname + "," + result.tname + "," + strconv.Itoa(len(result.snps)) + "," + strings.Join(result.snps, ";") + "\n"))
	}

	return nil
}

// Closest finds the single closest sequence by raw genetic distance to a query/queries. It writes the results
// to stdout or to file. Ties for distance are broken by genome completeness
func Closest(query, target io.Reader, out io.Writer, threads int) error {

	if threads == 0 {
		threads = runtime.NumCPU()
	} else if threads < runtime.NumCPU() {
		runtime.GOMAXPROCS(threads)
	}

	queries, err := fastaio.ReadEncodeAlignmentToList(query, false)
	if err != nil {
		return err
	}

	nQ := len(queries)

	fmt.Fprintf(os.Stderr, "number of sequences in query alignment: %d\n", nQ)

	QResultsArray := make([]resultsStruct, nQ)

	cErr := make(chan error)

	cTEFR := make(chan fastaio.EncodedFastaRecord, runtime.NumCPU())
	cTEFRscored := make(chan fastaio.EncodedFastaRecord, runtime.NumCPU())
	cTEFRdone := make(chan bool)
	cTEFRscoreddone := make(chan bool)
	cSplitDone := make(chan bool)

	cResults := make(chan resultsStruct)

	go fastaio.ReadEncodeAlignment(target, false, cTEFR, cErr, cTEFRdone)

	var wgScore sync.WaitGroup
	wgScore.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			scoreEncodedAlignment(cTEFR, cTEFRscored)
			wgScore.Done()
		}()
	}

	go splitInput(queries, cTEFRscored, cResults, cErr, cSplitDone)

	go func() {
		wgScore.Wait()
		cTEFRscoreddone <- true
	}()

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cTEFRdone:
			close(cTEFR)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cTEFRscoreddone:
			close(cTEFRscored)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cSplitDone:
			n--
		}
	}

	for i := 0; i < nQ; i++ {
		result := <-cResults
		QResultsArray[result.qidx] = result
	}

	err = writeClosest(QResultsArray, out)
	if err != nil {
		return err
	}

	return nil
}
