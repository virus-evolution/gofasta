package closest

import (
	"os"
	"fmt"
	"sort"
	"sync"
	"errors"
	"runtime"
	"strings"

	"github.com/cov-ert/gofasta/pkg/fastaio"
)

// this is defined elsewhere, but for reference:
// type resultsStruct struct {
// 	qname string
// 	qidx int
// 	tname string
// 	completeness int64
// 	distance float64
// 	snps []string
// }

type catchmentStruct struct {
	qname string
	qidx int
	catchment []resultsStruct
	furthestDistance float64 // this is distance for the least close of the current set of neighbours in catchment
	furthestCompleteness int64 // this is completeness for the least close of the current set of neighbours in catchment
}

func rearrangeCatchment(nS *catchmentStruct, catchmentSize int) {
	sort.SliceStable(nS.catchment, func(i, j int) bool {
		return nS.catchment[i].distance < nS.catchment[j].distance || (nS.catchment[i].distance == nS.catchment[j].distance && nS.catchment[i].completeness > nS.catchment[j].completeness )
	})
	nS.catchment = nS.catchment[0:catchmentSize]
	nS.furthestDistance = nS.catchment[catchmentSize - 1].distance
	nS.furthestCompleteness = nS.catchment[catchmentSize - 1].completeness
}

func findClosestN(query fastaio.EncodedFastaRecord, catchmentSize int, cIn chan fastaio.EncodedFastaRecord, cOut chan catchmentStruct) {

	neighbours := catchmentStruct{qname: query.ID, qidx: query.Idx}
	neighbours.catchment = make([]resultsStruct, 0)

	var rs resultsStruct

	var n int64
	var d int64
	var distance float64

	for target := range(cIn) {
		n = 0
		d = 0
		for i, tNuc := range(target.Seq) {
			if (query.Seq[i] & tNuc) < 16 {
				n += 1
				d += 1
			}
			if (query.Seq[i] & 8 == 8) && query.Seq[i] == tNuc {
				d += 1
			}
		}
		distance = float64(n) / float64(d)

		if len(neighbours.catchment) < catchmentSize {
			rs = resultsStruct{tname: target.ID, completeness: target.Score, distance: distance}
			neighbours.catchment = append(neighbours.catchment, rs)

			if len(neighbours.catchment) == catchmentSize {
				rearrangeCatchment(&neighbours, catchmentSize)
			}

		} else if distance < neighbours.furthestDistance {
			rs = resultsStruct{tname: target.ID, completeness: target.Score, distance: distance}
			neighbours.catchment = append(neighbours.catchment, rs)
			rearrangeCatchment(&neighbours, catchmentSize)

		} else if distance == neighbours.furthestDistance && target.Score > neighbours.furthestCompleteness {
			rs = resultsStruct{tname: target.ID, completeness: target.Score, distance: distance}
			neighbours.catchment = append(neighbours.catchment, rs)
			rearrangeCatchment(&neighbours, catchmentSize)
		}
	}

	// If the user specified a larger catchment than there are records in the target file,
	// they won't be sorted above, so do it here (need to modify the size argument passed
	// to the function):
	if len(neighbours.catchment) < catchmentSize {
		rearrangeCatchment(&neighbours, len(neighbours.catchment))
	}

	cOut<- neighbours
}

func splitInputN(queries []fastaio.EncodedFastaRecord, catchmentSize int, cIn chan fastaio.EncodedFastaRecord, cOut chan catchmentStruct, cErr chan error, cSplitDone chan bool) {

	nQ := len(queries)

	// make an array of channels, one for each query
	QChanArray := make([]chan fastaio.EncodedFastaRecord, nQ)
	for i := 0; i < nQ; i++ {
	   QChanArray[i] = make(chan fastaio.EncodedFastaRecord)
	}

	for i, q := range(queries) {
		go findClosestN(q, catchmentSize, QChanArray[i], cOut)
	}

	targetCounter := 0
	for EFR := range(cIn) {
		if targetCounter == 0 {
			if len(EFR.Seq) != len(queries[0].Seq) {
				cErr<- errors.New("query and target alignments are not the same width")
			}
		}
		targetCounter++

		for i, _ := range(QChanArray) {
			QChanArray[i]<- EFR
		}
	}

	fmt.Printf("number of sequences in target alignment: %d\n", targetCounter)

	for i, _ := range(QChanArray) {
		close(QChanArray[i])
	}

	cSplitDone<- true
}

func writeClosestN(results []catchmentStruct, filepath string) error {
	f, err := os.Create(filepath)
	if err != nil {
		return err
	}
	defer f.Close()

	_, err = f.WriteString("query,closest\n")
	if err != nil {
		return err
	}

	for _, result := range results {
		temp := make([]string, 0)
		for _, hit := range(result.catchment) {
			temp = append(temp, hit.tname)
		}
		f.WriteString(result.qname + "," + strings.Join(temp, "|") + "\n")
	}

	return nil
}

func ClosestN(catchmentSize int, queryFile string, targetFile string, outFile string, threads int) error {

	if threads == 0 {
		threads = runtime.NumCPU()
	} else if threads < runtime.NumCPU() {
		runtime.GOMAXPROCS(threads)
	}

	queries, err := fastaio.ReadEncodeAlignmentToList(queryFile)
	if err != nil {
		return err
	}

	nQ := len(queries)

	fmt.Printf("number of sequences in query alignment: %d\n", nQ)

	QResultsArray := make([]catchmentStruct, nQ)

	cErr := make(chan error)

	cTEFR := make(chan fastaio.EncodedFastaRecord, runtime.NumCPU())
	cTEFRscored := make(chan fastaio.EncodedFastaRecord, runtime.NumCPU())
	cTEFRdone := make(chan bool)
	cTEFRscoreddone := make(chan bool)
	cSplitDone := make(chan bool)

	cResults := make(chan catchmentStruct)

	go fastaio.ReadEncodeAlignment(targetFile, cTEFR, cErr, cTEFRdone)

	var wgScore sync.WaitGroup
	wgScore.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			scoreEncodedAlignment(cTEFR, cTEFRscored)
			wgScore.Done()
		}()
	}

	go splitInputN(queries, catchmentSize, cTEFRscored, cResults, cErr, cSplitDone)

	go func() {
		wgScore.Wait()
		cTEFRscoreddone<- true
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

	err = writeClosestN(QResultsArray, outFile)
	if err != nil {
		return err
	}

	return nil
}
