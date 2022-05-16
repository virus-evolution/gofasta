package closest

import (
	"errors"
	"fmt"
	"io"
	"math"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/fastaio"
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

// catchmentStruct contains information about the closest sequences to a particular query
type catchmentStruct struct {
	qname                string
	qidx                 int
	catchment            []resultsStruct
	furthestDistance     float64 // this is distance for the least close of the current set of neighbours in catchment
	furthestCompleteness int64   // this is completeness for the least close of the current set of neighbours in catchment
}

// rearrangeCatchment sorts a catchmentStruct so that the sequences in its catchment field are in order of
// genetic distance, with ties broken by genome completeness. It is called before the catchment is written to
// output, or if a new sequence is added to a catchmentStruct which is already at capacity
func rearrangeCatchment(nS *catchmentStruct, catchmentSize int) {
	sort.SliceStable(nS.catchment, func(i, j int) bool {
		return nS.catchment[i].distance < nS.catchment[j].distance || (nS.catchment[i].distance == nS.catchment[j].distance && nS.catchment[i].completeness > nS.catchment[j].completeness)
	})
	nS.catchment = nS.catchment[0:catchmentSize]
	nS.furthestDistance = nS.catchment[catchmentSize-1].distance
	nS.furthestCompleteness = nS.catchment[catchmentSize-1].completeness
}

// findClosestN finds the closest sequences by genetic distance to single a query sequence
func findClosestN(query fastaio.EncodedFastaRecord, catchmentSize int, maxdist float64, measure string, cIn chan fastaio.EncodedFastaRecord, cOut chan catchmentStruct) {

	neighbours := catchmentStruct{qname: query.ID, qidx: query.Idx}
	neighbours.catchment = make([]resultsStruct, 0)

	var rs resultsStruct

	var distance float64

	for target := range cIn {

		switch measure {
		case "raw":
			distance = rawDistance(query, target)
		case "snp":
			distance = snpDistance(query, target)
		case "tn93":
			distance = tn93Distance(query, target)
		}

		if maxdist != -1.0 {
			if distance > maxdist {
				continue
			}
		}

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
	if len(neighbours.catchment) < catchmentSize && len(neighbours.catchment) > 0 {
		rearrangeCatchment(&neighbours, len(neighbours.catchment))
	}

	cOut <- neighbours
}

// splitInputN fans out target sequences over an array of query sequences, so that each target is passed over each query.
func splitInputN(queries []fastaio.EncodedFastaRecord, catchmentSize int, maxdist float64, measure string, cIn chan fastaio.EncodedFastaRecord, cOut chan catchmentStruct, cErr chan error, cSplitDone chan bool) {

	nQ := len(queries)

	// make an array of channels, one for each query
	QChanArray := make([]chan fastaio.EncodedFastaRecord, nQ)
	for i := 0; i < nQ; i++ {
		QChanArray[i] = make(chan fastaio.EncodedFastaRecord)
	}

	for i, q := range queries {
		go findClosestN(q, catchmentSize, maxdist, measure, QChanArray[i], cOut)
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

// writeClosestN parses an array of catchmentStructs in order to write them, usually to stdout or file
func writeClosestN(results []catchmentStruct, w io.Writer) error {

	var err error

	_, err = w.Write([]byte("query,closest\n"))
	if err != nil {
		return err
	}

	for _, result := range results {
		temp := make([]string, 0)
		for _, hit := range result.catchment {
			temp = append(temp, hit.tname)
		}
		w.Write([]byte(result.qname + "," + strings.Join(temp, ";") + "\n"))
	}

	return nil
}

func writeClosestNTable(results []catchmentStruct, w io.Writer, measure string) error {

	var err error

	_, err = w.Write([]byte("query,target,distance\n"))
	if err != nil {
		return err
	}

	switch measure {
	case "snp":
		for _, result := range results {
			for _, hit := range result.catchment {
				w.Write([]byte(result.qname + "," + hit.tname + "," + strconv.Itoa(int(hit.distance)) + "\n"))
			}
		}
	default:
		for _, result := range results {
			for _, hit := range result.catchment {
				w.Write([]byte(result.qname + "," + hit.tname + "," + strconv.FormatFloat(hit.distance, 'f', 9, 64) + "\n"))
			}
		}
	}

	return nil
}

// ClosestN finds the closest sequence(s) by genetic distance to a query/queries. It writes the results
// to stdout or to file. Ties for distance are broken by genome completeness.
func ClosestN(catchmentSize int, maxdist float64, query, target io.Reader, measure string, out io.Writer, table bool, threads int) error {

	if threads == 0 {
		threads = runtime.NumCPU()
	} else if threads < runtime.NumCPU() {
		runtime.GOMAXPROCS(threads)
	}

	if maxdist != -1.0 && catchmentSize == 0 {
		catchmentSize = math.MaxInt
	}

	queries, err := fastaio.ReadEncodeAlignmentToList(query, false)
	if err != nil {
		return err
	}

	nQ := len(queries)

	fmt.Fprintf(os.Stderr, "number of sequences in query alignment: %d\n", nQ)

	QResultsArray := make([]catchmentStruct, nQ)

	cErr := make(chan error)

	cTEFR := make(chan fastaio.EncodedFastaRecord, runtime.NumCPU())
	cTEFRdone := make(chan bool)
	cSplitDone := make(chan bool)
	cResults := make(chan catchmentStruct)

	go fastaio.ReadEncodeScoreAlignment(target, false, cTEFR, cErr, cTEFRdone)

	go splitInputN(queries, catchmentSize, maxdist, measure, cTEFR, cResults, cErr, cSplitDone)

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
		case <-cSplitDone:
			n--
		}
	}

	for i := 0; i < nQ; i++ {
		result := <-cResults
		QResultsArray[result.qidx] = result
	}

	switch table {
	case true:
		err = writeClosestNTable(QResultsArray, out, measure)
	case false:
		err = writeClosestN(QResultsArray, out)
	}
	if err != nil {
		return err
	}

	return nil
}
