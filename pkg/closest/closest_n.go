package closest

import (
	"errors"
	"io"
	"math"
	"runtime"
	"sort"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/fasta"
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
		return lessDistance(nS.catchment[i].distance, nS.catchment[j].distance) || (nS.catchment[i].distance == nS.catchment[j].distance && nS.catchment[i].completeness > nS.catchment[j].completeness)
	})
	nS.catchment = nS.catchment[0:catchmentSize]
	nS.furthestDistance = nS.catchment[catchmentSize-1].distance
	nS.furthestCompleteness = nS.catchment[catchmentSize-1].completeness
}

// consider adds target to the catchment if there is room, or if it is closer than the furthest member
// (or as close and more complete). A distance that is undefined (NaN) or greater than maxdist (unless
// maxdist is -1) never qualifies.
func (neighbours *catchmentStruct) consider(target fasta.EncodedRecord, distance float64, catchmentSize int, maxdist float64) {
	if maxdist != -1.0 {
		// an undefined distance (NaN) can't be said to be within maxdist
		if math.IsNaN(distance) || distance > maxdist {
			return
		}
	}
	rs := resultsStruct{tname: target.ID, completeness: target.Score, distance: distance}
	switch {
	case len(neighbours.catchment) < catchmentSize:
		neighbours.catchment = append(neighbours.catchment, rs)
		if len(neighbours.catchment) == catchmentSize {
			rearrangeCatchment(neighbours, catchmentSize)
		}
	case lessDistance(distance, neighbours.furthestDistance),
		distance == neighbours.furthestDistance && target.Score > neighbours.furthestCompleteness:
		neighbours.catchment = append(neighbours.catchment, rs)
		rearrangeCatchment(neighbours, catchmentSize)
	}
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
func ClosestN(catchmentSize int, maxdist float64, query, target io.Reader, measure string, out io.Writer, table bool, threads int, strict bool) error {

	if catchmentSize < 0 {
		return errors.New("the number of neighbours to find must not be negative")
	}

	if threads == 0 {
		threads = runtime.NumCPU()
	} else if threads < runtime.NumCPU() {
		runtime.GOMAXPROCS(threads)
	}

	if maxdist != -1.0 && catchmentSize == 0 {
		catchmentSize = math.MaxInt
	}

	queries, err := loadQueries(query, strict)
	if err != nil {
		return err
	}

	nQ := len(queries)

	QResultsArray := make([]catchmentStruct, nQ)
	for i, q := range queries {
		QResultsArray[i] = catchmentStruct{qname: q.ID, qidx: q.Idx, catchment: make([]resultsStruct, 0)}
	}

	err = scanTargets(queries, target, measure, threads, strict, func(st scoredTarget) {
		for i := range QResultsArray {
			QResultsArray[i].consider(st.target, st.distances[i], catchmentSize, maxdist)
		}
	})
	if err != nil {
		return err
	}

	// a catchment that never filled up has not been sorted yet
	for i := range QResultsArray {
		if n := len(QResultsArray[i].catchment); n < catchmentSize && n > 0 {
			rearrangeCatchment(&QResultsArray[i], n)
		}
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
