package updown


import (
	"os"
	"fmt"
	"math"
	"sort"
	"errors"
	"strings"
	"path/filepath"

	"github.com/cov-ert/gofasta/pkg/fastaio"
)

/*
           Venn diagram of SNPs

Q: in query
T: in target
QT: in query and target

          _____             _____
   X‾‾‾‾‾‾     ‾‾‾‾xxxx‾‾‾‾‾     ‾‾‾‾X
  XX            XX     XXX           XXX
 XX           XXX        XXX            XX
X             X            X             |
|             |            |             |
|     Q       |     QT     |      T      |
|             |            |             |
|             X            X             |
XX            XXX         XX            XX
 XXX            XXX    XXX            XXX
   xx-----_____----XXXX----_____-----xx


Table of relationships

   Q   |  QT  |   T   |  interpretation           |   bucket for T
-------------------------------------------------------------------
   0       0      0      on a polytomy (with ref)     same
   0      >0      0      on a polytomy                same
  >0      >0      0      Q is child of T              up
  >0       0      0      Q is child of T              up
   0      >0     >0      T is child of Q              down
   0       0     >0      T is child of Q              down
  >0      >0     >0      Q and T are siblings         side
  >0       0     >0      Q and T are orphans          side?

*/

func index(vs []string, t string) int {
    for i, v := range vs {
        if v == t {
            return i
        }
    }
    return -1
}

func snpOverlap(vs []string, t string) bool {
    return index(vs, t) >= 0
}

func isSiteAmb(pos int, a []int) bool {
	for i := 0; i < len(a); i+=2 {
		if pos >= a[i] && pos <= a[i+1] {
			return true
		}
	}
	return false
}

// whichWay returns direction values of 0,1,2,3 = same,up,down,side respectively + SNP distance
// distance is -1 if the pair fails the ambiguity threshold test.
func whichWay(q, t updownLine, thresh float32) (int, int) {

	// table has 4 items: number of Q SNPs; no. QT SMPs; no. T SNPs; no of consequential ambiguous sites for this pair
	var table [4]int
	for i, qsnp := range(q.snps) {
		if isSiteAmb(q.snpsPos[i], t.ambs) {
			table[3]++
		} else if snpOverlap(t.snps, qsnp) {
			table[1]++
		} else {
			table[0]++
		}
	}
	for i, tsnp := range(t.snps) {
		if isSiteAmb(t.snpsPos[i], q.ambs) {
			table[3]++
		} else if !snpOverlap(q.snps, tsnp)  {
			table[2]++
		}
	}

	sum := 0
	for _, c := range(table) {
		sum+=c
	}
	if (float32(table[3]) / float32(sum)) > thresh {
		return 0, -1
	}

	// Table of relationships
	//
	//  t[0]  | t[1] | t[2]  |  interpretation           |   bucket for T
	// -------------------------------------------------------------------
	//    0       0      0      on a polytomy (with ref)     same
	//    0      >0      0      on a polytomy                same
	//   >0      >0      0      Q is child of T              up
	//   >0       0      0      Q is child of T              up
	//    0      >0     >0      T is child of Q              down
	//    0       0     >0      T is child of Q              down
	//   >0      >0     >0      Q and T are siblings         side
	//   >0       0     >0      Q and T are orphans          side?
	var direction int
	switch {
	case table[0] == 0 && table[2] == 0:
		direction = 0
	case table[0] > 0 && table[2] == 0:
		direction = 1
	case table[0] == 0 && table[2] > 0:
		direction = 2
	case table[0] > 0 && table[2] > 0:
		direction = 3
	}

	distance := table[0] + table[2]

	return direction, distance
}

type resultsStruct struct {
	qname string
	qidx int
	tname string
	distance int
	ambCount int
}

type updownCatchmentSubStruct struct {
	catchment []resultsStruct
	maxDist int
	minAmbig int
}

type updownCatchmentStruct struct {
	qname string
	qidx int
	same updownCatchmentSubStruct
	up updownCatchmentSubStruct
	down updownCatchmentSubStruct
	side updownCatchmentSubStruct
}

func rearrangeCatchment(nS *updownCatchmentSubStruct, catchmentSize int) {
	sort.SliceStable(nS.catchment, func(i, j int) bool {
		return nS.catchment[i].distance < nS.catchment[j].distance || (nS.catchment[i].distance == nS.catchment[j].distance && nS.catchment[i].ambCount < nS.catchment[j].ambCount )
	})
	nS.catchment = nS.catchment[0:catchmentSize]
	nS.maxDist = nS.catchment[catchmentSize - 1].distance
	nS.minAmbig = nS.catchment[catchmentSize - 1].ambCount
}

func findUpDownCatchment(q updownLine, sizetotal int, thresh float32, cIn chan updownLine, cOut chan updownCatchmentStruct) {

	neighbours := updownCatchmentStruct{qname: q.id, qidx: q.idx}
	neighbours.same = updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)}
	neighbours.up = updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)}
	neighbours.down = updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)}
	neighbours.side = updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)}

	var rs resultsStruct
	var distance int
	var direction int

	for target := range(cIn) {

		// return direction values of 0,1,2,3 = same,up,down,side respectively
		// distance is SNP-distance (int)
		direction, distance = whichWay(q, target, thresh)

		// if it doesn't pass the ambiguity threshold, distance will be -1
		if distance < 0 {
			continue
		}

		switch direction {
		case 0: // same
			if len(neighbours.same.catchment) < sizetotal {
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				neighbours.same.catchment = append(neighbours.same.catchment, rs)

				if len(neighbours.same.catchment) == sizetotal {
					rearrangeCatchment(&neighbours.same, sizetotal)
				}

			} else if distance < neighbours.same.maxDist {
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				neighbours.same.catchment = append(neighbours.same.catchment, rs)
				rearrangeCatchment(&neighbours.same, sizetotal)

			} else if distance == neighbours.same.maxDist && target.ambCount < neighbours.same.minAmbig {
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				neighbours.same.catchment = append(neighbours.same.catchment, rs)
				rearrangeCatchment(&neighbours.same, sizetotal)
			}
		case 1: // up
			if len(neighbours.up.catchment) < sizetotal {
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				neighbours.up.catchment = append(neighbours.up.catchment, rs)

				if len(neighbours.up.catchment) == sizetotal {
					rearrangeCatchment(&neighbours.up, sizetotal)
				}

			} else if distance < neighbours.up.maxDist {
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				neighbours.up.catchment = append(neighbours.up.catchment, rs)
				rearrangeCatchment(&neighbours.up, sizetotal)

			} else if distance == neighbours.up.maxDist && target.ambCount < neighbours.up.minAmbig {
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				neighbours.up.catchment = append(neighbours.up.catchment, rs)
				rearrangeCatchment(&neighbours.up, sizetotal)
			}
		case 2: // down
			if len(neighbours.down.catchment) < sizetotal {
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				neighbours.down.catchment = append(neighbours.down.catchment, rs)

				if len(neighbours.down.catchment) == sizetotal {
					rearrangeCatchment(&neighbours.down, sizetotal)
				}

			} else if distance < neighbours.down.maxDist {
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				neighbours.down.catchment = append(neighbours.down.catchment, rs)
				rearrangeCatchment(&neighbours.down, sizetotal)

			} else if distance == neighbours.down.maxDist && target.ambCount < neighbours.down.minAmbig {
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				neighbours.down.catchment = append(neighbours.down.catchment, rs)
				rearrangeCatchment(&neighbours.down, sizetotal)
			}
		case 3: // side
			if len(neighbours.side.catchment) < sizetotal {
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				neighbours.side.catchment = append(neighbours.side.catchment, rs)

				if len(neighbours.side.catchment) == sizetotal {
					rearrangeCatchment(&neighbours.side, sizetotal)
				}

			} else if distance < neighbours.side.maxDist {
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				neighbours.side.catchment = append(neighbours.side.catchment, rs)
				rearrangeCatchment(&neighbours.side, sizetotal)

			} else if distance == neighbours.side.maxDist && target.ambCount < neighbours.side.minAmbig {
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				neighbours.side.catchment = append(neighbours.side.catchment, rs)
				rearrangeCatchment(&neighbours.side, sizetotal)
			}
		}
	}

	if len(neighbours.same.catchment) < sizetotal && len(neighbours.same.catchment) > 0 {
		rearrangeCatchment(&neighbours.same, len(neighbours.same.catchment))
	}
	if len(neighbours.up.catchment) < sizetotal && len(neighbours.up.catchment) > 0 {
		rearrangeCatchment(&neighbours.up, len(neighbours.up.catchment))
	}
	if len(neighbours.down.catchment) < sizetotal && len(neighbours.down.catchment) > 0 {
		rearrangeCatchment(&neighbours.down, len(neighbours.down.catchment))
	}
	if len(neighbours.side.catchment) < sizetotal && len(neighbours.side.catchment) > 0 {
		rearrangeCatchment(&neighbours.side, len(neighbours.side.catchment))
	}

	cOut<- neighbours
}

func splitInput(queries []updownLine, sizetotal int, threshsnp float32, threshtarg int,
	cIn chan updownLine, cOut chan updownCatchmentStruct, cErr chan error, cSplitDone chan bool) {

	nQ := len(queries)

	// make an array of channels, one for each query
	QChanArray := make([]chan updownLine, nQ)
	for i := 0; i < nQ; i++ {
	   QChanArray[i] = make(chan updownLine)
	}

	for i, q := range(queries) {
		go findUpDownCatchment(q, sizetotal, threshsnp, QChanArray[i], cOut)
	}

	for udL := range(cIn) {
		if udL.ambCount > threshtarg {
			continue
		}
		for i, _ := range(QChanArray) {
			QChanArray[i]<- udL
		}
	}

	for i, _ := range(QChanArray) {
		close(QChanArray[i])
	}

	cSplitDone<- true
}

func allGreaterThanEqualTo4(a, b [4]int) bool {
	for i, _ := range(a) {
		if a[i] < b[i] {
			return false
		}
	}
	return true
}

// given the input options, and the final size of the neighbourhoods,
// return the size that the output data should be
func balance(sizetotal int, sizeIdeal [4]int, sizeObserved [4]int, nofill bool) [4]int {
	var size [4]int

	allGreaterThan := allGreaterThanEqualTo4(sizeObserved, sizeIdeal)

	if allGreaterThan {
		return sizeIdeal
	} else {
		var sizeAvail [4]int
		for i, _ := range(sizeObserved) {
			if sizeObserved[i] > sizeIdeal[i] {
				sizeAvail[i] = sizeObserved[i] - sizeIdeal[i]
			} else {
				sizeAvail[i] = 0
			}
		}
		for i, _ := range(sizeObserved) {
			if sizeObserved[i] >= sizeIdeal[i] {
				size[i] = sizeIdeal[i]
			} else {
				size[i] = sizeObserved[i]
			}
		}
		if nofill {
			return size
		}
		for n := 1; n > 0; {
			for i, _ := range(sizeObserved) {
				if sum4(sizeAvail) == 0 {
					n--
					break
				}
				if sizeObserved[i] > sizeIdeal[i] && sizeAvail[i] > 0 {
					size[i]++
					sizeAvail[i]--
				}
				if sum4(size) == sizetotal {
					n--
					break
				}
			}
		}
	}

	return size
}

func writeUpDownCatchment(fileOut string, results []updownCatchmentStruct, sizetotal int, sizeArray [4]int, nofill bool) error {
	// write to stdout for now...

	var err error
	f := os.Stdout

	if fileOut != "stdout" {
		f, err = os.Create(fileOut)
		if err != nil {
			return err
		}
	}

	defer f.Close()

	var temp []string
	var size [4]int
	var sizeObserved [4]int

	_, err = f.WriteString("query,closestsame,closestup,closestdown,closestside\n")
	if err != nil {
		return err
	}

	for _, result := range results {

		sizeObserved[0] = len(result.same.catchment)
		sizeObserved[1] = len(result.up.catchment)
		sizeObserved[2] = len(result.down.catchment)
		sizeObserved[3] = len(result.side.catchment)
		size = balance(sizetotal, sizeArray, sizeObserved, nofill)

		_, err = f.WriteString(result.qname + ",")
		if err != nil {
			return err
		}

		temp = make([]string, 0)
		for _, hit := range(result.same.catchment[0:size[0]]) {
			temp = append(temp, hit.tname)
		}
		_, err = f.WriteString(strings.Join(temp, ";") + ",")
		if err != nil {
			return err
		}

		temp = make([]string, 0)
		for _, hit := range(result.up.catchment[0:size[1]]) {
			temp = append(temp, hit.tname)
		}
		_, err = f.WriteString(strings.Join(temp, ";") + ",")
		if err != nil {
			return err
		}

		temp = make([]string, 0)
		for _, hit := range(result.down.catchment[0:size[2]]) {
			temp = append(temp, hit.tname)
		}
		_, err = f.WriteString(strings.Join(temp, ";") + ",")
		if err != nil {
			return err
		}

		temp = make([]string, 0)
		for _, hit := range(result.side.catchment[0:size[3]]) {
			temp = append(temp, hit.tname)
		}
		_, err = f.WriteString(strings.Join(temp, ";") + "\n")
		if err != nil {
			return err
		}
	}

	return nil
}

func allZero(s []int) bool {
	for _, n := range(s) {
		if n != 0 {
			return false
		}
	}
	return true
}

func checkArgs(query string, target string, reference string,
	sizetotal int, sizeup int, sizedown int, sizeside int, sizesame int) ([4]int, string, string, error) {

	var qtype string
	var ttype string

	switch filepath.Ext(query) {
	case ".csv":
		qtype = "csv"
	case ".fasta":
		qtype = "fasta"
	case ".fa":
		qtype = "fasta"
	default:
		return [4]int{}, qtype, ttype, errors.New("couldn't tell if --query was a .csv or a .fasta file")
	}

	switch filepath.Ext(target) {
	case ".csv":
		ttype = "csv"
	case ".fasta":
		ttype = "fasta"
	case ".fa":
		ttype = "fasta"
	default:
		return [4]int{}, qtype, ttype, errors.New("couldn't tell if --target was a .csv or a .fasta file")
	}

	if (qtype == "fasta" || ttype == "fasta") && len(reference) == 0 {
		return [4]int{}, qtype, ttype, errors.New("if your input files are fastas, you must provide a --reference")
	}

	// return bucket sizes for same,up,down,side,total respectively
	var sizeArray [4]int

	if sizetotal > 0 {
		sizeArray[1] = int(math.Floor(float64(sizetotal) / float64(4)))
		sizeArray[2] = int(math.Floor(float64(sizetotal) / float64(4)))
		sizeArray[3] = int(math.Floor(float64(sizetotal) / float64(4)))
		sizeArray[0] = sizetotal - (sizeArray[1] + sizeArray[2] + sizeArray[3])
	} else {
		sizeArray[0] = sizesame
		sizeArray[1] = sizeup
		sizeArray[2] = sizedown
		sizeArray[3] = sizeside
	}

	if sizetotal == 0 && allZero([]int{sizeup, sizedown, sizeside, sizesame}) {
		return sizeArray, qtype, ttype, errors.New("Please provide values to either --sizetotal or all or some of --sizeup, --sizedown, --sizeside and --sizesame")
	}

	if sizetotal != 0 && ! allZero([]int{sizeup, sizedown, sizeside, sizesame}) {
		return sizeArray, qtype, ttype, errors.New("warning: setting --size-total overrides --size-up, --size-down, --size-side and --size-same")
	}

	return sizeArray, qtype, ttype, nil
}

func sum4(a [4]int) int {
	var t int
	for _, n := range(a) {
		t+=n
	}
	return t
}

func TopRanking(query string, target string, outfile string, reference string,
	    sizetotal int, sizeup int, sizedown int, sizeside int, sizesame int,
	    threshsnp float32, threshtarg int, nofill bool) error {

	_ = fmt.Println

	sizeArray, q_in_type, t_in_type, err := checkArgs(query, target, reference, sizetotal, sizeup, sizedown, sizeside, sizesame)
	if err != nil {
		if err.Error() == "warning: setting --size-total overrides --size-up, --size-down, --size-side and --size-same" {
			os.Stderr.WriteString(err.Error() + "\n")
		} else {
			return err
		}
	}

	var refSeq []byte
	if q_in_type == "fasta" || t_in_type == "fasta" {
		temp, err := fastaio.ReadEncodeAlignmentToList(reference)
		if err != nil {
			return err
		}
		if len(temp) > 1 {
			return errors.New("More than one record in --reference")
		}
		refSeq = temp[0].Seq
	}

	cErr := make(chan error)
	var queries []updownLine

	switch q_in_type {
	case "csv":
		queries, err = readCSVToList(query)
		if err != nil {
			return err
		}
	case "fasta":

		queries, err = fastaToUDLslice(query, refSeq)
		if err != nil {
			return err
		}

	}

	nQ := len(queries)
	QResultsArray := make([]updownCatchmentStruct, nQ)

	cudL := make(chan updownLine)
	cReadDone := make(chan bool)

	cResults := make(chan updownCatchmentStruct)
	cSplitDone := make(chan bool)

	switch t_in_type {
	case "csv":
		go readCSVToChan(target, cudL, cErr, cReadDone)
	case "fasta":
		go readFastaToChan(target, refSeq, cudL, cErr, cReadDone)
	}

	go splitInput(queries,
		sum4(sizeArray), threshsnp, threshtarg,
		cudL, cResults, cErr, cSplitDone)

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cReadDone:
			close(cudL)
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

	err = writeUpDownCatchment(outfile, QResultsArray, sum4(sizeArray), sizeArray, nofill)
	if err != nil {
		return err
	}

	return nil
}
