/*
Package updown implements functions that leverage pseudo-tree aware
snp-distances between sequences. It relies on there being a suitable
sequence that one can use to represent the "root" (ancestor) of a
hypothetical phylogenetic tree, and for the distance between relevant
sequences to be small enough that maximum parsimony is a reasonable
assumption.
*/
package updown

import (
	"errors"
	"io"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/fastaio"
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

// index returns the index of a query string in an array of target strings, if it exists,
// else it returns -1
func index(vs []string, t string) int {
	for i, v := range vs {
		if v == t {
			return i
		}
	}
	return -1
}

// snpOverlap asks if one snp is present in an array of snps
func snpOverlap(vs []string, t string) bool {
	return index(vs, t) >= 0
}

// indexInt returns the index of a query integer in an array of target integers, if it exists,
// else it returns -1
func indexInt(vs []int, t int) int {
	for i, v := range vs {
		if v == t {
			return i
		}
	}
	return -1
}

// posOverlap asks if one position is present in an array of positions
func posOverlap(vs []int, t int) bool {
	return indexInt(vs, t) >= 0
}

// isSiteAmb asks if a position is within an ambiguity tract, given an array of pairs of start-stop coordinates
// of such tracts
func isSiteAmb(pos int, a []int) bool {
	for i := 0; i < len(a); i += 2 {
		if pos >= a[i] && pos <= a[i+1] {
			return true
		}
	}
	return false
}

// whichWay returns direction values of 0,1,2,3 = same,up,down,side respectively + SNP distance.
// distance is -1 if the pair fails the ambiguity threshold test.
func whichWay(q, t updownLine, thresh float32) (int, int) {

	// table has 4 items: number of Q SNPs; no. QT SNPs; no. T SNPs; no of consequential ambiguous sites for this pair
	var table [4]int

	// d is an array of unique SNP positions to get distance from
	d := make([]int, 0)

	for i, qsnp := range q.snps {
		if isSiteAmb(q.snpsPos[i], t.ambs) {
			table[3]++
		} else if snpOverlap(t.snps, qsnp) {
			table[1]++
		} else {
			table[0]++
			d = append(d, q.snpsPos[i])
		}
	}
	for i, tsnp := range t.snps {
		if isSiteAmb(t.snpsPos[i], q.ambs) {
			table[3]++
		} else if !snpOverlap(q.snps, tsnp) {
			table[2]++
			if !posOverlap(d, t.snpsPos[i]) {
				d = append(d, t.snpsPos[i])
			}
		}
	}

	sum := 0
	for _, c := range table {
		sum += c
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

	// this is wrong (when there are multiple hits):
	// distance := table[0] + table[2]

	distance := len(d)

	return direction, distance
}

// a resultsStruct contains information about the relationship between one query and one target.
type resultsStruct struct {
	qname    string // name of the query
	qidx     int    // query's position in the input file
	tname    string // name of the target
	distance int    // snp distance between query and target
	ambCount int    // number of non-ATGC characters in the target
}

// an updownCatchmentSubStruct contains an array of resultsStructs which are the current closest neighbours of one query
// in one direction. The fields maxDist and minAmbig are used when deciding whether or not the next target should replace
// the last item in catchment
type updownCatchmentSubStruct struct {
	catchment []resultsStruct
	maxDist   int
	minAmbig  int
	// nDists    int
}

// an updownCatchmentStruct contains the current lists of targets in each direction for one query sequence
type updownCatchmentStruct struct {
	qname string
	qidx  int
	same  updownCatchmentSubStruct
	up    updownCatchmentSubStruct
	down  updownCatchmentSubStruct
	side  updownCatchmentSubStruct
	// minDistArray [4]int
}

// rearrangeCatchment sorts the catchment slice in an updownCatchmentSubStruct, pops the last item from it, and updates
// the distance and ambiguity fields in preparation for deciding what to do with the next target. It should only be
// called when the updownCatchmentSubStruct is over its capacity (catchmentSize)
func rearrangeCatchment(nS *updownCatchmentSubStruct, catchmentSize int) {
	sort.SliceStable(nS.catchment, func(i, j int) bool {
		return nS.catchment[i].distance < nS.catchment[j].distance || (nS.catchment[i].distance == nS.catchment[j].distance && nS.catchment[i].ambCount < nS.catchment[j].ambCount)
	})
	nS.catchment = nS.catchment[0:catchmentSize]
	nS.maxDist = nS.catchment[catchmentSize-1].distance
	nS.minAmbig = nS.catchment[catchmentSize-1].ambCount
}

// // how many different snp distances are represented in this struct
// func catchmentNDists(nS *updownCatchmentSubStruct) int {
// 	m := make(map[int]int)
// 	for _, rS := range nS.catchment {
// 		m[rS.distance]++
// 	}

// 	count := 0
// 	for k := range m {
// 		_ = k
// 		count++
// 	}

// 	return count
// }

// a pushCatchmentSubStruct is used to keep a record of the targets at different snp distances away, in case there are no
// targets under/at the threshold distance defined on the command line and the user wants to expand the search to encompass
// the closest sequences in any particular direction. catchmentMap is a map with snp-distances as keys, and arrays of
// resultsStructs as values
type pushCatchmentSubStruct struct {
	catchmentMap map[int][]resultsStruct
	size         int // how many sequences in total
	nDists       int // how many snp distances in total
	maxDist      int // the furthest snp distance
}

// getMaxKey gets the largest value key from a map whose keys are snp-distances
func getMaxKey(m map[int][]resultsStruct) int {
	max := 0
	for k := range m {
		if k > max {
			max = k
		}
	}
	return max
}

// refactorPushCatchment is like rearrangeCatchment but for pushCatchmentSubStructs. It either adds targets to snp distances that
// already exist, or adds new snp distances if they don't and the map isn't at capacity, or if it is but the new snp distance is
// lower than any that already exist. In the latter case, it deletes the key and value associated with the further snp distance
func refactorPushCatchment(pC *pushCatchmentSubStruct, rS resultsStruct, nodeDistance int) {
	if _, ok := pC.catchmentMap[rS.distance]; ok {
		// if this sequence's distance is already present, then add to it
		pC.catchmentMap[rS.distance] = append(pC.catchmentMap[rS.distance], rS)
	} else if len(pC.catchmentMap) == nodeDistance {
		// if it isn't, and the map is at capacity, delete the farthest sequences
		// and add this one  elsewhere
		maxKey := getMaxKey(pC.catchmentMap)
		delete(pC.catchmentMap, maxKey)
		pC.catchmentMap[rS.distance] = []resultsStruct{rS}
		pC.maxDist = getMaxKey(pC.catchmentMap)
		pC.nDists = len(pC.catchmentMap)
	} else {
		// otherwise the map isn't at capacity and we can add this sequence with impunity
		pC.catchmentMap[rS.distance] = []resultsStruct{rS}
		pC.maxDist = getMaxKey(pC.catchmentMap)
		pC.nDists = len(pC.catchmentMap)
	}
}

// pushCatchment2Catchment converts a pushCatchmentSubStruct to an updownCatchmentSubStruct
// TO DO - sort each snp distances results structs as they come out of the map here?
func pushCatchment2Catchment(pC pushCatchmentSubStruct) updownCatchmentSubStruct {
	uC := updownCatchmentSubStruct{}
	uC.catchment = make([]resultsStruct, 0)
	for _, v := range pC.catchmentMap {
		for _, r := range v {
			uC.catchment = append(uC.catchment, r)
		}
	}
	return uC
}

// stringInArray returns true/false s is present in the slice sa
func stringInArray(s string, sa []string) bool {
	for i := range sa {
		if s == sa[i] {
			return true
		}
	}
	return false
}

// findUpDownCatchmentPushDistance does all the work for one query, by iterating over targets as they arrive and assigning then to
// the correct bins based on the results of whichWay. It maintains a set of pushCatchmentSubStructs in case the bins are empty under
// the user-defined snp distance thresholds from the command line, to push the distances out to.
func findUpDownCatchmentPushDistance(q updownLine, ignore []string, sizeArray [4]int, pushDist int, thresh float32, cIn chan updownLine, cOut chan updownCatchmentStruct) {

	var rs resultsStruct
	var distance int
	var direction int

	// the substruct for the sequences on the polytomy
	same := updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)}
	// the substructs for the queries that are the min dist(s) away
	pushup := pushCatchmentSubStruct{catchmentMap: make(map[int][]resultsStruct), nDists: 0}
	pushdown := pushCatchmentSubStruct{catchmentMap: make(map[int][]resultsStruct), nDists: 0}
	pushside := pushCatchmentSubStruct{catchmentMap: make(map[int][]resultsStruct), nDists: 0}

	// then we iterate over all the targets
	for target := range cIn {

		// skip ignored sequences
		if len(ignore) > 0 {
			if stringInArray(target.id, ignore) {
				continue
			}
		}

		// return direction values of 0,1,2,3 = same,up,down,side respectively
		// distance is SNP-distance (int)
		direction, distance = whichWay(q, target, thresh)

		// if it doesn't pass the pairwise ambiguity threshold, distance will be -1
		if distance < 0 {
			continue
		}

		switch direction {
		case 0: // same
			rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
			same.catchment = append(same.catchment, rs)
		case 1: // up
			// is the distance lower or are there not enough distances yet:
			if distance <= pushup.maxDist || pushup.nDists < pushDist {
				// the results struct:
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				// slot it in:
				refactorPushCatchment(&pushup, rs, pushDist)
			}
		case 2: // down
			// is the distance lower or are there not enough distances yet:
			if distance <= pushdown.maxDist || pushdown.nDists < pushDist {
				// the results struct:
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				// slot it in:
				refactorPushCatchment(&pushdown, rs, pushDist)
			}
		case 3: // side
			// is the distance lower or are there not enough distances yet:
			if distance <= pushside.maxDist || pushside.nDists < pushDist {
				// the results struct:
				rs = resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
				// slot it in:
				refactorPushCatchment(&pushside, rs, pushDist)
			}
		}
	}

	neighbours := updownCatchmentStruct{qname: q.id, qidx: q.idx}
	neighbours.same = same

	neighbours.up = pushCatchment2Catchment(pushup)
	sort.SliceStable(neighbours.up.catchment, func(i, j int) bool {
		return neighbours.up.catchment[i].distance < neighbours.up.catchment[j].distance || (neighbours.up.catchment[i].distance == neighbours.up.catchment[j].distance && neighbours.up.catchment[i].ambCount < neighbours.up.catchment[j].ambCount)
	})

	neighbours.down = pushCatchment2Catchment(pushdown)
	sort.SliceStable(neighbours.down.catchment, func(i, j int) bool {
		return neighbours.down.catchment[i].distance < neighbours.down.catchment[j].distance || (neighbours.down.catchment[i].distance == neighbours.down.catchment[j].distance && neighbours.down.catchment[i].ambCount < neighbours.down.catchment[j].ambCount)
	})

	neighbours.side = pushCatchment2Catchment(pushside)
	sort.SliceStable(neighbours.side.catchment, func(i, j int) bool {
		return neighbours.side.catchment[i].distance < neighbours.side.catchment[j].distance || (neighbours.side.catchment[i].distance == neighbours.side.catchment[j].distance && neighbours.side.catchment[i].ambCount < neighbours.side.catchment[j].ambCount)
	})

	// TO DO - balance the neighbours here if sizeArray is given, somehow. Probably will need a balance function specific to pushing

	cOut <- neighbours
}

// findUpDownCatchment does all the work for one query, by iterating over targets as they arrive and assigning then to
// the correct bins based on the results of whichWay. It can't do any pushing if any bins are empty after all targets
// have been processed.
func findUpDownCatchment(q updownLine, ignore []string, sizeArray [4]int, nofill bool, distArray [4]int, thresh float32, cIn chan updownLine, cOut chan updownCatchmentStruct) {

	neighbours := updownCatchmentStruct{qname: q.id, qidx: q.idx}
	neighbours.same = updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)}
	neighbours.up = updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)}
	neighbours.down = updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)}
	neighbours.side = updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)}

	var rs resultsStruct
	var distance int
	var direction int

	// set the total size for all bins to be the maximum out of any bin (for filling in)
	var sizetotal int
	check := false
	for _, n := range sizeArray {
		if n == math.MaxInt32 {
			check = true
		}
	}
	if check {
		sizetotal = math.MaxInt32
	} else {
		sizetotal = sum4(sizeArray)
	}

	// then we iterate over all the targets
	for target := range cIn {

		// skip ignored sequences
		if len(ignore) > 0 {
			if stringInArray(target.id, ignore) {
				continue
			}
		}

		// return direction values of 0,1,2,3 = same,up,down,side respectively
		// distance is SNP-distance (int)
		direction, distance = whichWay(q, target, thresh)

		// if it doesn't pass the pairwise ambiguity threshold, distance will be -1
		if distance < 0 {
			continue
		}

		// if the distance is too big, also skip it (for now)
		if distance > distArray[direction] {
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

	var sizeObserved [4]int
	sizeObserved[0] = len(neighbours.same.catchment)
	sizeObserved[1] = len(neighbours.up.catchment)
	sizeObserved[2] = len(neighbours.down.catchment)
	sizeObserved[3] = len(neighbours.side.catchment)

	size := balance(sizetotal, sizeArray, sizeObserved, nofill)

	neighbours.same.catchment = neighbours.same.catchment[0:size[0]]
	neighbours.up.catchment = neighbours.up.catchment[0:size[1]]
	neighbours.down.catchment = neighbours.down.catchment[0:size[2]]
	neighbours.side.catchment = neighbours.side.catchment[0:size[3]]

	cOut <- neighbours
}

// splitInput fans each target out over the array of queries
func splitInput(queries []updownLine, ignore []string, sizeArray [4]int, nofill bool, distArray [4]int, threshpair float32, threshtarg int,
	pushDistance int, cIn chan updownLine, cOut chan updownCatchmentStruct, cErr chan error, cSplitDone chan bool) {

	nQ := len(queries)

	// make an array of channels, one for each query
	QChanArray := make([]chan updownLine, nQ)
	for i := 0; i < nQ; i++ {
		QChanArray[i] = make(chan updownLine)
	}

	for i, q := range queries {
		switch {
		case pushDistance > 0:
			go findUpDownCatchmentPushDistance(q, ignore, sizeArray, pushDistance, threshpair, QChanArray[i], cOut)
		default:
			go findUpDownCatchment(q, ignore, sizeArray, nofill, distArray, threshpair, QChanArray[i], cOut)
		}

	}

	for udL := range cIn {
		if udL.ambCount > threshtarg {
			continue
		}
		for i, _ := range QChanArray {
			QChanArray[i] <- udL
		}
	}

	for i, _ := range QChanArray {
		close(QChanArray[i])
	}

	cSplitDone <- true
}

// allGreaterThanEqualTo4 returns true/false; a[i] is greater that or equal to b[i]; for every item in a
// it is used to check whether any balancing is required given the --size options when some bins fall
// short of their desired sizes
func allGreaterThanEqualTo4(a, b [4]int) bool {
	for i, _ := range a {
		if a[i] < b[i] {
			return false
		}
	}
	return true
}

// given the input options, and the final size of the neighbourhoods,
// balance returns the size that the output data should be
func balance(sizetotal int, sizeIdeal [4]int, sizeObserved [4]int, nofill bool) [4]int {
	var size [4]int

	allGreaterThan := allGreaterThanEqualTo4(sizeObserved, sizeIdeal)

	if allGreaterThan {
		return sizeIdeal
	} else {
		for i := range sizeObserved {
			if sizeObserved[i] >= sizeIdeal[i] {
				size[i] = sizeIdeal[i]
			} else {
				size[i] = sizeObserved[i]
			}
		}
		if nofill {
			return size
		}
		var sizeAvail [4]int
		for i := range sizeObserved {
			if sizeObserved[i] > sizeIdeal[i] {
				sizeAvail[i] = sizeObserved[i] - sizeIdeal[i]
			} else {
				sizeAvail[i] = 0
			}
		}
		for n := 1; n > 0; {
			for i := range sizeObserved {
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

// writeUpDownCatchment writes the catchments for each query to a file/stdout
func writeUpDownCatchment(w io.Writer, results []updownCatchmentStruct) error {

	var err error
	temp := make([]string, 0)

	_, err = w.Write([]byte("query,closestsame,closestup,closestdown,closestside\n"))
	if err != nil {
		return err
	}

	for _, result := range results {

		_, err = w.Write([]byte(result.qname + ","))
		if err != nil {
			return err
		}

		temp = make([]string, 0)
		for _, hit := range result.same.catchment {
			temp = append(temp, hit.tname)
		}
		_, err = w.Write([]byte(strings.Join(temp, ";") + ","))
		if err != nil {
			return err
		}

		temp = make([]string, 0)
		for _, hit := range result.up.catchment {
			temp = append(temp, hit.tname)
		}
		_, err = w.Write([]byte(strings.Join(temp, ";") + ","))
		if err != nil {
			return err
		}

		temp = make([]string, 0)
		for _, hit := range result.down.catchment {
			temp = append(temp, hit.tname)
		}
		_, err = w.Write([]byte(strings.Join(temp, ";") + ","))
		if err != nil {
			return err
		}

		temp = make([]string, 0)
		for _, hit := range result.side.catchment {
			temp = append(temp, hit.tname)
		}
		_, err = w.Write([]byte(strings.Join(temp, ";") + "\n"))
		if err != nil {
			return err
		}
	}

	return nil
}

// writeUpdownTable writes the output in table format, including SNP-distances
func writeUpdownTable(w io.Writer, results []updownCatchmentStruct) error {
	w.Write([]byte("query,direction,distance,target\n"))
	for _, result := range results {
		for _, neighbour := range result.same.catchment {
			_, err := w.Write([]byte(strings.Join([]string{result.qname, "same", strconv.Itoa(neighbour.distance), neighbour.tname}, ",") + "\n"))
			if err != nil {
				return err
			}
		}
		for _, neighbour := range result.up.catchment {
			_, err := w.Write([]byte(strings.Join([]string{result.qname, "up", strconv.Itoa(neighbour.distance), neighbour.tname}, ",") + "\n"))
			if err != nil {
				return err
			}
		}
		for _, neighbour := range result.down.catchment {
			_, err := w.Write([]byte(strings.Join([]string{result.qname, "down", strconv.Itoa(neighbour.distance), neighbour.tname}, ",") + "\n"))
			if err != nil {
				return err
			}
		}
		for _, neighbour := range result.side.catchment {
			_, err := w.Write([]byte(strings.Join([]string{result.qname, "side", strconv.Itoa(neighbour.distance), neighbour.tname}, ",") + "\n"))
			if err != nil {
				return err
			}
		}
	}

	return nil
}

// true false every item in s == 0
func allZero(s []int) bool {
	for _, n := range s {
		if n != 0 {
			return false
		}
	}
	return true
}

// checkArgs sanity checks the command line options
func checkArgs(sizetotal int, sizeup int, sizedown int, sizeside int, sizesame int,
	distall int, distup int, distdown int, distside int) ([4]int, [4]int, error) {

	if sizetotal == 0 && allZero([]int{sizeup, sizedown, sizeside, sizesame}) && allZero([]int{distup, distdown, distside, distall}) {
		return [4]int{}, [4]int{}, errors.New(`
Please provide values to either --size-total,
or all or some of --size-up, --size-down, --size-side and --size-same,

                            * or *

to --dist-all,
or all or some of --dist-up, --dist-down, --dist-side,

and optionally to the --size- arguments (to constrain the output by both distance and sample size)`)
	}

	if sizetotal != 0 && !allZero([]int{sizeup, sizedown, sizeside, sizesame}) {
		return [4]int{}, [4]int{}, errors.New("warning: setting --size-total overrides --size-up, --size-down, --size-side and --size-same")
	}

	// bucket sizes for same,up,down,side,total respectively
	var sizeArray [4]int

	switch {
	case sizetotal > 0: // size total is set - try to split evenly
		sizeArray[1] = int(math.Floor(float64(sizetotal) / float64(4)))
		sizeArray[2] = int(math.Floor(float64(sizetotal) / float64(4)))
		sizeArray[3] = int(math.Floor(float64(sizetotal) / float64(4)))
		sizeArray[0] = sizetotal - (sizeArray[1] + sizeArray[2] + sizeArray[3])

	case !allZero([]int{sizeup, sizedown, sizeside, sizesame}): // bucket sizes are set, use these:
		sizeArray[0] = sizesame
		sizeArray[1] = sizeup
		sizeArray[2] = sizedown
		sizeArray[3] = sizeside

		// an easter egg - set any of these to -1 on the command line and you will get all of the seqs in that bucket
		for i, n := range []int{sizesame, sizeup, sizedown, sizeside} {
			if n == -1 {
				sizeArray[i] = math.MaxInt32
			}
		}

	default:
		// get everything (but in order, eventually, because it will be sorted, or a filter on distance measure will be applied).
		sizeArray[0] = math.MaxInt32
		sizeArray[1] = math.MaxInt32
		sizeArray[2] = math.MaxInt32
		sizeArray[3] = math.MaxInt32
	}

	if distall != 0 && !allZero([]int{distup, distdown, distside}) {
		return [4]int{}, [4]int{}, errors.New("warning: setting --dist-all overrides --dist-up, --dist-down and --dist-side")
	}

	var distArray [4]int

	switch {
	case distall > 0:
		distArray[0] = 0
		distArray[1] = distall
		distArray[2] = distall
		distArray[3] = distall
	case !allZero([]int{distup, distdown, distside}):
		distArray[0] = 0
		distArray[1] = distup
		distArray[2] = distdown
		distArray[3] = distside
	default:
		distArray[0] = math.MaxInt32
		distArray[1] = math.MaxInt32
		distArray[2] = math.MaxInt32
		distArray[3] = math.MaxInt32
	}

	return sizeArray, distArray, nil
}

// the sum of all the integers in a
func sum4(a [4]int) int {
	var t int
	for _, n := range a {
		t += n
	}
	return t
}

// TopRanking finds pseudo-tree-aware catchments for query sequences, given a large database of target sequences, the closest
// of which should be returned in the output. Targets are split into bins depending on whether they are likely direct ancestors of,
// direct descendants of, polyphyletic with, or exactly the same as, the query
func TopRanking(query, target, reference io.Reader, out io.Writer, table bool,
	q_in_type, t_in_type string, ignoreArray []string,
	sizetotal int, sizeup int, sizedown int, sizeside int, sizesame int,
	distall int, distup int, distdown int, distside int,
	threshpair float32, threshtarg int, nofill bool, pushdist int) error {

	sizeArray, distArray, err := checkArgs(sizetotal, sizeup, sizedown, sizeside, sizesame, distall, distup, distdown, distside)
	if err != nil {
		switch {
		case err.Error() == "warning: setting --size-total overrides --size-up, --size-down, --size-side and --size-same":
			os.Stderr.WriteString(err.Error() + "\n")
		case err.Error() == "warning: setting --dist-all overrides --dist-up, --dist-down and --dist-side":
			os.Stderr.WriteString(err.Error() + "\n")
		default:
			return err
		}
	}

	var refSeq []byte
	if q_in_type == "fasta" || t_in_type == "fasta" {
		temp, err := fastaio.ReadEncodeAlignmentToList(reference, false)
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
		queries, err = readCSVToUDLList(query)
		if err != nil {
			return err
		}
	case "fasta":
		queries, err = fastaToUDLList(query, refSeq)
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
		go readCSVToUDLChan(target, cudL, cErr, cReadDone)
	case "fasta":
		go readFastaToUDLChan(target, refSeq, cudL, cErr, cReadDone)
	}

	go splitInput(queries, ignoreArray,
		sizeArray, nofill, distArray, threshpair, threshtarg, pushdist,
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

	if table {
		err = writeUpdownTable(out, QResultsArray)
	} else {
		err = writeUpDownCatchment(out, QResultsArray)
	}
	if err != nil {
		return err
	}

	return nil
}
