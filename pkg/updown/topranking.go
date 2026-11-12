/*
Package updown implements functions that calculate pseudo-tree aware
snp-distances between sequences. It relies on there being a suitable
sequence that one can use to represent the "root" (ancestor) of a
hypothetical phylogenetic tree, and for the distance between relevant
sequences to be small enough that multiple hits are unlikely
*/
package updown

import (
	"context"
	"errors"
	"io"
	"math"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
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

// isSiteAmb asks if a position is within an ambiguity tract, given an array of pairs of start-stop coordinates
// of such tracts. Co-ordinates are 1-based, inclusive
func isSiteAmb(pos int, a []int) bool {
	for i := 0; i < len(a); i += 2 {
		if pos >= a[i] && pos <= a[i+1] {
			return true
		}
	}
	return false
}

// is pos present in list, true/false
// sort.Search* is the go standard library's implementation of binary search.
// Should scale O(log(N))
func posOverlapBinarySearch(list []int, pos int) bool {

	idx := sort.SearchInts(list, pos)

	if idx < len(list) && list[idx] == pos {
		return true
	}

	return false
}

// is snp present in list, true/false
func snpOverlapBinarySearch(list []string, snp string) bool {

	idx := sort.SearchStrings(list, snp)

	if idx < len(list) && list[idx] == snp {
		return true
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
		// is the site ambiguous in the target
		if isSiteAmb(q.snpsPos[i], t.ambs) {
			table[3]++
			// if this exact snp (allele and position) is present in the target
		} else if snpOverlapBinarySearch(t.snpsSorted, qsnp) {
			table[1]++
		} else {
			// otherwise we +1 the Q snps bin and append the snp to d
			table[0]++
			d = append(d, q.snpsPos[i])
		}
	}
	d_plus := 0
	for i, tsnp := range t.snps {
		// is the site ambiguous in the query
		if isSiteAmb(t.snpsPos[i], q.ambs) {
			table[3]++
			// if there is no overlap between the exact snp (allele and position) and the query's snps
		} else if !snpOverlapBinarySearch(q.snpsSorted, tsnp) {
			// we +1 the T snps bin
			table[2]++
			// if this position has not already been added to the list of differences we use to calculate d, then add it here too (otherwise two mutations historically at the same position, but only counts 1 for distance)
			if !posOverlapBinarySearch(d, t.snpsPos[i]) {
				d_plus++
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

	distance := len(d) + d_plus

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
// NB - could use binary search here but this function is only ever called to check whether a sequence ID is
// present in a list of IDs to ignore, which is likely always going to be small (?!), so slice iteration seems reasonable?
func stringInArray(s string, sa []string) bool {
	for i := range sa {
		if s == sa[i] {
			return true
		}
	}
	return false
}

// byDistance sorts a catchment by snp distance, then by number of ambiguities, keeping ties in their order
func byDistance(c []resultsStruct) {
	sort.SliceStable(c, func(i, j int) bool {
		return c[i].distance < c[j].distance || (c[i].distance == c[j].distance && c[i].ambCount < c[j].ambCount)
	})
}

// consider adds rs to the catchment if it has fewer than size members, or if rs is closer than the
// furthest member (or as close, with fewer ambiguities)
func (nS *updownCatchmentSubStruct) consider(rs resultsStruct, size int) {
	switch {
	case len(nS.catchment) < size:
		nS.catchment = append(nS.catchment, rs)
		if len(nS.catchment) == size {
			rearrangeCatchment(nS, size)
		}
	case rs.distance < nS.maxDist,
		rs.distance == nS.maxDist && rs.ambCount < nS.minAmbig:
		nS.catchment = append(nS.catchment, rs)
		rearrangeCatchment(nS, size)
	}
}

// a search collects the neighbours of one query, one target at a time
type search interface {
	// add considers a target that is in direction (0,1,2,3 = same,up,down,side) at snp distance
	add(target updownLine, direction, distance int)
	// finish returns the neighbours
	finish() updownCatchmentStruct
}

// catchmentSearch assigns targets to the bins by direction, keeping the closest in each, within the
// distance limits for each bin. It can't do any pushing if any bins are empty after all targets have
// been processed.
type catchmentSearch struct {
	neighbours updownCatchmentStruct
	sizeArray  [4]int
	sizetotal  int // the maximum size of any bin (for filling in)
	distArray  [4]int
	nofill     bool
}

func newCatchmentSearch(q updownLine, sizeArray [4]int, nofill bool, distArray [4]int) *catchmentSearch {
	s := &catchmentSearch{sizeArray: sizeArray, distArray: distArray, nofill: nofill}
	s.neighbours = updownCatchmentStruct{qname: q.id, qidx: q.idx}
	s.neighbours.same = updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)}
	s.neighbours.up = updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)}
	s.neighbours.down = updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)}
	s.neighbours.side = updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)}
	s.sizetotal = sum4(sizeArray)
	for _, n := range sizeArray {
		if n == math.MaxInt32 {
			s.sizetotal = math.MaxInt32
		}
	}
	return s
}

func (s *catchmentSearch) bin(direction int) *updownCatchmentSubStruct {
	return [4]*updownCatchmentSubStruct{&s.neighbours.same, &s.neighbours.up, &s.neighbours.down, &s.neighbours.side}[direction]
}

func (s *catchmentSearch) add(target updownLine, direction, distance int) {
	// if the distance is too big, skip it
	if distance > s.distArray[direction] {
		return
	}
	s.bin(direction).consider(resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}, s.sizetotal)
}

func (s *catchmentSearch) finish() updownCatchmentStruct {
	var sizeObserved [4]int
	for direction := range sizeObserved {
		b := s.bin(direction)
		// a bin that never filled up has not been sorted yet
		if n := len(b.catchment); n < s.sizetotal && n > 0 {
			rearrangeCatchment(b, n)
		}
		sizeObserved[direction] = len(b.catchment)
	}

	size := balance(s.sizetotal, s.sizeArray, sizeObserved, s.nofill)
	for direction := range size {
		b := s.bin(direction)
		b.catchment = b.catchment[0:size[direction]]
	}
	return s.neighbours
}

// pushSearch assigns targets to the bins by direction, keeping, for the up, down and side bins, all the
// targets at the pushDist smallest snp distances found in that direction (so that a bin is not empty
// just because nothing is within a fixed distance).
type pushSearch struct {
	q              updownLine
	pushDist       int
	same           updownCatchmentSubStruct // the sequences on the polytomy
	up, down, side pushCatchmentSubStruct
}

func newPushSearch(q updownLine, pushDist int) *pushSearch {
	return &pushSearch{
		q:        q,
		pushDist: pushDist,
		same:     updownCatchmentSubStruct{catchment: make([]resultsStruct, 0)},
		up:       pushCatchmentSubStruct{catchmentMap: make(map[int][]resultsStruct)},
		down:     pushCatchmentSubStruct{catchmentMap: make(map[int][]resultsStruct)},
		side:     pushCatchmentSubStruct{catchmentMap: make(map[int][]resultsStruct)},
	}
}

func (s *pushSearch) add(target updownLine, direction, distance int) {
	rs := resultsStruct{tname: target.id, ambCount: target.ambCount, distance: distance}
	if direction == 0 {
		s.same.catchment = append(s.same.catchment, rs)
		return
	}
	pC := [4]*pushCatchmentSubStruct{nil, &s.up, &s.down, &s.side}[direction]
	// is the distance lower, or are there not enough distances yet
	if distance <= pC.maxDist || pC.nDists < s.pushDist {
		refactorPushCatchment(pC, rs, s.pushDist)
	}
}

func (s *pushSearch) finish() updownCatchmentStruct {
	neighbours := updownCatchmentStruct{qname: s.q.id, qidx: s.q.idx, same: s.same}
	neighbours.up = pushCatchment2Catchment(s.up)
	byDistance(neighbours.up.catchment)
	neighbours.down = pushCatchment2Catchment(s.down)
	byDistance(neighbours.down.catchment)
	neighbours.side = pushCatchment2Catchment(s.side)
	byDistance(neighbours.side.catchment)
	// TO DO - balance the neighbours here if sizeArray is given? Would need a balance function specific to pushing
	return neighbours
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
	distall int, distup int, distdown int, distside int, distpush int) ([4]int, [4]int, error) {

	for _, v := range []int{sizetotal, sizeup, sizedown, sizeside, sizesame, distall, distup, distdown, distside, distpush} {
		if v < 0 {
			return [4]int{}, [4]int{}, errors.New("--size-* and --dist-* values must not be negative")
		}
	}

	if sizetotal == 0 && allZero([]int{sizeup, sizedown, sizeside, sizesame}) && distpush == 0 && allZero([]int{distup, distdown, distside, distall}) {
		return [4]int{}, [4]int{}, errors.New(`
Please provide values to either --size-total,
or all or some of --size-up, --size-down, --size-side and --size-same,

                            * or *

to --dist-all,
or all or some of --dist-up, --dist-down, --dist-side
or to --dist-push`)
	}

	if (sizetotal != 0 && !allZero([]int{sizeup, sizedown, sizeside, sizesame})) && !allZero([]int{distup, distdown, distside, distall}) && (distpush > 0) {
		return [4]int{}, [4]int{}, errors.New("you can't combine --dist* and --size* flags if you also invoke --dist-push")
	}

	if sizetotal != 0 && !allZero([]int{sizeup, sizedown, sizeside, sizesame}) {
		os.Stderr.WriteString("warning: setting --size-total overrides --size-up, --size-down, --size-side and --size-same")
	}

	if distall != 0 && !allZero([]int{distup, distdown, distside}) {
		os.Stderr.WriteString("warning: setting --dist-all overrides --dist-up, --dist-down and --dist-side")
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

// returns the sum of all the integers in an length-4 array of ints
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
	threshpair float32, threshtarg int, nofill bool, distpush int, strict bool) error {

	sizeArray, distArray, err := checkArgs(sizetotal, sizeup, sizedown, sizeside, sizesame, distall, distup, distdown, distside, distpush)
	if err != nil {
		return err
	}

	var refSeq []byte
	if q_in_type == "fasta" || t_in_type == "fasta" {
		temp, err := fasta.LoadEncodeAlignment(reference, false, false, false)
		if err != nil {
			return err
		}
		if len(temp) > 1 {
			return errors.New("More than one record in --reference")
		}
		refSeq = temp[0].Seq
	}

	var queries []updownLine

	switch q_in_type {
	case "csv":
		queries, err = readCSVToUDLList(query, strict)
		if err != nil {
			return err
		}
	case "fasta":
		queries, err = fastaToUDLList(query, refSeq, strict)
		if err != nil {
			return err
		}
	}

	searches := make([]search, len(queries))
	for i, q := range queries {
		if distpush > 0 {
			searches[i] = newPushSearch(q, distpush)
		} else {
			searches[i] = newCatchmentSearch(q, sizeArray, nofill, distArray)
		}
	}

	// a target is read as a fasta record or as a line of updown list's csv output
	type targetIn struct {
		record *fasta.Record
		row    csvRow
	}
	// a target with its direction and snp distance from each query
	type rankedTarget struct {
		line                 updownLine
		skip                 bool
		directions, distance []int
	}

	var next func() (targetIn, error)
	switch t_in_type {
	case "csv":
		src := newCSVSource(target, errors.New("bad header when parsing --target csv: is this the output of gofasta updown list?"))
		next = func() (targetIn, error) {
			row, err := src()
			return targetIn{row: row}, err
		}
	case "fasta":
		src := fasta.NewRecordSource(target)
		next = func() (targetIn, error) {
			record, err := src()
			return targetIn{record: &record}, err
		}
	default:
		return errors.New("unknown --target file type: " + t_in_type)
	}

	DA := encoding.MakeDecodingArray()
	work := func(in targetIn) (rankedTarget, error) {
		var t updownLine
		var err error
		if in.record != nil {
			t, err = lineFromRecord(refSeq, *in.record, DA)
		} else {
			t, err = parseCSVRow(in.row)
		}
		if err != nil {
			return rankedTarget{}, err
		}
		if t.ambCount > threshtarg || (len(ignoreArray) > 0 && stringInArray(t.id, ignoreArray)) {
			return rankedTarget{skip: true}, nil
		}
		rt := rankedTarget{line: t, directions: make([]int, len(queries)), distance: make([]int, len(queries))}
		for i, q := range queries {
			// direction values 0,1,2,3 = same,up,down,side; distance is -1 if the pair fails the ambiguity threshold
			rt.directions[i], rt.distance[i] = whichWay(q, t, threshpair)
		}
		return rt, nil
	}

	h := pipeline.NewHandler(strict, os.Stderr, pipeline.Noun{One: "target", Many: "targets"}, func(rt rankedTarget) error {
		if rt.skip {
			return nil
		}
		for i := range searches {
			if rt.distance[i] >= 0 {
				searches[i].add(rt.line, rt.directions[i], rt.distance[i])
			}
		}
		return nil
	})
	if err := pipeline.Run(context.Background(), runtime.NumCPU(), next, work, h.Sink); err != nil {
		return err
	}
	if err := h.Finish(); err != nil {
		return err
	}

	QResultsArray := make([]updownCatchmentStruct, len(queries))
	for i := range searches {
		QResultsArray[i] = searches[i].finish()
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
