package sam

import (
	"fmt"
	"sync"
	"sort"
	"os"
	"path"
	"strings"
	"strconv"

	"github.com/cov-ert/gofasta/pkg/genbank"
	"github.com/cov-ert/gofasta/pkg/fastaio"

	biogosam "github.com/biogo/hts/sam"
)


// alignPair is a struct for storing an aligned ref + query sequence pair
type alignPair struct {
	ref []byte
	query []byte
	refname string
	queryname string
	descriptor string
	featPosArray []int // parsed start/end positions of feature from annotation (can be >length(2) if Join())
	featType string
	featName string
	idx int // for retaining input order in the output
}

// for passing groups of alignPairs around with an index which is used to retain input
// order in the output
type alignPairs struct {
	aps []alignPair
	idx int
}

// store some information about a single insertion from cigars in a block of SAM records
type insertionOccurence struct {
	start int
	length int
	rowNumber int
}

type byStart []insertionOccurence

func (a byStart) Len() int           { return len(a) }
func (a byStart) Less(i, j int) bool { return a[i].start < a[j].start }
func (a byStart) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }


type alignedBlockInfo struct {
	seqpairArray []alignPair // an array of ref + query seqs from primary + secondary mappings
	cigarArray []biogosam.Cigar // the cigars for each mapping in seqpairArray
	posArray []int // the start POS from the sam file for each mapping in seqpairArray
}

func blockToSeqPair(alignedBlock alignedBlockInfo, ref []byte) alignPair {
	// get all insertions in order of occurence - earliest first
	//
	// and then collapse the sequences in the block down to one sequence
	//
	// The incoming sequences are all left-aligned. So for every insertion
	// in every cigar in cigarArray, need to gap all the OTHER imcoming sequences
	// except the one that it occurs in (because this has already been done)
	//
	// need to do this in order, left-most insertion first.

	insertions := make([]insertionOccurence, 0)

	for i, cigar := range(alignedBlock.cigarArray) {

		pos := alignedBlock.posArray[i]

		for _, op := range(cigar) {
			// populate orderedInsertions here...

			// fmt.Println(op.Type())

			if op.Type().String() == "I" {

				I := insertionOccurence{start: pos, length: op.Len(), rowNumber: i}

				insertions = append(insertions, I)
			}

			// if op.Type().Consumes().Query == 1 || op.Type().Consumes().Reference == 1 {
			if op.Type().Consumes().Reference == 1 {
				pos += op.Len()
			}
		}
	}

	refSeqArray := make([][]byte, len(alignedBlock.seqpairArray))
	queSeqArray := make([][]byte, len(alignedBlock.seqpairArray))

	if len(insertions) > 0 {
		sort.Sort(byStart(insertions))

		for _, insertion := range(insertions) {
			rowNumber := insertion.rowNumber
			for j, seqPair := range(alignedBlock.seqpairArray) {
				// don't reinsert - the insertion already exists in this one
				if j == rowNumber {
					refSeqArray[j] = seqPair.ref
					queSeqArray[j] = seqPair.query
					continue
				}

				// some lines will never get this far
				if insertion.start > len(seqPair.ref) {
					refSeqArray[j] = seqPair.ref
					queSeqArray[j] = seqPair.query
					continue
				}

				gaps := make([]byte, insertion.length)
				for i, _ := range gaps {
					gaps[i] = '-'
				}

				refSeqArray[j] = seqPair.ref[:insertion.start]
				refSeqArray[j] = append(refSeqArray[j], gaps...)
				refSeqArray[j] = append(refSeqArray[j], seqPair.ref[insertion.start:]...)

				queSeqArray[j] = seqPair.query[:insertion.start]
				queSeqArray[j] = append(queSeqArray[j], gaps...)
				queSeqArray[j] = append(queSeqArray[j], seqPair.query[insertion.start:]...)
			}
		}

	} else {
		for i, seqPair := range(alignedBlock.seqpairArray) {
			refSeqArray[i] = seqPair.ref
			queSeqArray[i] = seqPair.query
		}
	}

	max := 0

	for _, line := range(refSeqArray){
		if len(line) > max {
			max = len(line)
		}
	}

	RBlock := make([][]byte, 0)
	for _, line := range(refSeqArray){
		if len(line) < max {
			diff := max - len(line)
			stars := make([]byte, diff)
			for i, _ := range stars {
				stars[i] = '*'
			}
			line = append(line, stars...)
		}
		RBlock = append(RBlock, line)
		// fmt.Println(line)
	}

	R := checkAndGetFlattenedSeq(RBlock, "reference")

	QBlock := make([][]byte, 0)
	for _, line := range(queSeqArray){
		if len(line) < max {
			diff := max - len(line)
			stars := make([]byte, diff)
			for i, _ := range stars {
				stars[i] = '*'
			}
			line = append(line, stars...)
		}
		QBlock = append(QBlock, line)
	}

	Q := checkAndGetFlattenedSeq(QBlock, alignedBlock.seqpairArray[0].queryname)

	// extend the alignment to the ref length + the total number of
	// insertions relative to the reference...
	totalInsertionLength := 0
	for _, I := range(insertions) {
		totalInsertionLength += I.length
	}

	if len(R) < totalInsertionLength + len(ref) {
		diff := (totalInsertionLength + len(ref)) - len(R)
		extendRight := make([]byte, diff)
		for i, _ := range(extendRight) {
			extendRight[i] = '*'
		}
		Q = append(Q, extendRight...)
		R = append(R, ref[len(ref) - diff:]...)
	}

	Q = swapInNs(Q)

	return alignPair{ref: R, query: Q}
}

// blockToPairwiseAlignment should convert a block of SAM records that correspond
// to the same query sequence to a pairwise alignment between that query and the
// reference. It should return the pair of sequences (query + reference) aligned
// to each other - insertions in the query can be represented or not.
func blockToPairwiseAlignment(cSR chan samRecords, cPair chan alignPair, cErr chan error, ref []byte, omitIns bool) {

	for group := range(cSR) {

		// seqs is an array of seqs, one item for each line in the
		// block of sam lines for one query
		seqs := make([]alignPair, 0)
		cigars := make([]biogosam.Cigar, 0)
		positions := make([]int, 0)

		infoBlock := alignedBlockInfo{}

		// might not use this:
		Q := make([][]byte, 0)


		if !omitIns {
			// populate it
			for _, line := range(group.records) {
				seq, refseq, err := getOneLinePlusRef(line, ref, !omitIns)
				if err != nil {
					cErr<- err
				}
				seqs = append(seqs, alignPair{ref: refseq, query: seq, queryname: line.Name})
				cigars = append(cigars, line.Cigar)
				positions = append(positions, line.Pos)
			}

			infoBlock.seqpairArray = seqs
			infoBlock.cigarArray = cigars
			infoBlock.posArray = positions

			pair := blockToSeqPair(infoBlock, ref)
			pair.refname = string(group.records[0].Ref.Name())
			pair.queryname = group.records[0].Name
			pair.idx = group.idx
			cPair<- pair

		} else {
			qname := group.records[0].Name

			for _, line := range(group.records) {
				seq, _, err := getOneLinePlusRef(line, ref, !omitIns)
				if err != nil {
					cErr<- err
				}
				Q = append(Q, seq)
			}

			Qflat := swapInNs(checkAndGetFlattenedSeq(Q, qname))

			pair := alignPair{ref: ref, query: Qflat}
			pair.queryname = group.records[0].Name
			pair.refname = string(group.records[0].Ref.Name())
			pair.idx = group.idx

			cPair<- pair
		}
	}

	return
}

func parsePositions(position string) ([]int, error) {
	var A []int
	if position[0:4] == "join" {
		A = make([]int, 0)
		position = strings.TrimLeft(position, "join(")
		position = strings.TrimRight(position, ")")
		ranges := strings.Split(position, ",")
		for _, x := range(ranges) {
			y := strings.Split(x, "..")
			for _, z := range(y) {
				temp, err := strconv.Atoi(z)
				if err != nil {
					return []int{}, err
				}
				A = append(A, temp)
			}
		}
	} else {
		A = make([]int, 0)
		y := strings.Split(position, "..")
		for _, z := range(y) {
			temp, err := strconv.Atoi(z)
			if err != nil {
				return []int{}, err
			}
			A = append(A, temp)
		}
	}

	return A, nil
}

func getFeaturesFromAnnotation(gb genbank.Genbank, annotation string) []genbank.GenbankFeature {

	FEATS := make([]genbank.GenbankFeature, 0)

	for _, F := range(gb.FEATURES) {
		if F.Feature == annotation {
			FEATS = append(FEATS, F)
		}
	}

	return FEATS
}

func getRefAdjustedPositions(seq []byte) []int {
	idx := make([]int, len(seq))
	pos := 0
	for i, nuc := range(seq) {
		if nuc != '-' {
			pos +=1
		}
		idx[i] = pos
	}

	return idx
}

func findOffsetPos(i int, a []int) int {
	var pos int
	for j, x := range(a) {
		if x == i {
			pos = j
			break
		}
	}
	return pos
}

// TODO - allow multiple feature types in features []genbank.GenbankFeature
func parseAlignmentByAnnotation(features []genbank.GenbankFeature, cPairIn chan alignPair, cPairOut chan alignPairs, cErr chan error) {

	// if no feature is specified on the command line:
	if len(features) == 0 {
		for pair := range(cPairIn) {
			pair.descriptor = pair.queryname
			cPairOut<- alignPairs{aps: []alignPair{pair}, idx: pair.idx}
		}
	// if one feature is specified on the command line:
	} else {
		ftype := features[0].Feature

		var anno string

		if ftype == "CDS" || ftype == "gene" {
			anno = "gene"
		}

		if ftype == "source" {
			anno = "organism"
		}

		for pair := range(cPairIn) {

			A := alignPairs{idx: pair.idx}

			idx := getRefAdjustedPositions(pair.ref)

			for _, feature := range(features) {

				subPair := alignPair{}

				subPair.refname = pair.refname
				subPair.queryname = pair.queryname
				subPair.featType = feature.Feature
				subPair.featName = feature.Info[anno]
				subPair.descriptor = pair.queryname + "." + feature.Feature + "." + strings.ReplaceAll(feature.Info[anno], " ", "_")

				positions, err := parsePositions(feature.Pos)
				if err != nil {
					cErr<- err
				}

				subPair.featPosArray = positions

				var newRef []byte
				var newQue []byte

				if len(positions) / 2 > 1 {
					for i := 0; i < len(positions); i += 2 {
						start := findOffsetPos(positions[i], idx)
						stop := findOffsetPos(positions[i + 1], idx) + 1
						newRef = append(newRef, pair.ref[start:stop]...)
						newQue = append(newQue, pair.query[start:stop]...)
					}
					subPair.ref = newRef
					subPair.query = newQue
				} else {
					start := findOffsetPos(positions[0], idx)
					stop := findOffsetPos(positions[1], idx) + 1
					newRef = pair.ref[start:stop]
					newQue = pair.query[start:stop]
					subPair.ref = newRef
					subPair.query = newQue
				}

				A.aps = append(A.aps, subPair)
			}

			cPairOut<- A
		}
	}

	return
}

func writePairwiseAlignment(p string, cPair chan alignPairs, cWriteDone chan bool, cErr chan error, omitRef bool) {

	_ = path.Join()

	var err error

	if p == "stdout" {
		for array := range cPair {
			for _, AP := range(array.aps) {
				if ! omitRef {
					_, err = fmt.Fprintln(os.Stdout, ">" + AP.refname)
					if err != nil {
						cErr <- err
					}
					_, err = fmt.Fprintln(os.Stdout, string(AP.ref))
					if err != nil {
						cErr <- err
					}
				}
				_, err = fmt.Fprintln(os.Stdout, ">" + AP.queryname)
				if err != nil {
					cErr <- err
				}
				_, err = fmt.Fprintln(os.Stdout, string(AP.query))
				if err != nil {
					cErr <- err
				}
			}
		}
	} else {
		os.MkdirAll(p, 0755)

		for array := range cPair {
			for _, AP := range(array.aps) {
				des := strings.ReplaceAll(AP.descriptor, "/", "_")
				des = strings.ReplaceAll(des, "|", "_")
				f, err := os.Create(path.Join(p, des + ".fasta"))
				if err != nil {
					cErr <- err
				}
				if ! omitRef {
					_, err = f.WriteString(">" + AP.refname + "\n")
					if err != nil {
						cErr <- err
					}
					_, err = f.WriteString(string(AP.ref) + "\n")
					if err != nil {
						cErr <- err
					}
				}
				_, err = f.WriteString(">" + AP.queryname + "\n")
				if err != nil {
					cErr <- err
				}
				_, err = f.WriteString(string(AP.query) + "\n")
				if err != nil {
					cErr <- err
				}
				f.Close()
			}
		}
	}
	cWriteDone<- true
}

// ToPairAlign converts a SAM file into pairwise fasta-format alignments
// optionally including the reference, optionally split by annotations,
// optionally skipping insertions relative to the reference
func ToPairAlign(samFile string, referenceFile string, genbankFile string, feat string, outpath string, omitRef bool, omitIns bool, threads int) error {

	gb, err := genbank.ReadGenBank(genbankFile)
	if err != nil {
		return err
	}

	// NB probably uncomment the below and use it for checks (e.g. for
	// reference length)
	// samHeader, err := getSamHeader(samFile)
	// if err != nil {
	// 	return err
	// }

	// refLen := samHeader.Refs()[0].Len()

	cErr := make(chan error)

	cRef := make(chan fastaio.FastaRecord)
	cRefDone := make(chan bool)

	go fastaio.ReadAlignment(referenceFile, cRef, cErr, cRefDone)

	var refSeq string
	// var refName string

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case FR := <-cRef:
			refSeq = FR.Seq
			// refName = FR.ID
		case <-cRefDone:
			close(cRef)
			n--
		}
	}

	// fmt.Println(refSeq)
	// fmt.Println()

	cSR := make(chan samRecords, threads)
	cSH := make(chan biogosam.Header)

	cPairAlign := make(chan alignPair)
	cPairParse := make(chan alignPairs)

	cReadDone := make(chan bool)
	cAlignWaitGroupDone := make(chan bool)
	cParseWaitGroupDone := make(chan bool)
	cWriteDone := make(chan bool)

	go groupSamRecords(samFile, cSH, cSR, cReadDone, cErr)

	_ = <-cSH

	go writePairwiseAlignment(outpath, cPairParse, cWriteDone, cErr, omitRef)

	var wgAlign sync.WaitGroup
	wgAlign.Add(threads)

	var wgParse sync.WaitGroup
	wgParse.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			blockToPairwiseAlignment(cSR, cPairAlign, cErr, []byte(refSeq), omitIns)
			wgAlign.Done()
		}()
	}

	features := getFeaturesFromAnnotation(gb, feat)

	for n := 0; n < threads; n++ {
		go func() {
			parseAlignmentByAnnotation(features, cPairAlign, cPairParse, cErr)
			wgParse.Done()
		}()
	}

	go func() {
		wgAlign.Wait()
		cAlignWaitGroupDone<- true
	}()

	go func() {
		wgParse.Wait()
		cParseWaitGroupDone<- true
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
		case <-cAlignWaitGroupDone:
			close(cPairAlign)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cParseWaitGroupDone:
			close(cPairParse)
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

	// fmt.Println(gb)
	// fmt.Println(samHeader)

	return nil
}
