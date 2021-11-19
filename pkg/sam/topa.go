package sam

import (
	"fmt"
	"io"
	"os"
	"path"
	"sort"
	"strings"
	"sync"

	"github.com/cov-ert/gofasta/pkg/fastaio"
	"github.com/cov-ert/gofasta/pkg/variants"

	biogosam "github.com/biogo/hts/sam"
)

// alignPair is a struct for storing an aligned reference & query sequence pair
type alignPair struct {
	ref          []byte
	query        []byte
	refname      string
	queryname    string
	descriptor   string
	featPosArray []int // parsed start/end positions of feature from annotation (can be >length(2) if Join())
	featType     string
	featName     string
	idx          int // for retaining input order in the output
}

// alignPairs is for passing groups of alignPair around with an index which is used to retain input
// order in the output
type alignPairs struct {
	aps []alignPair
	idx int
}

// insertionOccurence stores information about a single insertion from cigars in a block of SAM records,
// which is used to appropriate gap the other sam lines belonging to the same query sequence
type insertionOccurence struct {
	start     int
	length    int
	rowNumber int
}

type byStart []insertionOccurence

// the functions required by the sorting interface so that insertion structs can be started by their genomic position
func (a byStart) Len() int           { return len(a) }
func (a byStart) Less(i, j int) bool { return a[i].start < a[j].start }
func (a byStart) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }

// alignedBlockInfo stores arrays of sam fields for one queriy's sam lines -
// they should be equal length and the item at the same index in each array should
// originate from the same sam line
type alignedBlockInfo struct {
	seqpairArray []alignPair      // an array of ref + query seqs from primary + secondary mappings
	cigarArray   []biogosam.Cigar // the cigars for each mapping in seqpairArray
	posArray     []int            // the start POS from the sam file for each mapping in seqpairArray
}

// blockToSeqPair converts a block of sam records corresponding to the same query sequence to a pairwise
// alignment between reference and query, including insertions in the query relative to the reference
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

	for i, cigar := range alignedBlock.cigarArray {

		pos := alignedBlock.posArray[i]

		for _, op := range cigar {
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

	// copy the sequences into new slices
	refSeqArray := make([][]byte, len(alignedBlock.seqpairArray))
	queSeqArray := make([][]byte, len(alignedBlock.seqpairArray))
	for i := range alignedBlock.seqpairArray {
		refSeqArray[i] = alignedBlock.seqpairArray[i].ref
		queSeqArray[i] = alignedBlock.seqpairArray[i].query
	}

	// if there are insertions, we need to modify these new slices
	if len(insertions) > 0 {
		sort.Sort(byStart(insertions))

		// if we are going to insert multiple insertions into one pair then we will need to keep track
		// of the coordinate offset after the first one
		offsets := make([]int, len(alignedBlock.seqpairArray))

		// for every insertion
		for _, insertion := range insertions {
			// this is the pair it is already present in, which we will skip:
			rowNumber := insertion.rowNumber
			for j, seqPair := range alignedBlock.seqpairArray {
				// don't reinsert - the insertion already exists in this one
				if j == rowNumber {
					continue
				}

				// if the insertions starts after the (offset) length of this sequence,
				// we don't have to do anything to this pair here
				if insertion.start > len(alignedBlock.seqpairArray[j].ref)-offsets[j] {
					continue
				}

				// otherwise, we make a slice of gaps to insert into the slices
				gaps := make([]byte, insertion.length)
				for k := range gaps {
					gaps[k] = '-'
				}

				refSeqArray[j] = refSeqArray[j][:insertion.start+offsets[j]]
				refSeqArray[j] = append(refSeqArray[j], gaps...)
				refSeqArray[j] = append(refSeqArray[j], seqPair.ref[insertion.start+offsets[j]:]...)

				queSeqArray[j] = seqPair.query[:insertion.start+offsets[j]]
				queSeqArray[j] = append(queSeqArray[j], gaps...)
				queSeqArray[j] = append(queSeqArray[j], seqPair.query[insertion.start+offsets[j]:]...)

				// and we add the relevant offset to account for this insertion in future coordinates
				offsets[j] += insertion.length
			}
		}
	}

	max := 0

	for _, line := range refSeqArray {
		if len(line) > max {
			max = len(line)
		}
	}

	RBlock := make([][]byte, 0)
	for _, line := range refSeqArray {
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
	for _, line := range queSeqArray {
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
	for _, I := range insertions {
		totalInsertionLength += I.length
	}

	if len(R) < totalInsertionLength+len(ref) {
		diff := (totalInsertionLength + len(ref)) - len(R)
		extendRight := make([]byte, diff)
		for i, _ := range extendRight {
			extendRight[i] = '*'
		}
		Q = append(Q, extendRight...)
		R = append(R, ref[len(ref)-diff:]...)
	}

	Q = swapInNs(Q)

	return alignPair{ref: R, query: Q}
}

// blockToPairwiseAlignment should convert a block of SAM records that correspond
// to the same query sequence to a pairwise alignment between that query and the
// reference. It should return the pair of sequences (query + reference) aligned
// to each other - insertions in the query can be represented or not.
func blockToPairwiseAlignment(cSR chan samRecords, cPair chan alignPair, cErr chan error, ref []byte, omitIns bool) {

	for group := range cSR {

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
			for _, line := range group.records {
				seq, refseq, err := getOneLinePlusRef(line, ref, !omitIns)
				if err != nil {
					cErr <- err
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
			cPair <- pair

		} else {
			qname := group.records[0].Name

			for _, line := range group.records {
				seq, _, err := getOneLinePlusRef(line, ref, !omitIns)
				if err != nil {
					cErr <- err
				}
				Q = append(Q, seq)
			}

			Qflat := swapInNs(checkAndGetFlattenedSeq(Q, qname))

			pair := alignPair{ref: ref, query: Qflat}
			pair.queryname = group.records[0].Name
			pair.refname = string(group.records[0].Ref.Name())
			pair.idx = group.idx

			cPair <- pair
		}
	}

	return
}

// // getOffset returns an array of integers the same length as the reference sequence, which consist of the
// // the number of bases to add needed to convert each position to msa coordinates when insertions relative
// // to the refernce are included
// func getOffset(refseq []byte) []int {
// 	// we make an array of integers to offset the positions by.
// 	// this should be the same length as a degapped refseq?
// 	degappedLen := 0
// 	for _, nuc := range refseq {
// 		// if there is no alignment gap at this site, ++
// 		if nuc != '-' {
// 			degappedLen++
// 		}
// 	}
// 	// if there are no alignment gaps, we can return a slice of 0s
// 	if degappedLen == len(refseq) {
// 		return make([]int, len(refseq), len(refseq))
// 	}

// 	// otherwise, we get one offset:
// 	// 1) offsetRefCoord = the number of bases to add to convert each position to msa coordinates
// 	gapsum := 0
// 	offsetRefCoord := make([]int, degappedLen)
// 	for i, nuc := range refseq {
// 		if nuc == '-' {
// 			gapsum++
// 			continue
// 		}
// 		offsetRefCoord[i-gapsum] = gapsum
// 	}

// 	return offsetRefCoord
// }

// trimAlignment trims the alignment to user-specified coordinates in degapped reference space
func trimAlignment(trim bool, trimStart int, trimEnd int, cPairIn chan alignPair, cPairOut chan alignPair, cErr chan error) {

	// if no trimming is specified on the command line, we take the whole sequence:
	if !trim {
		for pair := range cPairIn {
			pair.descriptor = pair.queryname
			cPairOut <- pair
		}
		// if one feature is specified on the command line:
	} else {
		for pair := range cPairIn {

			offsetRefCoord, _ := variants.GetOffsets(pair.ref)

			adjTrimStart := trimStart + offsetRefCoord[trimStart]
			adjTrimEnd := trimEnd + offsetRefCoord[trimEnd]

			pair.query = pair.query[adjTrimStart:adjTrimEnd]
			pair.ref = pair.ref[adjTrimStart:adjTrimEnd]

			cPairOut <- pair
		}
	}
}

// writePairwiseAlignment writes the pairwise alignments between reference and queries to a directory, p, one fasta
// file per query
func writePairwiseAlignment(p string, cPair chan alignPair, cWriteDone chan bool, cErr chan error, omitRef bool) {

	_ = path.Join()

	var err error

	if p == "stdout" {
		for AP := range cPair {
			if !omitRef {
				_, err = fmt.Fprintln(os.Stdout, ">"+AP.refname)
				if err != nil {
					cErr <- err
				}
				_, err = fmt.Fprintln(os.Stdout, string(AP.ref))
				if err != nil {
					cErr <- err
				}
			}
			_, err = fmt.Fprintln(os.Stdout, ">"+AP.queryname)
			if err != nil {
				cErr <- err
			}
			_, err = fmt.Fprintln(os.Stdout, string(AP.query))
			if err != nil {
				cErr <- err
			}
		}
	} else {
		os.MkdirAll(p, 0755)

		for AP := range cPair {
			// forward slashes are illegal in unix filenames (so is ascii NUL ?)
			des := strings.ReplaceAll(AP.queryname, "/", "_")
			// unix filenames must be <= 255 chars, (account for ".fasta")
			if len(des) > 249 {
				fmt.Fprintf(os.Stderr, "Filename too long, truncating \"%s\" to: \"%s\"\n", des, des[0:249])
				des = des[0:249]
			}
			f, err := os.Create(path.Join(p, des+".fasta"))
			if err != nil {
				cErr <- err
			}
			if !omitRef {
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
	cWriteDone <- true
}

// // checkArgsPairAlign sanity checks the trimming and padding arguments (given the length of the ref seq)
// func checkArgsPairAlign(refLen int, trim bool, trimstart int, trimend int) error {

// 	if trim {
// 		if trimstart > refLen-2 || trimstart < 1 {
// 			return errors.New("error parsing trimming coordinates: check or include --trimstart")
// 		}
// 		if trimend > refLen-1 || trimend < 1 {
// 			return errors.New("error parsing trimming coordinates: check or include --trimend")
// 		}
// 		if trimstart >= trimend {
// 			return errors.New("error parsing trimming coordinates: check trimstart and trimend")
// 		}
// 	}

// 	return nil
// }

// ToPairAlign converts a SAM file into pairwise fasta-format alignments, optionally including the reference,
// optionally skipping insertions relative to the reference, optionally trimmed to coordinates in (degapped-)reference space
func ToPairAlign(samIn, ref io.Reader, outpath string, trim bool, trimStart int, trimEnd int, omitRef bool, omitIns bool, threads int) error {

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

	go fastaio.ReadAlignment(ref, cRef, cErr, cRefDone)

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

	err := checkArgs(len(refSeq), trim, trimStart, trimEnd)
	if err != nil {
		return err
	}

	cSR := make(chan samRecords, threads)
	cSH := make(chan biogosam.Header)

	cPairAlign := make(chan alignPair)
	cPairTrim := make(chan alignPair)

	cReadDone := make(chan bool)
	cAlignWaitGroupDone := make(chan bool)
	cTrimWaitGroupDone := make(chan bool)
	cWriteDone := make(chan bool)

	go groupSamRecords(samIn, cSH, cSR, cReadDone, cErr)

	_ = <-cSH

	go writePairwiseAlignment(outpath, cPairTrim, cWriteDone, cErr, omitRef)

	var wgAlign sync.WaitGroup
	wgAlign.Add(threads)

	var wgTrim sync.WaitGroup
	wgTrim.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			blockToPairwiseAlignment(cSR, cPairAlign, cErr, []byte(refSeq), omitIns)
			wgAlign.Done()
		}()
	}

	for n := 0; n < threads; n++ {
		go func() {
			trimAlignment(trim, trimStart, trimEnd, cPairAlign, cPairTrim, cErr)
			wgTrim.Done()
		}()
	}

	go func() {
		wgAlign.Wait()
		cAlignWaitGroupDone <- true
	}()

	go func() {
		wgTrim.Wait()
		cTrimWaitGroupDone <- true
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
		case <-cTrimWaitGroupDone:
			close(cPairTrim)
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
