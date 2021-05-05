package sam

import (
	"errors"
	"io"
	"os"
	"unicode"
	"unicode/utf8"

	biogosam "github.com/biogo/hts/sam"
)

// samRecords is a struct that carries a group of sam lines (belonging to the
// same sequence, probably) and a integer index which is used to keep track of
// the order of the input when we parallelise
type samRecords struct {
	records []biogosam.Record
	idx	int
}

// getOneLine processes one non-header line of a SAM file into an aligned sequence
func getOneLine(samLine biogosam.Record, refLen int, includeInsertions bool) ([]byte, error) {

	var lambda_dict map[string]func(int, int, int, []byte) (int, int, []byte)

	if includeInsertions {
		lambda_dict = getCigarOperationMapWithInsertions()
	} else {
		lambda_dict = getCigarOperationMapNoInsertions()
	}

	// QNAME := samLine.Name

	POS := samLine.Pos

	if POS < 0 {
		return []byte{}, errors.New("unmapped read")
	}

	SEQ := samLine.Seq.Expand()

	CIGAR := samLine.Cigar

	newSeqArray := make([]byte, POS)
	for i, _ := range newSeqArray {
		newSeqArray[i] = '*'
	}

	qstart := 0
	rstart := POS

	for _, op := range CIGAR {
		// fmt.Println(op.Type().String())
		// fmt.Println(op.Len())

		operation := op.Type().String()
		size := op.Len()

		new_qstart, new_rstart, extension := lambda_dict[operation](qstart, rstart, size, SEQ)

		newSeqArray = append(newSeqArray, extension...)

		qstart = new_qstart
		rstart = new_rstart

	}

	if ! includeInsertions {
		rightpad := make([]byte, refLen-len(newSeqArray))
		for i, _ := range rightpad {
			rightpad[i] = '*'
		}

		newSeqArray = append(newSeqArray, rightpad...)
	}

	// fmt.Println(string(newSeqArray))

	return newSeqArray, nil
}

// getOneLinePlusRef processes one non-header line of a SAM file into a pair of
// (aligned) sequences which are the query and the reference
func getOneLinePlusRef(samLine biogosam.Record, reference []byte, includeInsertions bool) ([]byte, []byte, error) {

	var lambda_dict map[string]func(int, int, int, []byte, []byte) (int, int, []byte, []byte)

	if includeInsertions {
		lambda_dict = getCigarOperationMapWithInsertionsWithRef()
	} else {
		lambda_dict = getCigarOperationMapNoInsertionsWithRef()
	}

	// QNAME := samLine.Name

	POS := samLine.Pos

	if POS < 0 {
		return []byte{}, []byte{}, errors.New("unmapped read")
	}

	SEQ := samLine.Seq.Expand()

	CIGAR := samLine.Cigar

	newSeqArray := make([]byte, POS)
	for i, _ := range newSeqArray {
		newSeqArray[i] = '*'
	}

	newRefSeqArray := make([]byte, POS)
	for i, _ := range newRefSeqArray {
		newRefSeqArray[i] = reference[i]
	}

	qstart := 0
	rstart := POS

	for _, op := range CIGAR {
		// fmt.Println(op.Type().String())
		// fmt.Println(op.Len())

		operation := op.Type().String()
		size := op.Len()

		new_qstart, new_rstart, extension, refextension := lambda_dict[operation](qstart, rstart, size, SEQ, reference)

		newSeqArray = append(newSeqArray, extension...)
		newRefSeqArray = append(newRefSeqArray, refextension...)

		// fmt.Println(extension)
		// fmt.Println(refextension)
		// fmt.Println()


		qstart = new_qstart
		rstart = new_rstart

	}

	if ! includeInsertions {
		rightpad := make([]byte, len(reference)-len(newSeqArray))
		// fmt.Println(len(rightpad))
		for i, _ := range rightpad {
			rightpad[i] = '*'
		}

		newSeqArray = append(newSeqArray, rightpad...)
	}

	// fmt.Println(string(newSeqArray))

	return newSeqArray, newRefSeqArray, nil
}

// getSetFromSlice returns the Set of bytes from an array of bytes. It is used
// to get the set of nucleotides at an alignment column within one query sequence
// in order to flatten the site properly.
func getSetFromSlice(s []byte) []byte {

	s_out := make([]byte, 0)

	m := make(map[byte]bool)

	for _, b := range s {
		m[b] = true
	}

	for key, value := range m {
		if value {
			s_out = append(s_out, key)
		}
	}

	return s_out
}

// getNucFromSite flattens a site to a single nucleotide when a query sequence
// has secondary mappings (multiple records/lines) in the SAM file.
// * If there is more than one alphabetic character at this site, an N is returned.
// * Alphabetic characters override '-'s and '*'s
func getNucFromSite(s []byte, qname string) byte {

	check := 0

	ss := getSetFromSlice(s)

	for _, e := range ss {
		r, _ := utf8.DecodeRune([]byte{e})
		if unicode.IsLetter(r) {
			check++
		}
	}

	if check > 1 {
		os.Stderr.WriteString("ambiguous overlapping alignment: " + qname + ": " + string(ss) + "\n")
		return 'N'
	}

	var m byte
	for i, e := range ss {
		if i == 0 || e > m {
			m = e
		}
	}

	return m
}

// checkAndGetFlattenedSeq applies getNucFromSite over all sites in a block
// of SAM records to get a single flattened sequence for one query
func checkAndGetFlattenedSeq(block [][]byte, qname string) []byte {

	seq := make([]byte, len(block[0]))
	site := make([]byte, len(block))

	for j, _ := range block[0] {
		for i, _ := range block {
			site[i] = block[i][j]
		}
		nuc := getNucFromSite(site, qname)
		seq[j] = nuc
	}

	return seq
}

// getSeqFromBlock wraps the above functions to get a sequence from one query's
// SAM records - if there is only one line (only a primary mapping) it
// returns that aligned sequence without needing to do any flattening
func getSeqFromBlock(records []biogosam.Record, refLen int, includeInsertions bool) ([]byte, error) {

	qname := records[0].Name

	block := make([][]byte, len(records))
	for i, _ := range block {
		block[i] = make([]byte, refLen)
	}

	for i, line := range records {
		temp, err := getOneLine(line, refLen, includeInsertions)
		if err != nil {
			return []byte{}, err
		}
		block[i] = temp
	}

	var seq []byte

	if len(block) > 1 {
		seq = checkAndGetFlattenedSeq(block, qname)
	} else {
		seq = block[0]
	}

	// fmt.Println(seq)
	return seq, nil
}

// swapInNs replaces unmapped positions with Ns
func swapInNs(seq []byte) []byte {
	for i, L := range seq {
		if L == '*' {
			seq[i] = 'N'
		}
	}
	return seq
}

// swapInGapNss replaces internal unmapped positions with Ns and external
// unmapped positions with '-'s
func swapInGapsNs(seq []byte) []byte {
	firstLetter := true
	var firstLetterIndx int
	var lastLetterIndx int

	for i, L := range seq {
		r, _ := utf8.DecodeRune([]byte{L})
		if unicode.IsLetter(r) {

			if firstLetter {
				firstLetterIndx = i
				firstLetter = false
			}

			lastLetterIndx = i

		}
	}

	for i, L := range seq {
		if i < firstLetterIndx {
			if L == '*' {
				seq[i] = '-'
			}
		}
		if i > firstLetterIndx && i < lastLetterIndx {
			if L == '*' {
				seq[i] = 'N'
			}
		}
		if i > lastLetterIndx {
			if L == '*' {
				seq[i] = '-'
			}
		}
	}

	return seq
}

// // getSamHeader uses Biogo/sam to return the header of a SAM file
// func getSamHeader(infile string) (biogosam.Header, error) {
//
// 	var err error
// 	f := os.Stdin
//
// 	if len(infile) > 0 {
// 		f, err = os.Open(infile)
// 		if err != nil {
// 			return biogosam.Header{}, err
// 		}
// 	}
//
// 	defer f.Close()
//
// 	s, err := biogosam.NewReader(f)
// 	if err != nil {
// 		return biogosam.Header{}, err
// 	}
//
// 	header := *s.Header()
//
// 	return header, nil
// }

// groupSamRecords yields blocks of SAM records that correspond to the same query
// sequence (to a channel)
func groupSamRecords(infile string, cHeader chan biogosam.Header, chnl chan samRecords, cdone chan bool, cerr chan error) {

	var err error
	f := os.Stdin

	if len(infile) > 0 {
		f, err = os.Open(infile)
		if err != nil {
			cerr <- err
		}
	}

	defer f.Close()

	s, err := biogosam.NewReader(f)
	if err != nil {
		cerr <- err
	}

	cHeader<- *s.Header()

	// this counter will be used to preserve order in input and output:
	counter := 0

	first := true
	samLineGroup := samRecords{idx: counter}
	var previous string

	for {

		rec, err := s.Read()

		if err == io.EOF {

			break

		} else if err != nil {

			cerr <- err

		} else {
			// if this read is unmapped, then skip it.
			// the third bit (== 4) in the sam flag is set if the read is unmapped,
			// can use the rightshift method to check this:
			if ((rec.Flags >> 2) & 1) == 1 {
				os.Stderr.WriteString("skipping unmapped read: " + rec.Name + "\n")
				continue
			}

			// if this mapping is secondary, then skip it.
			// the 9th bit (== 256) in the sam flag is set if the mapping is secondary,
			// can use the rightshift method to check this:
			if ((rec.Flags >> 8) & 1) == 1 {
				os.Stderr.WriteString("ignoring secondary mapping: " + rec.Name + "\n")
				continue
			}

			if first {
				samLineGroup.records = append(samLineGroup.records, *rec)
				first = false
				previous = rec.Name
				continue
			}

			if rec.Name != previous {
				chnl <- samLineGroup
				counter++

				samLineGroup = samRecords{idx: counter}
				samLineGroup.records = append(samLineGroup.records, *rec)
				previous = rec.Name
				continue
			}

			samLineGroup.records = append(samLineGroup.records, *rec)
			previous = rec.Name

		}

	}

	if len(samLineGroup.records) > 0 {
		chnl <- samLineGroup
	}

	cdone <- true
}
