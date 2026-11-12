/*
Package sam provides functionality to handle sam format files.

It was written to process pairwise alignments between assembled SARS-CoV-2
genomes from minimap2.

Its main purpose is to generate fasta format alignments. It also has a routine
to list mutations between pairs of sequences.
*/
package sam

import (
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"unicode"
	"unicode/utf8"

	biogosam "github.com/biogo/hts/sam"
)

// samRecords is a struct that carries a group of sam lines (belonging to the
// same sequence, probably) and a integer index which is used to keep track of
// the order of the input when we parallelise
type samRecords struct {
	records []biogosam.Record
	idx     int
}

// getOneLine processes one non-header line of a SAM file into an aligned sequence
func getOneLine(samLine biogosam.Record, refLen int, includeInsertions bool) ([]byte, error) {

	var lambda_dict map[string]func(int, int, int, []byte) (int, int, []byte)

	if includeInsertions {
		lambda_dict = getCigarOperationMapWithInsertions()
	} else {
		lambda_dict = getCigarOperationMapNoInsertions()
	}

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

		operation := op.Type().String()
		size := op.Len()

		new_qstart, new_rstart, extension := lambda_dict[operation](qstart, rstart, size, SEQ)

		newSeqArray = append(newSeqArray, extension...)

		qstart = new_qstart
		rstart = new_rstart

	}

	if !includeInsertions {
		rightpad := make([]byte, refLen-len(newSeqArray))
		for i, _ := range rightpad {
			rightpad[i] = '*'
		}

		newSeqArray = append(newSeqArray, rightpad...)
	}

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

		operation := op.Type().String()
		size := op.Len()

		new_qstart, new_rstart, extension, refextension := lambda_dict[operation](qstart, rstart, size, SEQ, reference)

		newSeqArray = append(newSeqArray, extension...)
		newRefSeqArray = append(newRefSeqArray, refextension...)

		qstart = new_qstart
		rstart = new_rstart

	}

	if !includeInsertions {
		rightpad := make([]byte, len(reference)-len(newSeqArray))

		for i, _ := range rightpad {
			rightpad[i] = '*'
		}

		newSeqArray = append(newSeqArray, rightpad...)
	}

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
func getNucFromSite(s []byte, qname string, pos int) byte {

	check := 0

	ss := getSetFromSlice(s)

	for _, e := range ss {
		r, _ := utf8.DecodeRune([]byte{e})
		if unicode.IsLetter(r) {
			check++
		}
	}

	if check > 1 {
		os.Stderr.WriteString("ambiguous overlapping alignment: " + qname + ": " + strconv.Itoa(pos+1) + ": " + string(ss) + "\n")
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
		nuc := getNucFromSite(site, qname, j)
		seq[j] = nuc
	}

	return seq
}

// getSeqFromBlock wraps the other functions in this file to get a sequence from one query's
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

// The biogo sam parser panics on some malformed input (e.g. a CIGAR string that ends in a number); these
// wrappers turn its panics into errors.

func newSamReader(r io.Reader) (s *biogosam.Reader, err error) {
	defer func() {
		if p := recover(); p != nil {
			err = fmt.Errorf("error parsing the sam header: %v", p)
		}
	}()
	return biogosam.NewReader(r)
}

func readSamRecord(s *biogosam.Reader) (rec *biogosam.Record, err error) {
	defer func() {
		if p := recover(); p != nil {
			err = fmt.Errorf("error parsing a sam record: %v", p)
		}
	}()
	return s.Read()
}

// validateRecord checks that a mapped sam record is internally consistent and fits the reference, so that
// the code that builds alignments from CIGAR strings can index its sequence without going out of range.
func validateRecord(rec *biogosam.Record, refLen int) error {

	bad := func(msg string) error {
		return errors.New("bad sam record " + rec.Name + ": " + msg)
	}

	if rec.Ref == nil {
		return bad("no reference sequence (RNAME is *) but not flagged as unmapped")
	}
	if rec.Pos < 0 {
		return bad("negative position")
	}

	for _, op := range rec.Cigar {
		switch op.Type() {
		case biogosam.CigarMatch, biogosam.CigarEqual, biogosam.CigarMismatch, biogosam.CigarInsertion,
			biogosam.CigarSoftClipped, biogosam.CigarDeletion, biogosam.CigarSkipped, biogosam.CigarHardClipped,
			biogosam.CigarPadded:
		default:
			return bad("unsupported CIGAR operation " + op.Type().String())
		}
	}
	refSpan, queryLen := rec.Cigar.Lengths()
	refEnd := rec.Pos + refSpan

	if queryLen != rec.Seq.Length {
		return bad("the CIGAR string describes " + strconv.Itoa(queryLen) + " query bases but SEQ has " + strconv.Itoa(rec.Seq.Length))
	}
	if refEnd > refLen {
		return bad("the alignment extends beyond the end of the reference (" + strconv.Itoa(refEnd) + " > " + strconv.Itoa(refLen) + ")")
	}

	return nil
}

// newSamSource reads the header of a sam file and returns it, with a function that returns the records
// for the next query each time it is called, and io.EOF after the last query. The records for a query are
// consecutive records with the same QNAME; unmapped and secondary records, and records mapped to any
// reference other than the first @SQ, are skipped (with a message on stderr), every other record is
// validated, and overlaps between a query's records are resolved.
func newSamSource(sam io.Reader) (biogosam.Header, func() (samRecords, error), error) {

	s, err := newSamReader(sam)
	if err != nil {
		return biogosam.Header{}, nil, err
	}

	if len(s.Header().Refs()) == 0 {
		return biogosam.Header{}, nil, errors.New("no reference sequence (@SQ line) in the sam header")
	}
	refLen := s.Header().Refs()[0].Len()

	// readUsable returns the next record that is not skipped, or io.EOF
	readUsable := func() (*biogosam.Record, error) {
		for {
			rec, err := readSamRecord(s)
			if err != nil {
				return nil, err
			}
			// the third bit (== 4) in the sam flag is set if the read is unmapped
			if ((rec.Flags >> 2) & 1) == 1 {
				os.Stderr.WriteString("skipping unmapped read: " + rec.Name + "\n")
				continue
			}
			// the 9th bit (== 256) in the sam flag is set if the mapping is secondary
			if ((rec.Flags >> 8) & 1) == 1 {
				os.Stderr.WriteString("ignoring secondary mapping: " + rec.Name + "\n")
				continue
			}
			// only one reference sequence is supported: records mapped to any other @SQ are skipped
			if rec.Ref != nil && rec.Ref.ID() != 0 {
				os.Stderr.WriteString("skipping record mapped to a reference other than the first @SQ (" + rec.Ref.Name() + "): " + rec.Name + "\n")
				continue
			}
			if err := validateRecord(rec, refLen); err != nil {
				return nil, err
			}
			return rec, nil
		}
	}

	var pending *biogosam.Record // the first record of the next query, once it has been read
	finished := false
	counter := 0

	next := func() (samRecords, error) {
		group := samRecords{idx: counter}
		for !finished {
			rec := pending
			pending = nil
			if rec == nil {
				rec, err = readUsable()
				if err == io.EOF {
					finished = true
					break
				}
				if err != nil {
					return samRecords{}, err
				}
			}
			if len(group.records) > 0 && rec.Name != group.records[0].Name {
				pending = rec
				break
			}
			group.records = append(group.records, *rec)
		}
		if len(group.records) == 0 {
			return samRecords{}, io.EOF
		}
		group.records = resolveOverlaps(group.records)
		counter++
		return group, nil
	}

	return *s.Header(), next, nil
}
