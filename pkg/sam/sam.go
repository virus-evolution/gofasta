/*
Package sam provides functionality to handle sam format files.

It was written to process pairwise alignments between assembled SARS-CoV-2
genomes from minimap2.

Its main purpose is to generate fasta format alignments. It also has a routine
to list mutations between pairs of sequences.
*/
package sam

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"sync"
	"unicode"
	"unicode/utf8"

	biogosam "github.com/biogo/hts/sam"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
)

// samRecords is a struct that carries a group of sam lines (belonging to the
// same sequence, probably) and a integer index which is used to keep track of
// the order of the input when we parallelise
type samRecords struct {
	name    string // the query's name (QNAME)
	records []biogosam.Record
	idx     int
	err     error // why one of the query's lines could not be parsed, if one could not
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

// parseSamHeader parses the header lines of a sam file. The biogo sam parser panics on some malformed
// input (as well as returning errors), so panics are turned into errors.
func parseSamHeader(text []byte) (h *biogosam.Header, err error) {
	defer func() {
		if p := recover(); p != nil {
			err = fmt.Errorf("error parsing the sam header: %v", p)
		}
	}()
	h, err = biogosam.NewHeader(nil, nil)
	if err != nil {
		return nil, err
	}
	if err := h.UnmarshalText(text); err != nil {
		return nil, err
	}
	return h, nil
}

// parseSamRecord parses one alignment line of a sam file, turning panics in the biogo parser into errors
func parseSamRecord(h *biogosam.Header, line []byte) (rec biogosam.Record, err error) {
	defer func() {
		if p := recover(); p != nil {
			err = fmt.Errorf("error parsing a sam record: %v", p)
		}
	}()
	err = rec.UnmarshalSAM(h, line)
	return rec, err
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

// samSource reads a sam file one query at a time (see newSamSource). The messages about skipped lines are
// kept with the query they were read with, so that they can be written in input order with the
// pipeline's own warnings (see withNotes).
type samSource struct {
	header biogosam.Header
	next   func() (samRecords, error)

	mu    sync.Mutex
	notes map[int][]string // by the index of the query being read when the line was skipped
	tail  []string         // lines skipped after the last query
}

func (s *samSource) note(idx int, msg string) {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.notes[idx] = append(s.notes[idx], msg)
}

// withNotes wraps a sink so that the messages about the lines skipped while reading each query are written
// to stderr before the query's result is handled.
func withNotes[T any](s *samSource, sink func(pipeline.Result[T]) error) func(pipeline.Result[T]) error {
	return func(r pipeline.Result[T]) error {
		s.mu.Lock()
		msgs := s.notes[r.Idx]
		delete(s.notes, r.Idx)
		s.mu.Unlock()
		for _, m := range msgs {
			os.Stderr.WriteString(m)
		}
		return sink(r)
	}
}

// writeTail writes the messages about the lines skipped after the last query
func (s *samSource) writeTail() {
	s.mu.Lock()
	defer s.mu.Unlock()
	for _, m := range s.tail {
		os.Stderr.WriteString(m)
	}
	s.tail = nil
}

// newSamSource reads the header of a sam file and returns a samSource whose next function returns the
// records for the next query each time it is called, and io.EOF after the last query. The records for a query are
// consecutive lines with the same QNAME; unmapped and secondary records, and records mapped to any
// reference other than the first @SQ, are skipped (with a message on stderr), and blank lines are
// ignored. A line that can't be parsed is recorded as the error of its query (whose other lines are still
// read, so that the next query starts in the right place). The records are not checked here: see
// prepareQuery. An error from the returned function means the file could not be read.
func newSamSource(sam io.Reader) (*samSource, error) {

	br := bufio.NewReader(sam)

	// readLine returns the next line without its line ending (a final line need not have one), or io.EOF
	readLine := func() ([]byte, error) {
		line, err := br.ReadBytes('\n')
		if err == io.EOF && len(line) > 0 {
			err = nil
		}
		if err != nil {
			return nil, err
		}
		line = bytes.TrimSuffix(line, []byte{'\n'})
		line = bytes.TrimSuffix(line, []byte{'\r'})
		return line, nil
	}

	if _, err := br.Peek(1); err != nil {
		return nil, err
	}
	var headerText []byte
	for {
		p, err := br.Peek(1)
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if p[0] != '@' {
			break
		}
		line, err := readLine()
		if err != nil {
			return nil, err
		}
		headerText = append(append(headerText, line...), '\n')
	}
	h, err := parseSamHeader(headerText)
	if err != nil {
		return nil, err
	}
	if len(h.Refs()) == 0 {
		return nil, errors.New("no reference sequence (@SQ line) in the sam header")
	}

	src := &samSource{header: *h, notes: make(map[int][]string)}
	counter := 0     // the index of the query being read
	started := false // whether a line of that query has been read
	// a skipped line is listed with the query that comes after it in the file: the one being read if
	// none of its lines has been read yet, otherwise (as far as we know yet) the next one
	noteIdx := func() int {
		if started {
			return counter + 1
		}
		return counter
	}

	// a parsed line: its QNAME, and the record or why it could not be parsed
	type samLine struct {
		name string
		rec  biogosam.Record
		err  error
	}

	// readUsable returns the next line that is not skipped, or io.EOF
	readUsable := func() (samLine, error) {
		for {
			line, err := readLine()
			if err != nil {
				return samLine{}, err
			}
			if len(line) == 0 {
				continue
			}
			name, _, _ := bytes.Cut(line, []byte{'\t'})
			rec, err := parseSamRecord(h, line)
			if err != nil {
				return samLine{name: string(name), err: err}, nil
			}
			// the third bit (== 4) in the sam flag is set if the read is unmapped
			if ((rec.Flags >> 2) & 1) == 1 {
				src.note(noteIdx(), "skipping unmapped read: "+rec.Name+"\n")
				continue
			}
			// the 9th bit (== 256) in the sam flag is set if the mapping is secondary
			if ((rec.Flags >> 8) & 1) == 1 {
				src.note(noteIdx(), "ignoring secondary mapping: "+rec.Name+"\n")
				continue
			}
			// only one reference sequence is supported: records mapped to any other @SQ are skipped
			if rec.Ref != nil && rec.Ref.ID() != 0 {
				src.note(noteIdx(), "skipping record mapped to a reference other than the first @SQ ("+rec.Ref.Name()+"): "+rec.Name+"\n")
				continue
			}
			return samLine{name: string(name), rec: rec}, nil
		}
	}

	var pending *samLine // the first line of the next query, once it has been read
	finished := false

	next := func() (samRecords, error) {
		group := samRecords{idx: counter}
		started = false
		for !finished {
			var l samLine
			if pending != nil {
				l = *pending
				pending = nil
			} else {
				l, err = readUsable()
				if err == io.EOF {
					finished = true
					break
				}
				if err != nil {
					return samRecords{}, err
				}
			}
			if started && l.name != group.name {
				pending = &l
				break
			}
			started = true
			group.name = l.name
			if l.err != nil {
				if group.err == nil {
					group.err = l.err
				}
				continue
			}
			group.records = append(group.records, l.rec)
		}
		if !started {
			src.mu.Lock()
			src.tail = append(src.tail, src.notes[counter]...)
			delete(src.notes, counter)
			src.mu.Unlock()
			return samRecords{}, io.EOF
		}
		counter++
		return group, nil
	}

	src.next = next
	return src, nil
}

// prepareQuery checks the records for one query against a reference of length refLen and resolves the
// overlaps between them. The error, if any, carries the query's name.
func prepareQuery(group samRecords, refLen int) (samRecords, error) {
	if group.err != nil {
		return samRecords{}, pipeline.Named(group.name, group.err)
	}
	for i := range group.records {
		if err := validateRecord(&group.records[i], refLen); err != nil {
			return samRecords{}, pipeline.Named(group.name, err)
		}
	}
	group.records = resolveOverlaps(group.records)
	return group, nil
}

// queryNoun names the records of the sam pipelines (one query's records) in warnings
var queryNoun = pipeline.Noun{One: "query", Many: "queries"}
