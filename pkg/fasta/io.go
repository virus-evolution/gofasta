package fasta

import (
	"bufio"
	"bytes"
	"errors"
	"io"
)

var (
	errBadlyFormedFasta = errors.New("badly formed fasta file")
	errEmptyFasta       = errors.New("empty fasta file")
	errDiffLenSeqs      = errors.New("different length sequences in input file: is this an alignment?")
)

type Reader struct {
	*bufio.Reader
}

func NewReader(f io.Reader) *Reader {
	return &Reader{bufio.NewReader(f)}
}

// Read reads one fasta record from the underlying reader. The final record is
// returned with error = nil, and the next call to Read() returns an empty Record
// struct and error = io.EOF.
func (r *Reader) Read() (Record, error) {

	var (
		buffer, line, peek []byte
		fields             [][]byte
		err                error
		FR                 Record
	)

	first := true

	for {
		if first {
			// "ReadBytes reads until the first occurrence of delim in the input,
			// returning a slice containing the data up to and including the delimiter.
			// If ReadBytes encounters an error before finding a delimiter,
			// it returns the data read before the error and the error itself (often io.EOF).
			// ReadBytes returns err != nil if and only if the returned data does not end in delim.
			// For simple uses, a Scanner may be more convenient."
			line, err = r.ReadBytes('\n')

			// A header line that is the last thing in the file without a trailing newline comes back
			// with err == io.EOF but non-empty data. It is kept, as a record with an empty sequence, the
			// same as a final header line that does end in a newline.
			// (The next Peek() below will see the EOF and end the record.)
			if err == io.EOF && len(line) > 0 {
				err = nil
			}

			if err != nil {
				return Record{}, err

				// if the header doesn't start with a > then something is also wrong
			} else if line[0] != '>' {
				return Record{}, errBadlyFormedFasta
			}

			drop := 0
			// Strip unix or dos newline characters from the header before setting the description.
			if line[len(line)-1] == '\n' {
				drop = 1
				if len(line) > 1 && line[len(line)-2] == '\r' {
					drop = 2
				}
				line = line[:len(line)-drop]
			}

			// split the header on whitespace
			fields = bytes.Fields(line[1:])
			// fasta ID (empty if the header has nothing after the ">")
			if len(fields) > 0 {
				FR.ID = string(fields[0])
			}
			// fasta description
			FR.Description = string(line[1:])

			// we are no longer on a header line
			first = false

		} else {
			// peek at the first next byte of the underlying reader, in order
			// to see if we've reached the end of this record (or the file)
			peek, err = r.Peek(1)

			// both these cases are fine if first = false, so we can exit the loop and return the fasta record
			if err == io.EOF || (err == nil && peek[0] == '>') {
				err = nil
				break

				// other errors are returned along with an empty fasta record
			} else if err != nil {
				return Record{}, err
			}

			// If we've got this far, this should be a sequence line.
			// The err from ReadBytes() may be io.EOF if the file ends before a newline character, but this is okay because it will
			// be caught when we peek in the next iteration of the while loop.
			line, err = r.ReadBytes('\n')
			if err != nil && err != io.EOF {
				return Record{}, err
			}

			drop := 0
			// Strip unix or dos newline characters from the sequence before appending it.
			if line[len(line)-1] == '\n' {
				drop = 1
				if len(line) > 1 && line[len(line)-2] == '\r' {
					drop = 2
				}
				line = line[:len(line)-drop]
			}

			buffer = append(buffer, line...)
		}
	}
	FR.Seq = string(buffer)

	return FR, err
}

// NewAlignmentSource returns a function that reads the next record of an alignment in fasta format each
// time it is called, setting its Idx. After the last record it returns io.EOF. It returns an error if the
// input has no records, or if a record's length differs from the first's.
func NewAlignmentSource(f io.Reader) func() (Record, error) {
	r := NewReader(f)
	counter := 0
	var width int
	return func() (Record, error) {
		record, err := r.Read()
		if err == io.EOF {
			if counter == 0 {
				return Record{}, errEmptyFasta
			}
			return Record{}, io.EOF
		}
		if err != nil {
			return Record{}, err
		}
		record.Idx = counter
		if counter == 0 {
			width = len(record.Seq)
		} else if width != len(record.Seq) {
			return Record{}, errDiffLenSeqs
		}
		counter++
		return record, nil
	}
}

// NewEncodedAlignmentSource is as NewAlignmentSource, but encodes each record (see Record.Encode;
// Record.EncodeHardGaps if hardGaps), optionally calculating its base content (atgc) and its
// completeness (score).
func NewEncodedAlignmentSource(f io.Reader, hardGaps, atgc, score bool) func() (EncodedRecord, error) {
	next := NewAlignmentSource(f)
	return func() (EncodedRecord, error) {
		record, err := next()
		if err != nil {
			return EncodedRecord{}, err
		}
		var encodedRecord EncodedRecord
		if hardGaps {
			encodedRecord, err = record.EncodeHardGaps()
		} else {
			encodedRecord, err = record.Encode()
		}
		if err != nil {
			return EncodedRecord{}, err
		}
		if score {
			encodedRecord.CalculateCompleteness()
		}
		if atgc {
			encodedRecord.CalculateBaseContent()
		}
		return encodedRecord, nil
	}
}

// StreamAlignment reads an alignment in fasta format to a channel of Record structs
func StreamAlignment(f io.Reader, cR chan Record, cErr chan error, cDone chan bool) {
	next := NewAlignmentSource(f)
	for {
		record, err := next()
		if err == io.EOF {
			break
		}
		if err != nil {
			cErr <- err
			return
		}
		cR <- record
	}
	cDone <- true
}

// StreamEncodeAlignment reads an alignment in fasta format to a channel of
// Record structs - converting the nucleotide sequence to EP's bitwise coding
// scheme optionally with or without hard gaps, scoring each record and getting
// ATGC counts
func StreamEncodeAlignment(
	f io.Reader,
	cER chan EncodedRecord,
	cErr chan error,
	cDone chan bool,
	hardGaps bool,
	atgc bool,
	score bool) {

	next := NewEncodedAlignmentSource(f, hardGaps, atgc, score)
	for {
		record, err := next()
		if err == io.EOF {
			break
		}
		if err != nil {
			cErr <- err
			return
		}
		cER <- record
	}
	cDone <- true
}

// LoadEncodeAlignment is as StreamEncodeAlignment but returns a slice of Records
// instead of passing each Record down a channel
func LoadEncodeAlignment(
	f io.Reader,
	hardGaps bool,
	atgc bool,
	score bool) ([]EncodedRecord, error) {

	cER := make(chan EncodedRecord)
	cErr := make(chan error)
	cDone := make(chan bool)
	encodedRecords := make([]EncodedRecord, 0)

	go StreamEncodeAlignment(f, cER, cErr, cDone, hardGaps, atgc, score)

	for n := 1; n > 0; {
		select {
		case encodedRecord := <-cER:
			encodedRecords = append(encodedRecords, encodedRecord)
		case err := <-cErr:
			return make([]EncodedRecord, 0), err
		case <-cDone:
			n--
		}
	}

	return encodedRecords, nil
}

// WriteAlignment reads Records from a channel and writes them to file or stdout,
// in the order in which they are present in the input file.
// It passes a true to a done channel when the channel of fasta records is empty
func WriteAlignment(cR chan Record, w io.Writer, cErr chan error, cDone chan bool) {
	outputMap := make(map[int]Record)
	counter := 0
	var err error
	for FR := range cR {
		outputMap[FR.Idx] = FR
		for {
			if record, ok := outputMap[counter]; ok {
				_, err = w.Write([]byte(">" + record.ID + "\n"))
				if err != nil {
					cErr <- err
					return
				}
				_, err = w.Write([]byte(record.Seq + "\n"))
				if err != nil {
					cErr <- err
					return
				}
				delete(outputMap, counter)
				counter++
			} else {
				break
			}
		}
	}
	cDone <- true
}

// WriteWrapAlignment is as WriteAlignment but with sequence lines wrapped to w
// characters in length
func WriteWrapAlignment(cR chan Record, w io.Writer, wrap int, cErr chan error, cDone chan bool) {
	outputMap := make(map[int]Record)
	var (
		counter, written int
		err              error
	)
	for R := range cR {
		outputMap[R.Idx] = R
		for {
			if record, ok := outputMap[counter]; ok {
				_, err = w.Write([]byte(">" + record.ID + "\n"))
				if err != nil {
					cErr <- err
					return
				}
				for {
					if written < len(record.Seq) {
						if written+wrap >= len(record.Seq) {
							_, err = w.Write([]byte(record.Seq[written:] + "\n"))
							if err != nil {
								cErr <- err
								return
							}
							written = written + wrap
						} else {
							_, err = w.Write([]byte(record.Seq[written:written+wrap] + "\n"))
							if err != nil {
								cErr <- err
								return
							}
							written = written + wrap
						}
					} else {
						break
					}
				}
				delete(outputMap, counter)
				counter++
				written = 0
			} else {
				break
			}
		}

	}
	cDone <- true
}

// getAlignmentDims gets the dimensions (i.e. the number of sequences and the
// width of the alignment in number of nucleotides) of an alignment in fasta
// format
func getAlignmentDims(f io.Reader) (int, int, error) {
	n := 0
	l := 0
	var err error
	r := NewReader(f)
	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return 0, 0, err
		}
		if n == 0 {
			l = len(record.Seq)
		}
		n++
	}
	return n, l, err
}
