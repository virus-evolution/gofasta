package fastaio

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"strings"

	"github.com/cov-ert/gofasta/pkg/encoding"
)

// FastaRecord is a simple struct for Fasta records
type FastaRecord struct {
	ID          string
	Description string
	Seq         string
	Idx         int
}

// EncodedFastaRecord is a struct for one Fasta record
type EncodedFastaRecord struct {
	ID          string
	Description string
	Seq         []byte
	Score       int64 // this is for e.g., genome completeness
	Idx         int
}

func getAlignmentDims(f io.Reader) (int, int, error) {
	n := 0
	l := 0

	s := bufio.NewScanner(f)

	for s.Scan() {
		line := s.Text()

		if string(line[0]) == ">" {
			n++
		}

		if n == 1 && string(line[0]) != ">" {
			l += len(line)
		}
	}

	err := s.Err()

	if err != nil {
		return 0, 0, err
	}

	return n, l, err
}

// ReadAlignment reads an alignment in fasta format to a channel
// of FastaRecord structs
func ReadAlignment(f io.Reader, chnl chan FastaRecord, chnlerr chan error, cdone chan bool) {

	var err error
	s := bufio.NewScanner(f)

	counter := 0

	first := true

	var id string
	var description string
	var seqBuffer string

	for s.Scan() {
		line := string(s.Text())

		if first {

			if string(line[0]) != ">" {
				chnlerr <- errors.New("badly formatted fasta file")
				return
			}

			description = line[1:]
			id = strings.Fields(description)[0]

			first = false

		} else if string(line[0]) == ">" {

			fr := FastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
			chnl <- fr
			counter++

			description = line[1:]
			id = strings.Fields(description)[0]
			seqBuffer = ""

		} else {
			seqBuffer = seqBuffer + strings.ToUpper(line)
		}

	}

	fr := FastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
	chnl <- fr

	err = s.Err()
	if err != nil {
		chnlerr <- err
		return
	}

	cdone <- true
}

// writeAlignmentOut reads fasta records from a channel and writes them to a single
// outfile, in the order in which they are present in the input file.
// It passes a true to a done channel when the channel of fasta records is empty
func WriteAlignment(ch chan FastaRecord, w io.Writer, cdone chan bool, cerr chan error) {

	outputMap := make(map[int]FastaRecord)

	counter := 0

	var err error

	for FR := range ch {

		outputMap[FR.Idx] = FR

		if fastarecord, ok := outputMap[counter]; ok {
			_, err = w.Write([]byte(">" + fastarecord.ID + "\n"))
			if err != nil {
				cerr <- err
			}
			_, err = w.Write([]byte(fastarecord.Seq + "\n"))
			if err != nil {
				cerr <- err
			}
			delete(outputMap, counter)
			counter++
		} else {
			continue
		}
	}

	for n := 1; n > 0; {
		if len(outputMap) == 0 {
			break
		}
		fastarecord := outputMap[counter]
		_, err = w.Write([]byte(">" + fastarecord.ID + "\n"))
		if err != nil {
			cerr <- err
		}
		_, err = w.Write([]byte(fastarecord.Seq + "\n"))
		if err != nil {
			cerr <- err
		}
		delete(outputMap, counter)
		counter++
	}

	cdone <- true
}

// // PopulateByteArrayGetNames reads a fasta format alignment into an array of
// // uint8 representations of nucleotides, and gets the names of the fasta records
// // at the same time.
// func PopulateByteArrayGetNames(infile string) ([][]uint8, []string, error) {
// 	height, width, err := getAlignmentDims(infile)

// 	if err != nil {
// 		return [][]uint8{}, []string{}, err
// 	}

// 	IDs := make([]string, 0)

// 	A := make([][]uint8, height)
// 	for i := 0; i < height; i++ {
// 		A[i] = make([]uint8, width)
// 	}

// 	byteDict := encoding.MakeByteDict()

// 	c := make(chan FastaRecord)
// 	cerr := make(chan error)
// 	cdone := make(chan bool)

// 	go ReadAlignment(infile, c, cerr, cdone)

// 	i := 0

// 	for n := 1; n > 0; {
// 		select {

// 		case err := <-cerr:
// 			return [][]uint8{}, []string{}, err

// 		case record := <-c:

// 			IDs = append(IDs, record.ID)

// 			if len(record.Seq) != width {
// 				return [][]uint8{}, []string{}, errors.New("different length sequences in input file: is this an alignment?")
// 			}

// 			for j, nuc := range record.Seq {

// 				if val, ok := byteDict[nuc]; ok {
// 					A[i][j] = val

// 				} else {
// 					return [][]uint8{}, []string{}, fmt.Errorf("invalid nucleotide in fasta file (%s)", string(nuc))
// 				}
// 			}

// 			i++

// 		case <-cdone:
// 			n--
// 		}
// 	}

// 	return A, IDs, nil
// }

// ReadEncodeAlignment reads an alignment in fasta format to a channel
// of encodedFastaRecord structs - converting sequence to EP's bitwise coding scheme
func ReadEncodeAlignment(f io.Reader, chnl chan EncodedFastaRecord, cErr chan error, cDone chan bool) {

	var err error

	coding := encoding.MakeEncodingArray()

	s := bufio.NewScanner(f)

	first := true

	var id string
	var description string
	var seqBuffer []byte
	var line []byte
	var nuc byte
	var width int

	var fr EncodedFastaRecord

	counter := 0

	for s.Scan() {
		line = s.Bytes()

		if first {

			if line[0] != '>' {
				cErr <- errors.New("badly formatted fasta file")
			}

			description = string(line[1:])
			id = strings.Fields(description)[0]

			first = false

		} else if line[0] == '>' {

			if counter == 0 {
				width = len(seqBuffer)
			} else if len(seqBuffer) != width {
				cErr <- errors.New("different length sequences in input file: is this an alignment?")
				return
			}

			fr = EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
			chnl <- fr
			counter++

			description = string(line[1:])
			id = strings.Fields(description)[0]
			seqBuffer = make([]byte, 0)

		} else {
			encodedLine := make([]byte, len(line))
			for i := range line {
				nuc = coding[line[i]]
				if nuc == 0 {
					cErr <- fmt.Errorf("invalid nucleotide in fasta file (%s)", string(line[i]))
					return
				}
				encodedLine[i] = nuc
			}
			seqBuffer = append(seqBuffer, encodedLine...)
		}
	}

	fr = EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
	chnl <- fr

	err = s.Err()
	if err != nil {
		cErr <- err
		return
	}

	cDone <- true
}

func ReadEncodeAlignmentToList(f io.Reader) ([]EncodedFastaRecord, error) {

	var err error

	records := make([]EncodedFastaRecord, 0)

	coding := encoding.MakeEncodingArray()

	s := bufio.NewScanner(f)

	first := true

	var id string
	var description string
	var seqBuffer []byte
	var line []byte
	var nuc byte
	var width int

	counter := 0

	for s.Scan() {
		line = s.Bytes()

		if first {

			if line[0] != '>' {
				return []EncodedFastaRecord{}, errors.New("badly formatted fasta file")
			}

			description = string(line[1:])
			id = strings.Fields(description)[0]

			first = false

		} else if line[0] == '>' {

			if counter == 0 {
				width = len(seqBuffer)
			} else if len(seqBuffer) != width {
				return []EncodedFastaRecord{}, errors.New("different length sequences in input file: is this an alignment?")
			}

			fr := EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
			records = append(records, fr)
			counter++

			description = string(line[1:])
			id = strings.Fields(description)[0]
			seqBuffer = make([]byte, 0)

		} else {
			encodedLine := make([]byte, len(line))
			for i := range line {
				nuc = coding[line[i]]
				if nuc == 0 {
					return []EncodedFastaRecord{}, fmt.Errorf("invalid nucleotide in fasta file (%s)", string(line[i]))
				}
				encodedLine[i] = nuc
			}
			seqBuffer = append(seqBuffer, encodedLine...)
		}
	}

	fr := EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
	records = append(records, fr)

	err = s.Err()
	if err != nil {
		return []EncodedFastaRecord{}, err
	}

	return records, nil
}

func ReadEncodeScoreAlignment(f io.Reader, chnl chan EncodedFastaRecord, cErr chan error, cDone chan bool) {

	var err error

	coding := encoding.MakeEncodingArray()
	scoring := encoding.MakeScoreArray()

	s := bufio.NewScanner(f)

	first := true

	var id string
	var description string
	var seqBuffer []byte
	var line []byte
	var nuc byte
	var score int64
	var width int

	counter := 0

	for s.Scan() {
		line = s.Bytes()

		if first {

			if line[0] != '>' {
				cErr <- errors.New("badly formatted fasta file")
			}

			description = string(line[1:])
			id = strings.Fields(description)[0]

			first = false

		} else if line[0] == '>' {

			if counter == 0 {
				width = len(seqBuffer)
			} else if len(seqBuffer) != width {
				cErr <- errors.New("different length sequences in input file: is this an alignment?")
				return
			}

			fr := EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Score: score, Idx: counter}
			chnl <- fr
			counter++

			description = string(line[1:])
			id = strings.Fields(description)[0]
			seqBuffer = make([]byte, 0)
			score = 0

		} else {
			encodedLine := make([]byte, len(line))
			for i := range line {
				nuc = coding[line[i]]
				if nuc == 0 {
					cErr <- fmt.Errorf("invalid nucleotide in fasta file (%s)", string(line[i]))
				}
				encodedLine[i] = nuc

				score += scoring[line[i]]
			}
			seqBuffer = append(seqBuffer, encodedLine...)
		}
	}

	fr := EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Score: score, Idx: counter}
	chnl <- fr

	err = s.Err()
	if err != nil {
		cErr <- err
		return
	}

	cDone <- true
}
