/*
Package fastaio provides functions for reading and writing
fasta format files
*/
package fastaio

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/encoding"
)

// A struct for one Fasta record
type FastaRecord struct {
	ID          string
	Description string
	Seq         string
	Idx         int
}

// A struct for one Fasta record whose sequence is encoded using EP's scheme
type EncodedFastaRecord struct {
	ID          string
	Description string
	Seq         []byte
	Score       int64 // this is for e.g., genome completeness
	Idx         int
	Count_A     int
	Count_T     int
	Count_G     int
	Count_C     int
}

// TO DO - test this
func (FR FastaRecord) Encode() EncodedFastaRecord {
	EFR := EncodedFastaRecord{ID: FR.ID, Description: FR.Description, Idx: FR.Idx}
	EA := encoding.MakeEncodingArray()
	seq := make([]byte, len(FR.Seq))
	for i, nuc := range FR.Seq {
		seq[i] = EA[nuc]
	}
	EFR.Seq = seq
	return EFR
}

// TO DO - test this
func (EFR EncodedFastaRecord) Decode() FastaRecord {
	FR := FastaRecord{ID: EFR.ID, Description: EFR.Description, Idx: EFR.Idx}
	DA := encoding.MakeDecodingArray()
	seq := ""
	for _, nuc := range EFR.Seq {
		seq = seq + DA[nuc]
	}
	FR.Seq = seq
	return FR
}

func (EEFT *EncodedFastaRecord) CalculateBaseContent() {
	count_A := 0
	count_T := 0
	count_G := 0
	count_C := 0

	for _, nuc := range EEFT.Seq {
		switch nuc {
		case 136:
			count_A++
		case 72:
			count_G++
		case 40:
			count_C++
		case 24:
			count_T++
		}
	}

	EEFT.Count_A = count_A
	EEFT.Count_C = count_C
	EEFT.Count_G = count_G
	EEFT.Count_T = count_T
}

// getAlignmentDims gets the dimensions (i.e. the number of sequences and the width of the alignment
// in numner of nucleotides) of an alignment in fasta format
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

// To do - treat ambiguities as the set of nucleotides they represent?
// could do this most easily with EP's coding scheme?
func Consensus(f io.Reader) (FastaRecord, error) {

	score := [256]byte{}
	for i := range score {
		score[i] = 17
	}
	score['A'] = 0
	score['C'] = 1
	score['G'] = 2
	score['T'] = 3
	score['R'] = 4
	score['M'] = 5
	score['W'] = 6
	score['S'] = 7
	score['K'] = 8
	score['Y'] = 9
	score['V'] = 10
	score['H'] = 11
	score['D'] = 12
	score['B'] = 13
	score['N'] = 14
	score['-'] = 15
	score['?'] = 16

	var err error

	cFR := make(chan FastaRecord)
	cErr := make(chan error)
	cReadDone := make(chan bool)
	cConsensusDone := make(chan bool)

	go ReadAlignment(f, cFR, cErr, cReadDone)

	var counts [][17]byte

	go func() {
		first := true

		for FR := range cFR {
			if first {
				counts = make([][17]byte, len(FR.Seq))
				first = false
			}
			for i := 0; i < len(FR.Seq); i++ {
				counts[i][score[FR.Seq[i]]]++
			}
		}

		cConsensusDone <- true
	}()

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return FastaRecord{}, err
		case <-cReadDone:
			close(cFR)
			n--
		}
	}

	<-cConsensusDone

	backTranslate := [17]string{"A", "C", "G", "T", "R", "M", "W", "S", "K", "Y", "V", "H", "D", "B", "N", "-", "?"}
	consensusSeq := ""
	for _, a := range counts {
		maxval := a[0]
		maxidx := 0
		for j := 1; j < len(a); j++ {
			if a[j] > maxval {
				maxval = a[j]
				maxidx = j
			}
		}
		consensusSeq = consensusSeq + backTranslate[maxidx]
	}

	consensus := FastaRecord{ID: "consensus", Seq: consensusSeq}

	return consensus, err
}

// ReadAlignment reads an alignment in fasta format to a channel of FastaRecord structs
func ReadAlignment(f io.Reader, chnl chan FastaRecord, cErr chan error, cdone chan bool) {

	var err error
	s := bufio.NewScanner(f)

	counter := 0

	first := true

	var id string
	var description string
	var seqBuffer string
	var width int

	for s.Scan() {
		line := s.Text()

		if first {

			if len(line) == 0 || string(line[0]) != ">" {
				cErr <- errors.New("badly formatted fasta file")
				return
			}

			description = line[1:]
			id = strings.Fields(description)[0]

			first = false

		} else if string(line[0]) == ">" {

			if counter == 0 {
				width = len(seqBuffer)
			} else if len(seqBuffer) != width {
				cErr <- errors.New("different length sequences in input file: is this an alignment?")
				return
			}

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

	if len(seqBuffer) > 0 {
		if counter > 0 && len(seqBuffer) != width {
			cErr <- errors.New("different length sequences in input file: is this an alignment?")
			return
		}
		fr := FastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
		chnl <- fr
		counter++
	}

	if counter == 0 {
		cErr <- errors.New("empty fasta file")
		return
	}

	err = s.Err()
	if err != nil {
		cErr <- err
		return
	}

	cdone <- true
}

// ReadEncodeAlignment reads an alignment in fasta format to a channel
// of EncodedFastaRecord structs - converting sequence to EP's bitwise coding scheme
func ReadEncodeAlignment(f io.Reader, hardGaps bool, chnl chan EncodedFastaRecord, cErr chan error, cDone chan bool) {

	var err error

	var coding [256]byte
	switch hardGaps {
	case true:
		coding = encoding.MakeEncodingArrayHardGaps()
	case false:
		coding = encoding.MakeEncodingArray()
	}

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

			if len(line) == 0 || line[0] != '>' {
				cErr <- errors.New("badly formatted fasta file")
				return
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

	if len(seqBuffer) > 0 {
		if counter > 0 && len(seqBuffer) != width {
			cErr <- errors.New("different length sequences in input file: is this an alignment?")
			return
		}
		fr = EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
		chnl <- fr
		counter++
	}

	if counter == 0 {
		cErr <- errors.New("empty fasta file")
		return
	}

	err = s.Err()
	if err != nil {
		cErr <- err
		return
	}

	cDone <- true
}

// ReadEncodeAlignmentToList is as ReadEncodeAlignment but returns an array of EncodedFastaRecords instead
// of passing each EncodedFastaRecord to a channel
func ReadEncodeAlignmentToList(f io.Reader, hardGaps bool) ([]EncodedFastaRecord, error) {

	var err error

	records := make([]EncodedFastaRecord, 0)

	var coding [256]byte
	switch hardGaps {
	case true:
		coding = encoding.MakeEncodingArrayHardGaps()
	case false:
		coding = encoding.MakeEncodingArray()
	}

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

			if len(line) == 0 || line[0] != '>' {
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

	if len(seqBuffer) > 0 {
		if counter > 0 && len(seqBuffer) != width {
			return []EncodedFastaRecord{}, errors.New("different length sequences in input file: is this an alignment?")
		}
		fr := EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
		records = append(records, fr)
		counter++
	}

	if counter == 0 {
		return []EncodedFastaRecord{}, errors.New("empty fasta file")
	}

	err = s.Err()
	if err != nil {
		return []EncodedFastaRecord{}, err
	}

	return records, nil
}

// WriteAlignment reads fasta records from a channel and writes them to file or stdout,
// in the order in which they are present in the input file.
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
