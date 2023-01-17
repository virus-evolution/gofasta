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

	"github.com/virus-evolution/gofasta/pkg/alphabet"
	"github.com/virus-evolution/gofasta/pkg/encoding"
)

// A struct for one Fasta record
type FastaRecord struct {
	ID          string
	Description string
	Seq         string
	Idx         int
}

// Convert a FastaRecord to an EncodedFastaRecord
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

// TO DO - test
func (FR FastaRecord) Degap() FastaRecord {
	NFR := FastaRecord{ID: FR.ID, Description: FR.Description, Idx: FR.Idx}
	t := ""
	for _, char := range FR.Seq {
		if char != '-' {
			t = t + string(char)
		}
	}
	NFR.Seq = t
	return NFR
}

// TO DO - test
func (FR FastaRecord) Complement() FastaRecord {
	NFR := FastaRecord{ID: FR.ID, Description: FR.Description, Idx: FR.Idx}
	CA := alphabet.MakeCompArray()
	ba := make([]byte, len(FR.Seq))
	for i := 0; i < len(FR.Seq); i++ {
		ba[i] = CA[FR.Seq[i]]
	}
	NFR.Seq = string(ba)
	return NFR
}

// TO DO - test
func (FR FastaRecord) ReverseComplement() FastaRecord {
	NFR := FR.Complement()
	temp := []byte(NFR.Seq)
	for i, j := 0, len(temp)-1; i < j; i, j = i+1, j-1 {
		temp[i], temp[j] = temp[j], temp[i]
	}
	NFR.Seq = string(temp)
	return NFR
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

// Convert an EncodedFastaRecord to a FastaRecord
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

// Calculate the ATGC content of an EncodedFastaRecord
func (EFR *EncodedFastaRecord) CalculateBaseContent() {

	var counting [256]int

	for _, nuc := range EFR.Seq {
		counting[nuc]++
	}

	EFR.Count_A = counting[136]
	EFR.Count_T = counting[24]
	EFR.Count_G = counting[72]
	EFR.Count_C = counting[40]
}

// TO DO - test
func (EFR EncodedFastaRecord) Complement() EncodedFastaRecord {
	NFR := EncodedFastaRecord{ID: EFR.ID, Description: EFR.Description, Idx: EFR.Idx}
	CA := alphabet.MakeEncodedCompArray()
	NFR.Seq = make([]byte, len(EFR.Seq))
	for i := 0; i < len(EFR.Seq); i++ {
		NFR.Seq[i] = CA[EFR.Seq[i]]
	}
	return NFR
}

// TO DO - test
func (EFR EncodedFastaRecord) ReverseComplement() EncodedFastaRecord {
	NEFR := EFR.Complement()
	for i, j := 0, len(NEFR.Seq)-1; i < j; i, j = i+1, j-1 {
		NEFR.Seq[i], NEFR.Seq[j] = NEFR.Seq[j], NEFR.Seq[i]
	}
	return NEFR
}

// getAlignmentDims gets the dimensions (i.e. the number of sequences and the width of the alignment
// in numner of nucleotides) of an alignment in fasta format
func getAlignmentDims(f io.Reader) (int, int, error) {
	n := 0
	l := 0

	s := bufio.NewScanner(f)
	s.Buffer(make([]byte, 0), 1024*1024)

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

// To DO - treat ambiguities as the set of nucleotides they represent?
// could do this most easily with EP's coding scheme?

// Generate a consensus nucleotide sequence
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
	s.Buffer(make([]byte, 0), 1024*1024)

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
// of EncodedFastaRecord structs - converting the nucleotide sequence to EP's bitwise coding scheme
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
	s.Buffer(make([]byte, 0), 1024*1024)

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
					cErr <- fmt.Errorf("invalid nucleotide in fasta file (\"%s\")", string(line[i]))
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

// ReadEncodeScoreAlignment reads an alignment in fasta format to a channel
// of EncodedFastaRecord structs - converting the nucleotide sequence to EP's bitwise coding scheme
// and additionally scoring each record and getting ATGC counts
func ReadEncodeScoreAlignment(f io.Reader, hardGaps bool, chnl chan EncodedFastaRecord, cErr chan error, cDone chan bool) {

	var err error

	var coding [256]byte
	switch hardGaps {
	case true:
		coding = encoding.MakeEncodingArrayHardGaps()
	case false:
		coding = encoding.MakeEncodingArray()
	}

	scoring := encoding.MakeEncodedScoreArray()

	s := bufio.NewScanner(f)
	s.Buffer(make([]byte, 0), 1024*1024)

	first := true

	var id string
	var description string
	var seqBuffer []byte
	var line []byte
	var nuc byte
	var width int
	var score int64

	var fr EncodedFastaRecord
	var counting [256]int

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

			fr = EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter, Score: score}
			fr.Count_A = counting[136]
			fr.Count_T = counting[24]
			fr.Count_G = counting[72]
			fr.Count_C = counting[40]
			chnl <- fr
			counter++

			description = string(line[1:])
			id = strings.Fields(description)[0]
			seqBuffer = make([]byte, 0)
			score = 0
			for i := range counting {
				counting[i] = 0
			}

		} else {
			encodedLine := make([]byte, len(line))
			for i := range line {
				nuc = coding[line[i]]
				if nuc == 0 {
					cErr <- fmt.Errorf("invalid nucleotide in fasta file (\"%s\")", string(line[i]))
					return
				}
				encodedLine[i] = nuc
				score += scoring[nuc]
				counting[nuc]++
			}
			seqBuffer = append(seqBuffer, encodedLine...)
		}
	}

	if len(seqBuffer) > 0 {
		if counter > 0 && len(seqBuffer) != width {
			cErr <- errors.New("different length sequences in input file: is this an alignment?")
			return
		}
		fr = EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter, Score: score}
		fr.Count_A = counting[136]
		fr.Count_T = counting[24]
		fr.Count_G = counting[72]
		fr.Count_C = counting[40]
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

// ReadEncodeAlignmentToList is as ReadEncodeAlignment but returns a slice of EncodedFastaRecords instead
// of passing each EncodedFastaRecord down a channel
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
	s.Buffer(make([]byte, 0), 1024*1024)

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
					return []EncodedFastaRecord{}, fmt.Errorf("invalid nucleotide in fasta file (\"%s\")", string(line[i]))
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

// WriteAlignment reads FastaRecords from a channel and writes them to file or stdout,
// in the order in which they are present in the input file.
// It passes a true to a done channel when the channel of fasta records is empty
func WriteAlignment(ch chan FastaRecord, w io.Writer, cdone chan bool, cerr chan error) {

	outputMap := make(map[int]FastaRecord)

	counter := 0

	var err error

	for FR := range ch {

		outputMap[FR.Idx] = FR

		for {
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
				break
			}
		}

	}

	cdone <- true
}

// TO DO - test
func WriteWrapAlignment(ch chan FastaRecord, w io.Writer, wrap int, cdone chan bool, cerr chan error) {

	outputMap := make(map[int]FastaRecord)

	var (
		counter, written int
		err              error
	)

	for FR := range ch {

		outputMap[FR.Idx] = FR

		for {
			if fastarecord, ok := outputMap[counter]; ok {
				_, err = w.Write([]byte(">" + fastarecord.ID + "\n"))
				if err != nil {
					cerr <- err
				}
				for {
					if written < len(fastarecord.Seq) {
						if written+wrap >= len(fastarecord.Seq) {
							_, err = w.Write([]byte(fastarecord.Seq[written:] + "\n"))
							if err != nil {
								cerr <- err
							}
							written = written + wrap
						} else {
							_, err = w.Write([]byte(fastarecord.Seq[written:written+wrap] + "\n"))
							if err != nil {
								cerr <- err
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

	cdone <- true
}
