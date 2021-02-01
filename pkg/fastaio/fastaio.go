package fastaio

import (
	"bufio"
	"errors"
	"fmt"
	"github.com/cov-ert/gofasta/pkg/encoding"
	"os"
	"strings"
)

// FastaRecord is a simple struct for Fasta records
type FastaRecord struct {
	ID          string
	Description string
	Seq         string
	Idx	    int
}

// EncodedFastaRecord is a struct for one Fasta record
type EncodedFastaRecord struct {
	ID          string
	Description string
	Seq         []byte
	Idx int
}

func getAlignmentDims(infile string) (int, int, error) {
	n := 0
	l := 0

	f, err := os.Open(infile)
	if err != nil {
		return 0, 0, err
	}
	defer f.Close()

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

	err = s.Err()

	if err != nil {
		return 0, 0, err
	}

	return n, l, err
}

// ReadAlignment reads an alignment in fasta format to a channel
// of FastaRecord structs
func ReadAlignment(infile string, chnl chan FastaRecord, chnlerr chan error, cdone chan bool) {

	f, err := os.Open(infile)

	if err != nil {
		chnlerr <- err
	}

	defer f.Close()

	s := bufio.NewScanner(f)

	first := true

	var id string
	var description string
	var seqBuffer string

	for s.Scan() {
		line := string(s.Text())

		if first {

			if string(line[0]) != ">" {
				chnlerr <- errors.New("badly formatted fasta file")
			}

			description = line[1:]
			id = strings.Fields(description)[0]

			first = false

		} else if string(line[0]) == ">" {

			fr := FastaRecord{ID: id, Description: description, Seq: seqBuffer}
			chnl <- fr

			description = line[1:]
			id = strings.Fields(description)[0]
			seqBuffer = ""

		} else {
			seqBuffer = seqBuffer + strings.ToUpper(line)
		}

	}

	fr := FastaRecord{ID: id, Description: description, Seq: seqBuffer}
	chnl <- fr

	err = s.Err()
	if err != nil {
		chnlerr <- err
	}

	cdone <- true
}

// PopulateByteArrayGetNames reads a fasta format alignment into an array of
// uint8 representations of nucleotides, and gets the names of the fasta records
// at the same time.
func PopulateByteArrayGetNames(infile string) ([][]uint8, []string, error) {
	height, width, err := getAlignmentDims(infile)

	if err != nil {
		return [][]uint8{}, []string{}, err
	}

	IDs := make([]string, 0)

	A := make([][]uint8, height)
	for i := 0; i < height; i++ {
		A[i] = make([]uint8, width)
	}

	byteDict := encoding.MakeByteDict()

	c := make(chan FastaRecord)
	cerr := make(chan error)
	cdone := make(chan bool)

	go ReadAlignment(infile, c, cerr, cdone)

	i := 0

	for n := 1; n > 0; {
		select {

		case err := <-cerr:
			return [][]uint8{}, []string{}, err

		case record := <-c:

			IDs = append(IDs, record.ID)

			if len(record.Seq) != width {
				return [][]uint8{}, []string{}, errors.New("different length sequences in input file: is this an alignment?")
			}

			for j, nuc := range record.Seq {

				if val, ok := byteDict[nuc]; ok {
					A[i][j] = val

				} else {
					return [][]uint8{}, []string{}, fmt.Errorf("invalid nucleotide in fasta file (%s)", string(nuc))
				}
			}

			i++

		case <-cdone:
			n--
		}
	}

	return A, IDs, nil
}


// ReadEncodeAlignment reads an alignment in fasta format to a channel
// of encodedFastaRecord structs - converting sequence to EP's bitwise coding scheme
func ReadEncodeAlignment(inFile string, chnl chan EncodedFastaRecord, cErr chan error, cDone chan bool) {

	var err error
	var f *os.File

	if inFile != "stdin" {
		f, err = os.Open(inFile)
		if err != nil {
			cErr <- err
		}
	} else {
		f = os.Stdin
	}

	defer f.Close()

	encoding := encoding.MakeEncodingArray()

	s := bufio.NewScanner(f)

	first := true

	var id string
	var description string
	var seqBuffer []byte
	var line []byte
	var nuc byte

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

			fr := EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
			chnl <- fr
			counter++

			description = string(line[1:])
			id = strings.Fields(description)[0]
			seqBuffer = make([]byte, 0)

		} else {
			encodedLine := make([]byte, len(line))
			for i := range(line) {
				nuc = encoding[line[i]]
				if nuc == 0 {
					cErr<- fmt.Errorf("invalid nucleotide in fasta file (%s)", string(nuc))
				}
				encodedLine[i] = nuc
			}
			seqBuffer = append(seqBuffer, encodedLine...)
		}
	}

	fr := EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
	chnl <- fr

	err = s.Err()
	if err != nil {
		cErr <- err
	}

	cDone <- true
}
