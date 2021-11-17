package fastaio

import (
	"bytes"
	"testing"

	"github.com/cov-ert/gofasta/pkg/encoding"
)

func TestGetAlignmentDims(t *testing.T) {
	alignmentData := []byte(
		`>Target1
ATGATC
>Target2
ATGATG
>Target3
ATTTTC
`)

	alignment := bytes.NewReader(alignmentData)

	n, l, err := getAlignmentDims(alignment)
	if err != nil {
		t.Error(err)
	}

	if n != 3 {
		t.Errorf("wrong number of sequences in TestGetAlignmentDims()")
	}

	if l != 6 {
		t.Errorf("wrong number of sequences in TestGetAlignmentDims()")
	}
}

func TestReadAlignment(t *testing.T) {
	alignmentData := []byte(
		`>Target1
ATGATC
>Target2
ATGATG
>Target3
ATTTTC
`)

	alignment := bytes.NewReader(alignmentData)

	cErr := make(chan error)
	cFR := make(chan FastaRecord)
	cReadDone := make(chan bool)

	go ReadAlignment(alignment, cFR, cErr, cReadDone)

	out := new(bytes.Buffer)
	cWriteDone := make(chan bool)

	go WriteAlignment(cFR, out, cWriteDone, cErr)

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			t.Error(err)
		case <-cReadDone:
			close(cFR)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			t.Error(err)
		case <-cWriteDone:
			n--
		}
	}

	if string(out.Bytes()) != `>Target1
ATGATC
>Target2
ATGATG
>Target3
ATTTTC
` {
		t.Errorf("problem in TestReadAlignment()")
	}

}

func TestReadEncodeAlignment(t *testing.T) {
	alignmentData := []byte(
		`>Target1
ATGATC
>Target2
ATGATG
>Target3
ATTTTC
`)

	alignment := bytes.NewReader(alignmentData)

	cErr := make(chan error)
	cFR := make(chan EncodedFastaRecord)
	cReadDone := make(chan bool)

	out := new(bytes.Buffer)

	go ReadEncodeAlignment(alignment, false, cFR, cErr, cReadDone)

	go func() {
		for n := 1; n > 0; {
			select {
			case err := <-cErr:
				t.Error(err)
			case <-cReadDone:
				close(cFR)
			}
		}
	}()

	DA := encoding.MakeDecodingArray()

	for FR := range cFR {
		out.Write([]byte(">" + FR.ID + "\n"))
		for _, nuc := range FR.Seq {
			out.Write([]byte(DA[nuc]))
		}
		out.Write([]byte("\n"))
	}

	if string(out.Bytes()) != `>Target1
ATGATC
>Target2
ATGATG
>Target3
ATTTTC
` {
		t.Errorf("problem in TestReadEncodeAlignment()")
	}

}

func TestReadEncodeAlignmentToList(t *testing.T) {
	alignmentData := []byte(
		`>Target1
ATGATC
>Target2
ATGATG
>Target3
ATTTTC
`)

	alignment := bytes.NewReader(alignmentData)
	out := new(bytes.Buffer)

	EFRs, err := ReadEncodeAlignmentToList(alignment, false)
	if err != nil {
		t.Error(err)
	}

	DA := encoding.MakeDecodingArray()

	for _, FR := range EFRs {
		out.Write([]byte(">" + FR.ID + "\n"))
		for _, nuc := range FR.Seq {
			out.Write([]byte(DA[nuc]))
		}
		out.Write([]byte("\n"))
	}

	if string(out.Bytes()) != `>Target1
ATGATC
>Target2
ATGATG
>Target3
ATTTTC
` {
		t.Errorf("problem in TestReadEncodeAlignmentToList()")
	}

}

// func TestReadEncodeScoreAlignment(t *testing.T) {
// 	alignmentData := []byte(
// 		`>Target1
// ATGATN
// >Target2
// ATGANN
// >Target3
// ATTR-N
// `)

// 	alignment := bytes.NewReader(alignmentData)

// 	cErr := make(chan error)
// 	cFR := make(chan EncodedFastaRecord)
// 	cReadDone := make(chan bool)

// 	alOut := new(bytes.Buffer)
// 	scoreOut := make([]int64, 0)

// 	go ReadEncodeScoreAlignment(alignment, cFR, cErr, cReadDone)

// 	go func() {
// 		for n := 1; n > 0; {
// 			select {
// 			case err := <-cErr:
// 				t.Error(err)
// 			case <-cReadDone:
// 				close(cFR)
// 			}
// 		}
// 	}()

// 	DA := encoding.MakeDecodingArray()

// 	for FR := range cFR {
// 		alOut.Write([]byte(">" + FR.ID + "\n"))
// 		for _, nuc := range FR.Seq {
// 			alOut.Write([]byte(DA[nuc]))
// 		}
// 		alOut.Write([]byte("\n"))

// 		scoreOut = append(scoreOut, FR.Score)
// 	}

// 	if string(alOut.Bytes()) != `>Target1
// ATGATN
// >Target2
// ATGANN
// >Target3
// ATTR-N
// ` {
// 		t.Errorf("problem in TestReadEncodeAlignment() (alignment)")
// 	}

// 	test := true
// 	for i := range scoreOut {
// 		if scoreOut[i] != []int64{63, 54, 48}[i] {
// 			test = false
// 			break
// 		}
// 	}
// 	if !test {
// 		t.Errorf("problem in TestReadEncodeAlignment() (score)")
// 	}
// }
