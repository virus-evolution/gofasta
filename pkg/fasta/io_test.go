package fasta

import (
	"bytes"
	"fmt"
	"reflect"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/encoding"
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

func TestStreamAlignment(t *testing.T) {
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
	cFR := make(chan Record)
	cReadDone := make(chan bool)

	go StreamAlignment(alignment, cFR, cErr, cReadDone)

	out := new(bytes.Buffer)
	cWriteDone := make(chan bool)

	go WriteAlignment(cFR, out, cErr, cWriteDone)

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
		t.Errorf("problem in TestStreamAlignment()")
		fmt.Print(string(out.Bytes()))
	}

}

func TestStreamAlignmentShortSeqs(t *testing.T) {
	alignmentData := []byte(`>Target1
ATGATC
>TargetShort1
ATTAT
>Target2
ATGATC
`)

	alignment := bytes.NewReader(alignmentData)

	cErr := make(chan error)
	cFR := make(chan Record, 10)
	cReadDone := make(chan bool)

	go StreamAlignment(alignment, cFR, cErr, cReadDone)

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			if err.Error() != "different length sequences in input file: is this an alignment?" {
				t.Error(err)
			}
			n--
		case <-cReadDone:
			close(cFR)
			t.Error("reached end of alignment without throwing a sequence length error")
			n--
		}
	}
}

func TestStreamEncodeAlignment(t *testing.T) {
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
	cFR := make(chan EncodedRecord)
	cReadDone := make(chan bool)

	out := new(bytes.Buffer)

	go StreamEncodeAlignment(alignment, cFR, cErr, cReadDone, false, false, false)

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
		t.Errorf("problem in TestStreamEncodeAlignment()")
	}

}

func TestStreamEncodeAlignmentShortSeqs(t *testing.T) {
	alignmentData := []byte(`>Target1
ATGATC
>TargetShort1
ATTAT
>Target2
ATGATC
`)

	alignment := bytes.NewReader(alignmentData)

	cErr := make(chan error)
	cFR := make(chan EncodedRecord, 10)
	cReadDone := make(chan bool)

	go StreamEncodeAlignment(alignment, cFR, cErr, cReadDone, false, false, false)

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			if err.Error() != "different length sequences in input file: is this an alignment?" {
				t.Error(err)
			}
			n--
		case <-cReadDone:
			close(cFR)
			t.Error("reached end of alignment without throwing a sequence length error")
			n--
		}
	}
}

func TestLoadEncodeAlignment(t *testing.T) {
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

	EFRs, err := LoadEncodeAlignment(alignment, false, false, false)
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
		t.Errorf("problem in TestStreamEncodeAlignmentToList()")
		fmt.Print(string(out.Bytes()))
	}
}

func TestStreamEncodeAlignmentToListShortSeqs(t *testing.T) {
	alignmentData := []byte(
		`>Target1
ATGATC
>Target2
ATGAT
>Target3
ATTTTC
`)

	alignment := bytes.NewReader(alignmentData)

	_, err := LoadEncodeAlignment(alignment, false, false, false)
	if err.Error() != "different length sequences in input file: is this an alignment?" {
		t.Error(err)
	}
}

func TestReadEncodeScoreAlignment(t *testing.T) {
	alignmentData := []byte(
		`>Target1
ATGATN
>Target2
ATGANN
>Target3
ATTR-N
`)

	alignment := bytes.NewReader(alignmentData)

	cErr := make(chan error)
	cFR := make(chan EncodedRecord)
	cReadDone := make(chan bool)

	alOut := new(bytes.Buffer)
	scoreOut := make([]int64, 0)

	go StreamEncodeAlignment(alignment, cFR, cErr, cReadDone, false, true, true)

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

	FRs := make([]EncodedRecord, 0)

	for FR := range cFR {
		FRs = append(FRs, FR)
		alOut.Write([]byte(">" + FR.ID + "\n"))
		for _, nuc := range FR.Seq {
			alOut.Write([]byte(DA[nuc]))
		}
		alOut.Write([]byte("\n"))

		scoreOut = append(scoreOut, FR.Score)
	}

	if string(alOut.Bytes()) != `>Target1
ATGATN
>Target2
ATGANN
>Target3
ATTR-N
` {
		t.Errorf("problem in TestReadEncodeScoreAlignment() (alignment)")
	}

	test := true
	for i := range scoreOut {
		if scoreOut[i] != []int64{63, 54, 48}[i] {
			test = false
			break
		}
	}
	if !test {
		t.Errorf("problem in TestReadEncodeScoreAlignment() (score)")
	}

	counts := [][]int{{2, 2, 1, 0}, {2, 1, 1, 0}, {1, 2, 0, 0}}
	for i := range FRs {
		for j := 0; j < 4; j++ {
			var count int
			switch j {
			case 0:
				count = FRs[i].Count_A
			case 1:
				count = FRs[i].Count_T
			case 2:
				count = FRs[i].Count_G
			case 3:
				count = FRs[i].Count_C
			}
			if counts[i][j] != count {
				t.Errorf("problem in TestReadEncodeScoreAlignment() (counts)")
			}
		}
	}
}

func TestWriteWrapAlignment(t *testing.T) {
	source := []Record{
		Record{ID: "Seq3", Seq: "ATGATGATG", Idx: 2},
		Record{ID: "Seq2", Seq: "ATGATGATG", Idx: 1},
		Record{ID: "Seq1", Seq: "ATGATGATG", Idx: 0},
	}

	sink := bytes.NewBuffer(make([]byte, 0))

	cFR := make(chan Record)
	cErr := make(chan error)
	cDone := make(chan bool)

	go WriteWrapAlignment(cFR, sink, 80, cErr, cDone)

	for _, record := range source {
		cFR <- record
	}
	close(cFR)

	select {
	case err := <-cErr:
		t.Error(err)
	case <-cDone:
	}

	desiredResult := []byte(`>Seq1
ATGATGATG
>Seq2
ATGATGATG
>Seq3
ATGATGATG
`)

	if !reflect.DeepEqual(sink.Bytes(), desiredResult) {
		t.Errorf("Problem in TestWriteWrapAlignment()")
	}

	sink.Reset()
	cFR = make(chan Record)

	go WriteWrapAlignment(cFR, sink, 5, cErr, cDone)

	for _, record := range source {
		cFR <- record
	}
	close(cFR)

	select {
	case err := <-cErr:
		t.Error(err)
	case <-cDone:
	}

	desiredResult = []byte(`>Seq1
ATGAT
GATG
>Seq2
ATGAT
GATG
>Seq3
ATGAT
GATG
`)

	if !reflect.DeepEqual(sink.Bytes(), desiredResult) {
		t.Errorf("Problem in TestWriteWrapAlignment()")
	}
}
