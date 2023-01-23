package fastaio

import (
	"bytes"
	"fmt"
	"reflect"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/encoding"
)

func TestEncodeDecode(t *testing.T) {
	FR := FastaRecord{ID: "Seq1", Description: "Seq1", Idx: 0, Seq: "ATGCRMWSKYVHDBN-?"}

	desiredResult := EncodedFastaRecord{ID: "Seq1", Description: "Seq1", Idx: 0, Seq: []byte{136, 24, 72, 40, 192, 160, 144, 96, 80, 48, 224, 176, 208, 112, 240, 244, 242}}

	if !reflect.DeepEqual(desiredResult, FR.Encode()) {
		t.Errorf("problem in TestEncodeDecode - Encode")
	}

	if !reflect.DeepEqual(FR, FR.Encode().Decode()) {
		t.Errorf("problem in TestEncodeDecode - Decode")
	}
}

func TestDegap(t *testing.T) {
	in := FastaRecord{Seq: "ACGT-ACGT-ACGT"}
	out := FastaRecord{Seq: "ACGTACGTACGT"}

	if !reflect.DeepEqual(in.Degap(), out) {
		t.Errorf("Problem in TestDegap()")
	}
}

func TestFRComplement(t *testing.T) {
	inn := FastaRecord{Seq: "ACGTRYSWKMBDHVN-?acgtryswkmbdhvn-?"}
	out := FastaRecord{Seq: "TGCAYRSWMKVHDBN-?tgcayrswmkvhdbn-?"}

	if !reflect.DeepEqual(inn.Complement(), out) {
		t.Errorf("Problem in TestFRComplement()")
	}
}

func TestFRReverseComplement(t *testing.T) {
	inn := FastaRecord{Seq: "ACGTRYSWKMBDHVN-?acgtryswkmbdhvn-?"}
	out := FastaRecord{Seq: "?-nbdhvkmwsryacgt?-NBDHVKMWSRYACGT"}

	if !reflect.DeepEqual(inn.ReverseComplement(), out) {
		t.Errorf("Problem in TestFRReverseComplement()")
	}
}

func TestEFRComplement(t *testing.T) {
	inn := FastaRecord{Seq: "ACGTRYSWKMBDHVN-?acgtryswkmbdhvn-?"}.Encode()
	out := FastaRecord{Seq: "TGCAYRSWMKVHDBN-?TGCAYRSWMKVHDBN-?"}.Encode()

	if !reflect.DeepEqual(inn.Complement(), out) {
		t.Errorf("Problem in TestEFRComplement()")
	}
}

func TestEFRReverseComplement(t *testing.T) {
	inn := FastaRecord{Seq: "ACGTRYSWKMBDHVN-?acgtryswkmbdhvn-?"}.Encode()
	out := FastaRecord{Seq: "?-NBDHVKMWSRYACGT?-NBDHVKMWSRYACGT"}.Encode()

	if !reflect.DeepEqual(inn.ReverseComplement(), out) {
		t.Errorf("Problem in TestEFRReverseComplement()")
	}
}

func TestCalculateBaseContent(t *testing.T) {
	FR := FastaRecord{ID: "Seq1", Description: "Seq1", Idx: 0, Seq: "ATGCATGATA"}
	EFR := FR.Encode()
	EFR.CalculateBaseContent()

	if EFR.Count_A != 4 {
		t.Errorf("probem in TestCalculateBaseContent()")
	}

	if EFR.Count_T != 3 {
		t.Errorf("probem in TestCalculateBaseContent()")
	}

	if EFR.Count_G != 2 {
		t.Errorf("probem in TestCalculateBaseContent()")
	}

	if EFR.Count_C != 1 {
		t.Errorf("probem in TestCalculateBaseContent()")
	}
}

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

func TestReadAlignmentShortSeqs(t *testing.T) {
	alignmentData := []byte(`>Target1
ATGATC
>TargetShort1
ATTAT
>Target2
ATGATC
`)

	alignment := bytes.NewReader(alignmentData)

	cErr := make(chan error)
	cFR := make(chan FastaRecord, 10)
	cReadDone := make(chan bool)

	go ReadAlignment(alignment, cFR, cErr, cReadDone)

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

func TestReadEncodeAlignmentShortSeqs(t *testing.T) {
	alignmentData := []byte(`>Target1
ATGATC
>TargetShort1
ATTAT
>Target2
ATGATC
`)

	alignment := bytes.NewReader(alignmentData)

	cErr := make(chan error)
	cFR := make(chan EncodedFastaRecord, 10)
	cReadDone := make(chan bool)

	go ReadEncodeAlignment(alignment, false, cFR, cErr, cReadDone)

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

func TestReadEncodeAlignmentToListShortSeqs(t *testing.T) {
	alignmentData := []byte(
		`>Target1
ATGATC
>Target2
ATGAT
>Target3
ATTTTC
`)

	alignment := bytes.NewReader(alignmentData)

	_, err := ReadEncodeAlignmentToList(alignment, false)
	if err.Error() != "different length sequences in input file: is this an alignment?" {
		t.Error(err)
	}
}

func TestConsensus(t *testing.T) {
	alignmentData := []byte(
		`>Target1
ATGATN
>Target2
ATGANN
>Target3
ATTR-N
`)

	alignment := bytes.NewReader(alignmentData)

	consensus, err := Consensus(alignment)
	if err != nil {
		t.Error(err)
	}

	if consensus.Seq != "ATGATN" {
		t.Errorf("problem in TestConsensus()")
		fmt.Println(consensus.Seq)
	}

	if consensus.ID != "consensus" {
		t.Errorf("problem in TestConsensus()")
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
	cFR := make(chan EncodedFastaRecord)
	cReadDone := make(chan bool)

	alOut := new(bytes.Buffer)
	scoreOut := make([]int64, 0)

	go ReadEncodeScoreAlignment(alignment, false, cFR, cErr, cReadDone)

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

	FRs := make([]EncodedFastaRecord, 0)

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
		t.Errorf("problem in TestReadEncodeAlignment() (alignment)")
	}

	test := true
	for i := range scoreOut {
		if scoreOut[i] != []int64{63, 54, 48}[i] {
			test = false
			break
		}
	}
	if !test {
		t.Errorf("problem in TestReadEncodeAlignment() (score)")
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
				t.Errorf("problem in TestReadEncodeAlignment() (counts)")
			}
		}
	}
}

func TestWriteWrapAlignment(t *testing.T) {
	source := []FastaRecord{
		FastaRecord{ID: "Seq3", Seq: "ATGATGATG", Idx: 2},
		FastaRecord{ID: "Seq2", Seq: "ATGATGATG", Idx: 1},
		FastaRecord{ID: "Seq1", Seq: "ATGATGATG", Idx: 0},
	}

	sink := bytes.NewBuffer(make([]byte, 0))

	cFR := make(chan FastaRecord)
	cErr := make(chan error)
	cDone := make(chan bool)

	go WriteWrapAlignment(cFR, sink, 80, cDone, cErr)

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
	cFR = make(chan FastaRecord)

	go WriteWrapAlignment(cFR, sink, 5, cDone, cErr)

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
