package fasta

import (
	"reflect"
	"testing"
)

func TestEncodeDecode(t *testing.T) {
	FR := Record{ID: "Seq1", Description: "Seq1", Idx: 0, Seq: "ATGCRMWSKYVHDBN-?"}

	desiredResult := EncodedRecord{ID: "Seq1", Description: "Seq1", Idx: 0, Seq: []byte{136, 24, 72, 40, 192, 160, 144, 96, 80, 48, 224, 176, 208, 112, 240, 244, 242}}
	actualResult, err := FR.Encode()
	if err != nil {
		t.Error(err)
	}
	if !reflect.DeepEqual(desiredResult, actualResult) {
		t.Errorf("problem in TestEncodeDecode - Encode")
	}

	if !reflect.DeepEqual(FR, actualResult.Decode()) {
		t.Errorf("problem in TestEncodeDecode - Decode")
	}

	FR = Record{ID: "Seq1", Description: "Seq1", Idx: 0, Seq: "atgcrmwskyvhdbn"}
	_, err = FR.Encode()
	if err != nil {
		t.Error(err)
	}
}

func TestEncodeError(t *testing.T) {
	FR := Record{ID: "Seq1", Description: "Seq1", Idx: 0, Seq: "ATGCRMWSKYVHDBN-?J"}
	_, err := FR.Encode()
	if err.Error() != "invalid nucleotide in fasta record (\"J\")" {
		t.Error(err)
	}

	FR = Record{ID: "Seq1", Description: "Seq1", Idx: 0, Seq: "ATGCRMWSKYVHDBN-?Å"}
	_, err = FR.Encode()
	if err.Error() != "invalid nucleotide in fasta record (\"Å\")" {
		t.Error(err)
	}
}

func TestDegap(t *testing.T) {
	in := Record{Seq: "ACGT-ACGT-ACGT"}
	out := Record{Seq: "ACGTACGTACGT"}

	if !reflect.DeepEqual(in.Degap(), out) {
		t.Errorf("problem in TestDegap()")
	}
}

func TestRecordComplement(t *testing.T) {
	inn := Record{Seq: "ACGTRYSWKMBDHVN-?acgtryswkmbdhvn-?"}
	out := Record{Seq: "TGCAYRSWMKVHDBN-?tgcayrswmkvhdbn-?"}

	if !reflect.DeepEqual(inn.Complement(), out) {
		t.Errorf("problem in TestFRComplement()")
	}
}

func TestRecordReverseComplement(t *testing.T) {
	inn := Record{Seq: "ACGTRYSWKMBDHVN-?acgtryswkmbdhvn-?"}
	out := Record{Seq: "?-nbdhvkmwsryacgt?-NBDHVKMWSRYACGT"}

	if !reflect.DeepEqual(inn.ReverseComplement(), out) {
		t.Errorf("problem in TestFRReverseComplement()")
	}
}

func TestEncodedRecordComplement(t *testing.T) {
	inn, err := Record{Seq: "ACGTRYSWKMBDHVN-?acgtryswkmbdhvn-?"}.Encode()
	if err != nil {
		t.Error(err)
	}
	out, err := Record{Seq: "TGCAYRSWMKVHDBN-?TGCAYRSWMKVHDBN-?"}.Encode()
	if err != nil {
		t.Error(err)
	}

	if !reflect.DeepEqual(inn.Complement(), out) {
		t.Errorf("problem in TestEFRComplement()")
	}
}

func TestEncodedRecordReverseComplement(t *testing.T) {
	inn, err := Record{Seq: "ACGTRYSWKMBDHVN-?acgtryswkmbdhvn-?"}.Encode()
	if err != nil {
		t.Error(err)
	}
	out, err := Record{Seq: "?-NBDHVKMWSRYACGT?-NBDHVKMWSRYACGT"}.Encode()
	if err != nil {
		t.Error(err)
	}

	if !reflect.DeepEqual(inn.ReverseComplement(), out) {
		t.Errorf("problem in TestEFRReverseComplement()")
	}
}

func TestCalculateBaseContent(t *testing.T) {
	FR := Record{ID: "Seq1", Description: "Seq1", Idx: 0, Seq: "ATGCATGATA"}
	EFR, err := FR.Encode()
	if err != nil {
		t.Error(err)
	}
	EFR.CalculateBaseContent()

	if EFR.Count_A != 4 {
		t.Errorf("problem in TestCalculateBaseContent()")
	}

	if EFR.Count_T != 3 {
		t.Errorf("problem in TestCalculateBaseContent()")
	}

	if EFR.Count_G != 2 {
		t.Errorf("problem in TestCalculateBaseContent()")
	}

	if EFR.Count_C != 1 {
		t.Errorf("problem in TestCalculateBaseContent()")
	}
}

func TestScore(t *testing.T) {
	FR := Record{ID: "Seq1", Description: "Seq1", Idx: 0, Seq: "ATGCATGATA"}
	EFR, err := FR.Encode()
	if err != nil {
		t.Error(err)
	}
	EFR.CalculateCompleteness()

	if EFR.Score != 120 {
		t.Errorf("problem in TestScore() (1)")
	}

	FR = Record{ID: "Seq1", Description: "Seq1", Idx: 0, Seq: "----"}
	EFR, err = FR.Encode()
	if err != nil {
		t.Error(err)
	}
	EFR.CalculateCompleteness()

	if EFR.Score != 12 {
		t.Errorf("problem in TestScore() (2)")
	}

	FR = Record{ID: "Seq1", Description: "Seq1", Idx: 0, Seq: "NNnn----"}
	EFR, err = FR.Encode()
	if err != nil {
		t.Error(err)
	}
	EFR.CalculateCompleteness()

	if EFR.Score != 24 {
		t.Errorf("problem in TestScore() (3)")
	}

	FR = Record{ID: "Seq1", Description: "Seq1", Idx: 0, Seq: "NNnn----ATGC"}
	EFR, err = FR.Encode()
	if err != nil {
		t.Error(err)
	}
	EFR.CalculateCompleteness()

	if EFR.Score != 72 {
		t.Errorf("problem in TestScore() (4)")
	}
}
