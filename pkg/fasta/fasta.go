package fasta

import (
	"fmt"

	"github.com/virus-evolution/gofasta/pkg/alphabet"
	"github.com/virus-evolution/gofasta/pkg/encoding"
)

// A struct for one Fasta record
type Record struct {
	ID          string
	Description string
	Seq         string
	Idx         int
}

// A struct for one Fasta record whose sequence is encoded using EP's scheme
type EncodedRecord struct {
	ID          string
	Description string
	Seq         []byte
	Idx         int
	Score       int64 // this is for genome completeness
	Count_A     int
	Count_T     int
	Count_G     int
	Count_C     int
}

func (FR Record) encode(hardGaps bool) (EncodedRecord, error) {
	var EA [256]byte
	if hardGaps {
		EA = encoding.MakeEncodingArrayHardGaps()
	} else {
		EA = encoding.MakeEncodingArray()
	}
	EFR := EncodedRecord{ID: FR.ID, Description: FR.Description, Idx: FR.Idx}
	seq := make([]byte, len(FR.Seq))
	for i, nuc := range FR.Seq {
		if nuc > 127 || EA[nuc] == 0 {
			return EncodedRecord{}, fmt.Errorf("invalid nucleotide in fasta record (\"%c\")", nuc)
		}
		seq[i] = EA[nuc]
	}
	EFR.Seq = seq
	return EFR, nil
}

// Convert a Record to an EncodedRecord
func (FR Record) Encode() (EncodedRecord, error) {
	return FR.encode(false)
}

// Convert a Record to an EncodedRecord with hard gaps
func (FR Record) EncodeHardGaps() (EncodedRecord, error) {
	return FR.encode(true)
}

// Strip the gaps from a Record's sequence, returning a new Record
func (FR Record) Degap() Record {
	NFR := Record{ID: FR.ID, Description: FR.Description, Idx: FR.Idx}
	t := ""
	for _, char := range FR.Seq {
		if char != '-' {
			t = t + string(char)
		}
	}
	NFR.Seq = t
	return NFR
}

// Complement a Record's sequence, returning a new Record
func (FR Record) Complement() Record {
	NFR := Record{ID: FR.ID, Description: FR.Description, Idx: FR.Idx}
	CA := alphabet.MakeCompArray()
	ba := make([]byte, len(FR.Seq))
	for i := 0; i < len(FR.Seq); i++ {
		ba[i] = CA[FR.Seq[i]]
	}
	NFR.Seq = string(ba)
	return NFR
}

// Reverse complement a Record's sequence, returning a new Record
func (FR Record) ReverseComplement() Record {
	NFR := FR.Complement()
	temp := []byte(NFR.Seq)
	for i, j := 0, len(temp)-1; i < j; i, j = i+1, j-1 {
		temp[i], temp[j] = temp[j], temp[i]
	}
	NFR.Seq = string(temp)
	return NFR
}

// Convert an EncodedRecord to a Record, returning a new Record
func (EFR EncodedRecord) Decode() Record {
	FR := Record{ID: EFR.ID, Description: EFR.Description, Idx: EFR.Idx}
	DA := encoding.MakeDecodingArray()
	seq := ""
	for _, nuc := range EFR.Seq {
		seq = seq + DA[nuc]
	}
	FR.Seq = seq
	return FR
}

// Calculate the ATGC content of an EncodedRecord in place
func (EFR *EncodedRecord) CalculateBaseContent() {
	var counting [256]int
	for _, nuc := range EFR.Seq {
		counting[nuc]++
	}
	EFR.Count_A = counting[136]
	EFR.Count_T = counting[24]
	EFR.Count_G = counting[72]
	EFR.Count_C = counting[40]
}

// Score an EncodedRecord for completeness in place
func (EFR *EncodedRecord) CalculateCompleteness() {
	var score int64
	scoring := encoding.MakeEncodedScoreArray()
	for _, nuc := range EFR.Seq {
		score += scoring[nuc]
	}
	EFR.Score = score
}

// Complement an EncodedRecord's sequence, returning a new EncodedRecord
func (EFR EncodedRecord) Complement() EncodedRecord {
	NFR := EncodedRecord{ID: EFR.ID, Description: EFR.Description, Idx: EFR.Idx}
	CA := alphabet.MakeEncodedCompArray()
	NFR.Seq = make([]byte, len(EFR.Seq))
	for i := 0; i < len(EFR.Seq); i++ {
		NFR.Seq[i] = CA[EFR.Seq[i]]
	}
	return NFR
}

// Reverse complement an EncodedRecord's sequence, returning a new EncodedRecord
func (EFR EncodedRecord) ReverseComplement() EncodedRecord {
	NEFR := EFR.Complement()
	for i, j := 0, len(NEFR.Seq)-1; i < j; i, j = i+1, j-1 {
		NEFR.Seq[i], NEFR.Seq[j] = NEFR.Seq[j], NEFR.Seq[i]
	}
	return NEFR
}
