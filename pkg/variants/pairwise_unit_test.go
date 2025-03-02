package variants

import (
	"fmt"
	"reflect"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/fasta"
)

func TestGetIndelsPair(t *testing.T) {

	ref, err := fasta.Record{Seq: "ATG---ATGATGAT"}.Encode()
	if err != nil {
		t.Error(err)
	}
	refSeq := ref.Seq
	que, err := fasta.Record{Seq: "ATGATGAT--TG--"}.Encode()
	if err != nil {
		t.Error(err)
	}
	queSeq := que.Seq

	offsetRefCoord, offsetMSACoord := GetMSAOffsets(refSeq)

	indels := getIndelsPair(refSeq, queSeq, offsetRefCoord, offsetMSACoord)

	desiredResultV := []Variant{
		Variant{Position: 3, Changetype: "ins", Length: 3},
		Variant{Position: 6, Changetype: "del", Length: 2},
	}

	if !reflect.DeepEqual(desiredResultV, indels) {
		t.Errorf("Problem in TestGetIndelsPair()")
		fmt.Println(indels)
	}

	s := make([]string, 0)

	for _, v := range indels {
		temp, _ := FormatVariant(v, false)
		s = append(s, temp)
	}

	desiredResultS := []string{"ins:3:3", "del:6:2"}

	if !reflect.DeepEqual(desiredResultS, s) {
		t.Errorf("Problem in TestGetIndelsPair()")
		fmt.Println(s)
	}
}

func TestGetNucsPair(t *testing.T) {

	ref, err := fasta.Record{Seq: "ATGA-TGACC"}.Encode()
	if err != nil {
		t.Error(err)
	}
	refSeq := ref.Seq

	que, err := fasta.Record{Seq: "TTGA-TGSCS"}.Encode()
	if err != nil {
		t.Error(err)
	}
	queSeq := que.Seq

	offsetRefCoord, offsetMSACoord := GetMSAOffsets(refSeq)

	nucs := getNucsPair(refSeq, queSeq, []int{1, 2, 3, 4, 5, 6, 7, 8, 9}, offsetRefCoord, offsetMSACoord)

	desiredResultV := []Variant{
		Variant{RefAl: "A", QueAl: "T", Position: 1, Changetype: "nuc"},
		Variant{RefAl: "A", QueAl: "S", Position: 7, Changetype: "nuc"},
	}

	if !reflect.DeepEqual(desiredResultV, nucs) {
		t.Errorf("Problem in TestGetNucsPair()")
		fmt.Println(nucs)
	}

	s := make([]string, 0)

	for _, v := range nucs {
		temp, _ := FormatVariant(v, false)
		s = append(s, temp)
	}

	desiredResultS := []string{"nuc:A1T", "nuc:A7S"}

	if !reflect.DeepEqual(desiredResultS, s) {
		t.Errorf("Problem in TestGetNucsPair()")
		fmt.Println(s)
	}
}

func TestGetAAsPair(t *testing.T) {

	ref, err := fasta.Record{Seq: "ATG-TCTAGACCC"}.Encode()
	if err != nil {
		t.Error(err)
	}
	refSeq := ref.Seq

	que, err := fasta.Record{Seq: "ATGATGTAGAAAA"}.Encode()
	if err != nil {
		t.Error(err)
	}
	queSeq := que.Seq

	r := Region{Whichtype: "protein-coding", Name: "nspX", Start: 1, Stop: 12, Translation: "MSRP", Strand: 1, Positions: []int{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}}

	offsetRefCoord, offsetMSACoord := GetMSAOffsets(refSeq)

	AAs := getAAsPair(refSeq, queSeq, r, offsetRefCoord, offsetMSACoord)

	desiredResultV := []Variant{
		Variant{RefAl: "S", QueAl: "C", Position: 4, Changetype: "aa", SNPs: "nuc:C5G", Residue: 2, Feature: "nspX"},
		Variant{RefAl: "P", QueAl: "K", Position: 10, Changetype: "aa", SNPs: "nuc:C10A;nuc:C11A;nuc:C12A", Residue: 4, Feature: "nspX"},
	}

	if !reflect.DeepEqual(desiredResultV, AAs) {
		t.Errorf("Problem in TestGetAAsPair()")
		fmt.Println(AAs)
	}

	s := make([]string, 0)

	for _, v := range AAs {
		temp, _ := FormatVariant(v, false)
		s = append(s, temp)
	}

	desiredResultS := []string{"aa:nspX:S2C", "aa:nspX:P4K"}

	if !reflect.DeepEqual(desiredResultS, s) {
		t.Errorf("Problem in TestGetAAsPair()")
		fmt.Println(s)
	}

	s = make([]string, 0)

	for _, v := range AAs {
		temp, _ := FormatVariant(v, true)
		s = append(s, temp)
	}

	desiredResultS = []string{"aa:nspX:S2C(nuc:C5G)", "aa:nspX:P4K(nuc:C10A;nuc:C11A;nuc:C12A)"}

	if !reflect.DeepEqual(desiredResultS, s) {
		t.Errorf("Problem in TestGetAAsPair()")
		fmt.Println(s)
	}
}
