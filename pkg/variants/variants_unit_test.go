package variants

import (
	"bytes"
	"reflect"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/fastaio"
	"github.com/virus-evolution/gofasta/pkg/genbank"
)

// type Region struct {
// 	Whichtype   string // int(ergenic) or CDS
// 	Name        string // name of CDS, if it is one
// 	Start       int    // 1-based first position of region, inclusive
// 	Stop        int    // 1-based last position of region, inclusive
// 	Codonstarts []int  // a slice of the 1-based start positions of all its codons, if this region is a CDS
// 	Translation string // amino acid sequence of this region if it is a CDS
// }

// TO DO: maybe move this to the genbank package?
func TestGetRegions(t *testing.T) {
	genbankReader := bytes.NewReader(genbankDataShort)
	gb, err := genbank.ReadGenBank(genbankReader)
	if err != nil {
		t.Error(err)
	}

	regions, err := GetRegions(gb)
	if err != nil {
		t.Error(err)
	}

	var desiredResult = []Region{
		{Whichtype: "int", Start: 1, Stop: 5},
		{Whichtype: "CDS", Name: "gene1", Start: 6, Stop: 14, Translation: "MMM*", Codonstarts: []int{6, 9, 12}},
		{Whichtype: "int", Start: 15, Stop: 20},
	}

	if !reflect.DeepEqual(regions, desiredResult) {
		t.Errorf("problem in TestGetRegions")
	}
}

func TestFindReference(t *testing.T) {
	msaData := []byte(`>MN908947.3
ATGATGATG
>nottheref hey
AAAAAAAAA
>stillnottheref
CCCCCCCCC
`)

	msaReader := bytes.NewReader(msaData)

	ref, err := findReference(msaReader, "MN908947.3")
	if err != nil {
		t.Error(err)
	}

	desiredResult := fastaio.EncodedFastaRecord{ID: "MN908947.3", Description: "MN908947.3", Seq: []byte{136, 24, 72, 136, 24, 72, 136, 24, 72}, Idx: 0}
	if !reflect.DeepEqual(ref, desiredResult) {
		t.Errorf("problem in TestfindReference")
	}

	msaReader = bytes.NewReader(msaData)
	ref, err = findReference(msaReader, "nottheref")
	if err != nil {
		t.Error(err)
	}

	desiredResult = fastaio.EncodedFastaRecord{ID: "nottheref", Description: "nottheref hey", Seq: []byte{136, 136, 136, 136, 136, 136, 136, 136, 136}, Idx: 1}
	if !reflect.DeepEqual(ref, desiredResult) {
		t.Errorf("problem in TestfindReference")
	}
}

func TestGetMSAOffsets(t *testing.T) {
	msaData := []byte(`>MN908947.3
ATG-ATG-ATG
>nottheref hey
ATGAATGAATG
>stillnottheref
ATGAATGAA-G
`)

	msaReader := bytes.NewReader(msaData)

	ref, err := findReference(msaReader, "MN908947.3")
	if err != nil {
		t.Error(err)
	}

	// 1) refToMSA = the number of bases to add to convert each reference position to MSA coordinates
	// 2) MSAToRef = the number of bases to subtract to convert each MSA position to reference coordinates
	//			   - note that for MSAToRef positions that are gaps in the ref get a 0 in the MSAToRef slice
	refToMSA, MSAToRef := GetMSAOffsets(ref.Seq)

	if !reflect.DeepEqual(refToMSA, []int{0, 0, 0, 1, 1, 1, 2, 2, 2}) {
		t.Errorf("problem in TestGetMSAOffsets (refToMSA)")
	}

	if !reflect.DeepEqual(MSAToRef, []int{0, 0, 0, 0, 1, 1, 1, 0, 2, 2, 2}) {
		t.Errorf("problem in TestGetMSAOffsets (MSAToRef)")
	}
}

func TestGetVariantsPair(t *testing.T) {
	/*
		>seq1
		nuc:C2T, gene1:M1L, nuc:A19T
		>seq2
		del:6:3
		>seq3
		ins:14:1
		>seq4
		nuc:T13G, del:14:1, nuc:A15T, del:19:1, nuc:A20T
	*/
	msaData := []byte(`>reference
ACGTAATGATGATG-AAAAAA
>seq1
ATGTATTGATGATG-AAAATA
>seq2
ACGTA---ATGATG-AAAAAA
>seq3
ACGTAATGATGATGAAAAAAA
>seq4
ACGTAATGATGAG--TAAA-T
`)

	genbankReader := bytes.NewReader(genbankDataShort)
	gb, err := genbank.ReadGenBank(genbankReader)
	if err != nil {
		t.Error(err)
	}

	regions, err := GetRegions(gb)
	if err != nil {
		t.Error(err)
	}

	msaReader := bytes.NewReader(msaData)
	ref, err := findReference(msaReader, "reference")
	if err != nil {
		t.Error(err)
	}
	refToMSA, MSAToRef := GetMSAOffsets(ref.Seq)

	msaReader = bytes.NewReader(msaData)

	queries, err := fastaio.ReadEncodeAlignmentToList(msaReader, false)
	if err != nil {
		t.Error(err)
	}

	mutations, err := GetVariantsPair(ref.Seq, queries[1].Seq, "reference", queries[1].ID, 1, regions, refToMSA, MSAToRef)
	if err != nil {
		t.Error(err)
	}

	desiredResult := AnnoStructs{Queryname: "seq1", Vs: []Variant{
		{RefAl: "C", QueAl: "T", Position: 1, Changetype: "nuc"},
		{RefAl: "M", QueAl: "L", Position: 5, Residue: 0, Changetype: "aa", Feature: "gene1", SNPs: "nuc:A6T"},
		{RefAl: "A", QueAl: "T", Position: 18, Changetype: "nuc"},
	}, Idx: 1}
	if !reflect.DeepEqual(mutations, desiredResult) {
		t.Errorf("problem in TestGetVariantsPair (seq1)")
	}

	mutations, err = GetVariantsPair(ref.Seq, queries[2].Seq, "reference", queries[2].ID, 2, regions, refToMSA, MSAToRef)
	if err != nil {
		t.Error(err)
	}

	desiredResult = AnnoStructs{Queryname: "seq2", Vs: []Variant{
		{Position: 5, Length: 3, Changetype: "del"},
	}, Idx: 2}
	if !reflect.DeepEqual(mutations, desiredResult) {
		t.Errorf("problem in TestGetVariantsPair (seq2)")
	}

	mutations, err = GetVariantsPair(ref.Seq, queries[3].Seq, "reference", queries[3].ID, 3, regions, refToMSA, MSAToRef)
	if err != nil {
		t.Error(err)
	}

	desiredResult = AnnoStructs{Queryname: "seq3", Vs: []Variant{
		{Position: 14, Length: 1, Changetype: "ins"},
	}, Idx: 3}
	if !reflect.DeepEqual(mutations, desiredResult) {
		t.Errorf("problem in TestGetVariantsPair (seq3)")
	}

	mutations, err = GetVariantsPair(ref.Seq, queries[4].Seq, "reference", queries[4].ID, 4, regions, refToMSA, MSAToRef)
	if err != nil {
		t.Error(err)
	}

	desiredResult = AnnoStructs{Queryname: "seq4", Vs: []Variant{
		{Position: 12, RefAl: "T", QueAl: "G", Changetype: "nuc"},
		{Position: 13, Length: 1, Changetype: "del"},
		{Position: 14, RefAl: "A", QueAl: "T", Changetype: "nuc"},
		{Position: 18, Length: 1, Changetype: "del"},
		{Position: 19, RefAl: "A", QueAl: "T", Changetype: "nuc"},
	}, Idx: 4}
	if !reflect.DeepEqual(mutations, desiredResult) {
		t.Errorf("problem in TestGetVariantsPair (seq4)")
	}
}

func TestFormatVariant(t *testing.T) {
	msaData := []byte(`>reference
ACGTAATGATGATG-AAAAAA
>seq1
ATGTATTGATGATG-AAAATA
>seq2
ACGTA---ATGATG-AAAAAA
>seq3
ACGTAATGATGATGAAAAAAA
>seq4
ACGTAATGATGAG--TAAA-T
`)

	genbankReader := bytes.NewReader(genbankDataShort)
	gb, err := genbank.ReadGenBank(genbankReader)
	if err != nil {
		t.Error(err)
	}

	regions, err := GetRegions(gb)
	if err != nil {
		t.Error(err)
	}

	msaReader := bytes.NewReader(msaData)
	ref, err := findReference(msaReader, "reference")
	if err != nil {
		t.Error(err)
	}
	refToMSA, MSAToRef := GetMSAOffsets(ref.Seq)

	msaReader = bytes.NewReader(msaData)

	queries, err := fastaio.ReadEncodeAlignmentToList(msaReader, false)
	if err != nil {
		t.Error(err)
	}

	mutations, err := GetVariantsPair(ref.Seq, queries[1].Seq, "reference", queries[1].ID, 1, regions, refToMSA, MSAToRef)
	if err != nil {
		t.Error(err)
	}

	/*
		>seq1
		nuc:C2T, gene1:M1L, nuc:A19T
		>seq2
		del:6:3
		>seq3
		ins:14:1
		>seq4
		nuc:T13G, del:14:1, nuc:A15T, del:19:1, nuc:A20T
	*/
	desiredResult := []string{"nuc:C2T", "aa:gene1:M1L", "nuc:A19T"}
	for i, mutation := range mutations.Vs {
		mut, err := FormatVariant(mutation, false)
		if err != nil {
			t.Error(err)
		}
		if mut != desiredResult[i] {
			t.Errorf("problem in TestFormatVariant 1")
		}
	}

	desiredResult = []string{"nuc:C2T", "aa:gene1:M1L(nuc:A6T)", "nuc:A19T"}
	for i, mutation := range mutations.Vs {
		mut, err := FormatVariant(mutation, true)
		if err != nil {
			t.Error(err)
		}
		if mut != desiredResult[i] {
			t.Errorf("problem in TestFormatVariant 2")
		}
	}

	mutations, err = GetVariantsPair(ref.Seq, queries[4].Seq, "reference", queries[4].ID, 4, regions, refToMSA, MSAToRef)
	if err != nil {
		t.Error(err)
	}

	desiredResult = []string{"nuc:T13G", "del:14:1", "nuc:A15T", "del:19:1", "nuc:A20T"}
	for i, mutation := range mutations.Vs {
		mut, err := FormatVariant(mutation, false)
		if err != nil {
			t.Error(err)
		}
		if mut != desiredResult[i] {
			t.Errorf("problem in TestFormatVariant 3")
		}
	}
}

var genbankDataShort []byte

func init() {
	genbankDataShort = []byte(`LOCUS       TEST               20 bp ss-RNA     linear   VRL 21-MAR-1987
FEATURES             Location/Qualifiers
		source          1..20
						/organism="Not a real organism"
		5'UTR           1..5
		gene            6..15
						/gene="gene1"
		CDS             6..14
						/gene="gene1"
						/translation="MMM"
		3'UTR           15..20
ORIGIN      
		1 acgtaatgat gatgaaaaaa 
`)
}
