package variants

import (
	"bytes"
	"fmt"
	"reflect"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/genbank"
	"github.com/virus-evolution/gofasta/pkg/gff"
)

// type Region struct {
// 	Whichtype   string // only "protein-coding" for now
// 	Name        string // name of feature, if it has one
// 	Start       int    // 1-based 5'-most position of region on the forward strand, inclusive
// 	Stop        int    // 1-based 3'-most position of region on the forward strand, inclusive
// 	Translation string // amino acid sequence of this region if it is CDS
// 	Strand      int    // values in the set {-1, +1} only (and "0" for a mixture?!)
// 	Positions   []int  // all the (1-based, unadjusted) positions in order, on the reverse strand if needs be
// }

// TO DO: maybe move this to the genbank package?
func TestGetRegionsGenbank(t *testing.T) {
	genbankReader := bytes.NewReader(genbankDataShort)
	gb, err := genbank.ReadGenBank(genbankReader)
	if err != nil {
		t.Error(err)
	}

	cdsregions, intregions, err := RegionsFromGenbank(gb, 23)
	if err != nil {
		t.Error(err)
	}

	var desiredCDSResult = []Region{
		{Whichtype: "protein-coding", Name: "gene1", Strand: 1, Start: 6, Stop: 17, Translation: "MMM*", Positions: []int{6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17}},
	}

	if !reflect.DeepEqual(cdsregions, desiredCDSResult) {
		t.Errorf("problem in TestGetRegionsGenbank")
		fmt.Println(cdsregions)
	}

	var desiredInterResult = []int{1, 2, 3, 4, 5, 18, 19, 20, 21, 22, 23}
	if !reflect.DeepEqual(intregions, desiredInterResult) {
		t.Errorf("problem in TestGetRegionsGenbank")
		fmt.Println(intregions)
	}
}

func TestGetRegionsGFF(t *testing.T) {
	gffReader := bytes.NewReader(gffDataShort)
	GFF, err := gff.ReadGFF(gffReader)
	if err != nil {
		t.Error(err)
	}

	cdsregions, intregions, err := RegionsFromGFF(GFF, GFF.FASTA["somefakething"].Seq)
	if err != nil {
		t.Error(err)
	}

	desiredCDSResult := []Region{
		{Whichtype: "protein-coding", Name: "gene1", Strand: 1, Start: 6, Stop: 17, Translation: "MMM*", Positions: []int{6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17}},
	}

	if !reflect.DeepEqual(cdsregions, desiredCDSResult) {
		t.Errorf("problem in TestGetRegionsGFF()")
		fmt.Println(cdsregions)
	}

	desiredInterResult := []int{1, 2, 3, 4, 5, 18, 19, 20, 21, 22, 23}
	if !reflect.DeepEqual(intregions, desiredInterResult) {
		t.Errorf("problem in TestGetRegionsGFF()")
		fmt.Println(intregions)
	}

	gffReader = bytes.NewReader(gffDataShortRev)
	GFF, err = gff.ReadGFF(gffReader)
	if err != nil {
		t.Error(err)
	}

	cdsregions, intregions, err = RegionsFromGFF(GFF, GFF.FASTA["somefakething"].Seq)
	if err != nil {
		t.Error(err)
	}

	desiredCDSResult = []Region{
		{Whichtype: "protein-coding", Name: "gene1", Strand: -1, Start: 7, Stop: 18, Translation: "MMM*", Positions: []int{18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7}},
	}

	if !reflect.DeepEqual(cdsregions, desiredCDSResult) {
		t.Errorf("problem in TestGetRegionsGFF(rev)")
		fmt.Println(cdsregions)
	}

	desiredInterResult = []int{1, 2, 3, 4, 5, 6, 19, 20, 21, 22, 23}
	if !reflect.DeepEqual(intregions, desiredInterResult) {
		t.Errorf("problem in TestGetRegionsGFF(rev)")
		fmt.Println(intregions)
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

	desiredResult := fasta.EncodedRecord{ID: "MN908947.3", Description: "MN908947.3", Seq: []byte{136, 24, 72, 136, 24, 72, 136, 24, 72}, Idx: 0}
	if !reflect.DeepEqual(ref, desiredResult) {
		t.Errorf("problem in TestfindReference")
	}

	msaReader = bytes.NewReader(msaData)
	ref, err = findReference(msaReader, "nottheref")
	if err != nil {
		t.Error(err)
	}

	desiredResult = fasta.EncodedRecord{ID: "nottheref", Description: "nottheref hey", Seq: []byte{136, 136, 136, 136, 136, 136, 136, 136, 136}, Idx: 1}
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
ACGTAATGATGATGTAG-AAAAAA
>seq1
ATGTATTGATGATGTAG-AAAATA
>seq2
ACGTA---ATGATGTAG-AAAAAA
>seq3
ACGTAATGATGATGTAGAAAAAAA
>seq4
ACGTAATGATGAG-TAG-TAAA-T
`)

	genbankReader := bytes.NewReader(genbankDataShort)
	gb, err := genbank.ReadGenBank(genbankReader)
	if err != nil {
		t.Error(err)
	}

	msaReader := bytes.NewReader(msaData)
	ref, err := findReference(msaReader, "reference")
	if err != nil {
		t.Error(err)
	}
	refToMSA, MSAToRef := GetMSAOffsets(ref.Seq)
	refLenDegapped := len(ref.Decode().Degap().Seq)

	cdsregions, intregions, err := RegionsFromGenbank(gb, refLenDegapped)
	if err != nil {
		t.Error(err)
	}

	msaReader = bytes.NewReader(msaData)

	queries, err := fasta.LoadEncodeAlignment(msaReader, false, false, false)
	if err != nil {
		t.Error(err)
	}

	mutations, err := GetVariantsPair(ref.Seq, queries[1].Seq, "reference", queries[1].ID, 1, cdsregions, intregions, refToMSA, MSAToRef)
	if err != nil {
		t.Error(err)
	}

	desiredResult := AnnoStructs{Queryname: "seq1", Vs: []Variant{
		{RefAl: "C", QueAl: "T", Position: 2, Changetype: "nuc"},
		{RefAl: "M", QueAl: "L", Position: 6, Residue: 1, Changetype: "aa", Feature: "gene1", SNPs: "nuc:A6T"},
		{RefAl: "A", QueAl: "T", Position: 22, Changetype: "nuc"},
	}, Idx: 1}
	if !reflect.DeepEqual(mutations, desiredResult) {
		t.Errorf("problem in TestGetVariantsPair (seq1)")
		fmt.Println(mutations)
	}

	mutations, err = GetVariantsPair(ref.Seq, queries[2].Seq, "reference", queries[2].ID, 2, cdsregions, intregions, refToMSA, MSAToRef)
	if err != nil {
		t.Error(err)
	}

	desiredResult = AnnoStructs{Queryname: "seq2", Vs: []Variant{
		{Position: 6, Length: 3, Changetype: "del"},
	}, Idx: 2}
	if !reflect.DeepEqual(mutations, desiredResult) {
		t.Errorf("problem in TestGetVariantsPair (seq2)")
	}

	mutations, err = GetVariantsPair(ref.Seq, queries[3].Seq, "reference", queries[3].ID, 3, cdsregions, intregions, refToMSA, MSAToRef)
	if err != nil {
		t.Error(err)
	}

	desiredResult = AnnoStructs{Queryname: "seq3", Vs: []Variant{
		{Position: 17, Length: 1, Changetype: "ins"},
	}, Idx: 3}
	if !reflect.DeepEqual(mutations, desiredResult) {
		t.Errorf("problem in TestGetVariantsPair (seq3)")
	}

	mutations, err = GetVariantsPair(ref.Seq, queries[4].Seq, "reference", queries[4].ID, 4, cdsregions, intregions, refToMSA, MSAToRef)
	if err != nil {
		t.Error(err)
	}

	desiredResult = AnnoStructs{Queryname: "seq4", Vs: []Variant{
		{Position: 13, RefAl: "T", QueAl: "G", Changetype: "nuc"},
		{Position: 14, Length: 1, Changetype: "del"},
		{Position: 18, RefAl: "A", QueAl: "T", Changetype: "nuc"},
		{Position: 22, Length: 1, Changetype: "del"},
		{Position: 23, RefAl: "A", QueAl: "T", Changetype: "nuc"},
	}, Idx: 4}
	if !reflect.DeepEqual(mutations, desiredResult) {
		t.Errorf("problem in TestGetVariantsPair (seq4)")
	}
}

func TestFormatVariant(t *testing.T) {
	msaData := []byte(`>reference
ACGTAATGATGATGTAG-AAAAAA
>seq1
ATGTATTGATGATGTAG-AAAATA
>seq2
ACGTA---ATGATGTAG-AAAAAA
>seq3
ACGTAATGATGATGTAGAAAAAAA
>seq4
ACGTAATGATGAG-TAG-TAAA-T
`)

	genbankReader := bytes.NewReader(genbankDataShort)
	gb, err := genbank.ReadGenBank(genbankReader)
	if err != nil {
		t.Error(err)
	}

	msaReader := bytes.NewReader(msaData)
	ref, err := findReference(msaReader, "reference")
	if err != nil {
		t.Error(err)
	}
	refToMSA, MSAToRef := GetMSAOffsets(ref.Seq)
	refLenDegapped := len(ref.Decode().Degap().Seq)

	cdsregions, intregions, err := RegionsFromGenbank(gb, refLenDegapped)
	if err != nil {
		t.Error(err)
	}

	msaReader = bytes.NewReader(msaData)

	queries, err := fasta.LoadEncodeAlignment(msaReader, false, false, false)
	if err != nil {
		t.Error(err)
	}

	mutations, err := GetVariantsPair(ref.Seq, queries[1].Seq, "reference", queries[1].ID, 1, cdsregions, intregions, refToMSA, MSAToRef)
	if err != nil {
		t.Error(err)
	}

	/*
		>seq1
		nuc:C2T, gene1:M1L, nuc:A22T
		>seq2
		del:6:3
		>seq3
		ins:17:1
		>seq4
		nuc:T13G, del:14:1, nuc:A18T, del:22:1, nuc:A23T
	*/
	desiredResult := []string{"nuc:C2T", "aa:gene1:M1L", "nuc:A22T"}
	for i, mutation := range mutations.Vs {
		mut, err := FormatVariant(mutation, false)
		if err != nil {
			t.Error(err)
		}
		if mut != desiredResult[i] {
			t.Errorf("problem in TestFormatVariant 1")
		}
	}

	desiredResult = []string{"nuc:C2T", "aa:gene1:M1L(nuc:A6T)", "nuc:A22T"}
	for i, mutation := range mutations.Vs {
		mut, err := FormatVariant(mutation, true)
		if err != nil {
			t.Error(err)
		}
		if mut != desiredResult[i] {
			t.Errorf("problem in TestFormatVariant 2")
		}
	}

	mutations, err = GetVariantsPair(ref.Seq, queries[4].Seq, "reference", queries[4].ID, 4, cdsregions, intregions, refToMSA, MSAToRef)
	if err != nil {
		t.Error(err)
	}

	desiredResult = []string{"nuc:T13G", "del:14:1", "nuc:A18T", "del:22:1", "nuc:A23T"}
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
var gffDataShort []byte
var gffDataShortRev []byte

func init() {
	genbankDataShort = []byte(`LOCUS       TEST               23 bp ss-RNA     linear   VRL 21-MAR-1987
FEATURES             Location/Qualifiers
		source          1..23
						/organism="Not a real organism"
		5'UTR           1..5
		gene            6..17
						/gene="gene1"
		CDS             6..17
						/gene="gene1"
						/codon_start=1
						/translation="MMM"
		3'UTR           18..23
ORIGIN      
		1 acgtaatgat gatgtagaaa aaa 
`)

	gffDataShort = []byte(`##gff-version 3
##sequence-region somefakething 1 23
somefakething	RefSeq	region	1	23	.	+	.	ID=somefakething:1..23
somefakething	RefSeq	five_prime_UTR	1	5	.	+	.	ID=somefakething:1..5
somefakething	RefSeq	gene	6	17	.	+	.	ID=gene1
somefakething	RefSeq	CDS	6	17	.	+	0	ID=CDS-gene1;Parent=gene1;Name=gene1
somefakething	RefSeq	three_prime_UTR	18	23	.	+	.	ID=somefakething:18..23
##FASTA
>somefakething
acgtaatgatgatgtagaaaaaa
`)

	gffDataShortRev = []byte(`##gff-version 3
##sequence-region somefakething 1 23
somefakething	RefSeq	region	1	23	.	+	.	ID=somefakething:1..23
somefakething	RefSeq	five_prime_UTR	1	6	.	+	.	ID=somefakething:1..5
somefakething	RefSeq	gene	7	18	.	-	.	ID=gene1
somefakething	RefSeq	CDS	7	18	.	-	0	ID=CDS-gene1;Parent=gene1;Name=gene1
somefakething	RefSeq	three_prime_UTR	19	23	.	+	.	ID=somefakething:18..23
##FASTA
>somefakething
ttttttctacatcatcattacgt
`)
}
