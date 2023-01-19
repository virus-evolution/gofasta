package gff

import (
	"bytes"
	"fmt"
	"reflect"
	"testing"
)

func TestVersionStringFromHeader(t *testing.T) {
	gff := GFF{
		HeaderLines: []string{
			"gff-version 3",
			"sequence-region NC_045512.2 1 29903",
		},
	}

	err := gff.versionStringFromHeader()
	if err != nil {
		t.Error(err)
	}

	if gff.GFF_version != "3" {
		t.Error("Problem in TestVersionStringFromHeader() (1)")
		fmt.Println(gff.GFF_version)
	}

	data := []byte(`##gff-version 3.1.26
# gff-spec-version 1.26 (https://github.com/The-Sequence-Ontology/Specifications/blob/fe73505276dd324bf6a55773f3413fe2bed47af4/gff3.md)
##sequence-region NC_045512.2 1 29903
NC_045512.2	RefSeq	CDS	266	13468	.	+	0	ID=CDS-pp1ab
`)

	reader := bytes.NewReader(data)

	gff, err = ReadGFF(reader)
	if err != nil {
		t.Error(err)
	}

	if gff.GFF_version != "3.1.26" {
		t.Error("Problem in TestVersionStringFromHeader() (2)")
		fmt.Println(gff.GFF_version)
	}
}

func TestSetSequenceRegionsFromHeader(t *testing.T) {
	gff := GFF{
		HeaderLines: []string{
			"gff-version 3",
			"sequence-region NC_045512.2 1 29903",
		},
		SequenceRegions: make(map[string]SequenceRegion),
	}

	err := gff.setSequenceRegionsFromHeader()
	if err != nil {
		t.Error(err)
	}

	keys := make([]string, 0)
	for k, _ := range gff.SequenceRegions {
		keys = append(keys, k)
	}

	if !reflect.DeepEqual(keys, []string{"NC_045512.2"}) {
		t.Errorf("Problem in TestSetSequenceRegionsFromHeader()")
	} else {
		if !reflect.DeepEqual(gff.SequenceRegions["NC_045512.2"], SequenceRegion{Seqid: "NC_045512.2", Start: 1, End: 29903}) {
			t.Errorf("Problem in TestSetSequenceRegionsFromHeader()")
		}
	}

	data := []byte(`##gff-version 3.1.26
# gff-spec-version 1.26 (https://github.com/The-Sequence-Ontology/Specifications/blob/fe73505276dd324bf6a55773f3413fe2bed47af4/gff3.md)
##sequence-region NC_045512.2 1 29903
NC_045512.2	RefSeq	CDS	266	13468	.	+	0	ID=CDS-pp1ab
`)

	reader := bytes.NewReader(data)

	gff, err = ReadGFF(reader)
	if err != nil {
		t.Error(err)
	}

	keys = make([]string, 0)
	for k, _ := range gff.SequenceRegions {
		keys = append(keys, k)
	}

	if !reflect.DeepEqual(keys, []string{"NC_045512.2"}) {
		t.Errorf("Problem in TestSetSequenceRegionsFromHeader()")
	} else {
		if !reflect.DeepEqual(gff.SequenceRegions["NC_045512.2"], SequenceRegion{Seqid: "NC_045512.2", Start: 1, End: 29903}) {
			t.Errorf("Problem in TestSetSequenceRegionsFromHeader()")
		}
	}
}

func TestFeatureFromLine(t *testing.T) {
	lines := []string{
		"NC_045512.2	RefSeq	CDS	266	13468	.	+	0	ID=CDS-pp1ab",
		"NC_045512.2	RefSeq	CDS	13468	21555	.	+	0	ID=CDS-pp1ab",
		"NC_045512.2	RefSeq	mature_protein_region_of_CDS	266	805	.	+	.	ID=nsp1;Parent=CDS-pp1ab;Name=nsp1",
	}

	feats := make([]Feature, 0)
	for _, line := range lines {
		feat, err := featureFromLine(line)
		if err != nil {
			t.Error(err)
		}
		feats = append(feats, feat)
	}

	if !reflect.DeepEqual(feats[0], Feature{
		Seqid:  "NC_045512.2",
		Source: "RefSeq",
		Type:   "CDS",
		Start:  266,
		End:    13468,
		Score:  ".",
		Strand: "+",
		Phase:  0,
		Attributes: map[string][]string{
			"ID": []string{"CDS-pp1ab"},
		},
	}) {
		t.Errorf("Problem in TestFeatureFromLine()")
	}

	if !reflect.DeepEqual(feats[1], Feature{
		Seqid:  "NC_045512.2",
		Source: "RefSeq",
		Type:   "CDS",
		Start:  13468,
		End:    21555,
		Score:  ".",
		Strand: "+",
		Phase:  0,
		Attributes: map[string][]string{
			"ID": []string{"CDS-pp1ab"},
		},
	}) {
		t.Errorf("Problem in TestFeatureFromLine()")
	}

	if !reflect.DeepEqual(feats[2], Feature{
		Seqid:  "NC_045512.2",
		Source: "RefSeq",
		Type:   "mature_protein_region_of_CDS",
		Start:  266,
		End:    805,
		Score:  ".",
		Strand: "+",
		Phase:  0,
		Attributes: map[string][]string{
			"ID":     []string{"nsp1"},
			"Parent": []string{"CDS-pp1ab"},
			"Name":   []string{"nsp1"},
		},
	}) {
		t.Errorf("Problem in TestFeatureFromLine()")
	}
}

func TestPopulateIDMap(t *testing.T) {

	data := []byte(`##gff-version 3
# gff-spec-version 1.26 (https://github.com/The-Sequence-Ontology/Specifications/blob/fe73505276dd324bf6a55773f3413fe2bed47af4/gff3.md)
##sequence-region NC_045512.2 1 29903
NC_045512.2	RefSeq	CDS	266	13468	.	+	0	ID=CDS-pp1ab
NC_045512.2	RefSeq	CDS	13468	21555	.	+	0	ID=CDS-pp1ab
NC_045512.2	RefSeq	mature_protein_region_of_CDS	266	805	.	+	.	ID=nsp1;Parent=CDS-pp1ab;Name=nsp1
`)

	reader := bytes.NewReader(data)

	gff, err := ReadGFF(reader)
	if err != nil {
		t.Error(err)
	}

	if !reflect.DeepEqual(gff.IDmap, map[string][]int{
		"CDS-pp1ab": []int{0, 1},
		"nsp1":      []int{2},
	}) {
		t.Errorf("Problem in TestPopulateIDMap()")
		fmt.Println(gff.IDmap)
	}
}

func TestFeaturesOfType(t *testing.T) {
	data := []byte(`##gff-version 3
# gff-spec-version 1.26 (https://github.com/The-Sequence-Ontology/Specifications/blob/fe73505276dd324bf6a55773f3413fe2bed47af4/gff3.md)
##sequence-region NC_045512.2 1 29903
NC_045512.2	RefSeq	CDS	266	13468	.	+	0	ID=CDS-pp1ab
NC_045512.2	RefSeq	CDS	13468	21555	.	+	0	ID=CDS-pp1ab
NC_045512.2	RefSeq	mature_protein_region_of_CDS	266	805	.	+	.	ID=nsp1;Parent=CDS-pp1ab;Name=nsp1
`)

	reader := bytes.NewReader(data)

	gff, err := ReadGFF(reader)
	if err != nil {
		t.Error(err)
	}

	if !reflect.DeepEqual(gff.featuresOfType("CDS"), []Feature{
		Feature{
			Seqid:  "NC_045512.2",
			Source: "RefSeq",
			Type:   "CDS",
			Start:  266,
			End:    13468,
			Score:  ".",
			Strand: "+",
			Phase:  0,
			Attributes: map[string][]string{
				"ID": []string{"CDS-pp1ab"},
			},
		},
		Feature{
			Seqid:  "NC_045512.2",
			Source: "RefSeq",
			Type:   "CDS",
			Start:  13468,
			End:    21555,
			Score:  ".",
			Strand: "+",
			Phase:  0,
			Attributes: map[string][]string{
				"ID": []string{"CDS-pp1ab"},
			},
		},
	}) {
		t.Errorf("Problem in TestFeaturesOfType()")
	}

	if !reflect.DeepEqual(gff.featuresOfType("mature_protein_region_of_CDS"), []Feature{
		Feature{
			Seqid:  "NC_045512.2",
			Source: "RefSeq",
			Type:   "mature_protein_region_of_CDS",
			Start:  266,
			End:    805,
			Score:  ".",
			Strand: "+",
			Phase:  0,
			Attributes: map[string][]string{
				"ID":     []string{"nsp1"},
				"Parent": []string{"CDS-pp1ab"},
				"Name":   []string{"nsp1"},
			},
		},
	}) {
		t.Errorf("Problem in TestFeaturesOfType()")
	}
}

func TestFeaturesFromID(t *testing.T) {
	data := []byte(`##gff-version 3
# gff-spec-version 1.26 (https://github.com/The-Sequence-Ontology/Specifications/blob/fe73505276dd324bf6a55773f3413fe2bed47af4/gff3.md)
##sequence-region NC_045512.2 1 29903
NC_045512.2	RefSeq	CDS	266	13468	.	+	0	ID=CDS-pp1ab
NC_045512.2	RefSeq	CDS	13468	21555	.	+	0	ID=CDS-pp1ab
NC_045512.2	RefSeq	mature_protein_region_of_CDS	266	805	.	+	.	ID=nsp1;Parent=CDS-pp1ab;Name=nsp1
`)

	reader := bytes.NewReader(data)

	gff, err := ReadGFF(reader)
	if err != nil {
		t.Error(err)
	}

	if !reflect.DeepEqual(gff.featuresFromID("CDS-pp1ab"), []Feature{
		Feature{
			Seqid:  "NC_045512.2",
			Source: "RefSeq",
			Type:   "CDS",
			Start:  266,
			End:    13468,
			Score:  ".",
			Strand: "+",
			Phase:  0,
			Attributes: map[string][]string{
				"ID": []string{"CDS-pp1ab"},
			},
		},
		Feature{
			Seqid:  "NC_045512.2",
			Source: "RefSeq",
			Type:   "CDS",
			Start:  13468,
			End:    21555,
			Score:  ".",
			Strand: "+",
			Phase:  0,
			Attributes: map[string][]string{
				"ID": []string{"CDS-pp1ab"},
			},
		},
	}) {
		t.Errorf("Problem in TestFeaturesFromID()")
	}
}

func TestHasAttribute(t *testing.T) {
	F := Feature{
		Seqid:  "NC_045512.2",
		Source: "RefSeq",
		Type:   "CDS",
		Start:  13468,
		End:    21555,
		Score:  ".",
		Strand: "+",
		Phase:  0,
		Attributes: map[string][]string{
			"ID": []string{"CDS-pp1ab"},
		},
	}

	if !F.HasAttribute("ID") {
		t.Errorf("Problem in TestHasAttribute()")
	}
	if F.HasAttribute("Name") {
		t.Errorf("Problem in TestHasAttribute()")
	}

	F = Feature{
		Seqid:  "NC_045512.2",
		Source: "RefSeq",
		Type:   "mature_protein_region_of_CDS",
		Start:  266,
		End:    805,
		Score:  ".",
		Strand: "+",
		Phase:  0,
		Attributes: map[string][]string{
			"ID":     []string{"nsp1"},
			"Parent": []string{"CDS-pp1ab"},
			"Name":   []string{"nsp1"},
		},
	}

	if !F.HasAttribute("ID") {
		t.Errorf("Problem in TestHasAttribute()")
	}
	if !F.HasAttribute("Parent") {
		t.Errorf("Problem in TestHasAttribute()")
	}
	if !F.HasAttribute("Name") {
		t.Errorf("Problem in TestHasAttribute()")
	}
	if F.HasAttribute("Exon") {
		t.Errorf("Problem in TestHasAttribute()")
	}
}

func TestIsEscapedCorrectly(t *testing.T) {
	s := ">Seq1"
	err := isEscapedCorrectly(s)
	if !reflect.DeepEqual(err, errorBuilder(errGFFParsingSeqID, s)) {
		t.Errorf("Problem in TestIsEscapedCorrectly()")
	}

	s = "Seq1 i think"
	err = isEscapedCorrectly(s)
	if !reflect.DeepEqual(err, errorBuilder(errGFFParsingSeqID, s)) {
		t.Errorf("Problem in TestIsEscapedCorrectly()")
	}

}

func TestSeqidFromField(t *testing.T) {

}
