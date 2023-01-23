package genbank

import (
	"bytes"
	"reflect"
	"testing"
)

func TestHasAttribute(t *testing.T) {
	f := GenbankFeature{
		Feature:  "CDS",
		Location: Location{Representation: "1..6"},
		Info: map[string]string{
			"translation": "MM",
		},
	}

	if !f.HasAttribute("translation") {
		t.Errorf("Problem in TestHasAttribute()")
	}

	if f.HasAttribute("gene") {
		t.Errorf("Problem in TestHasAttribute()")
	}
}

func TestParseGenbankFEATURES(t *testing.T) {
	gf := genbankField{
		header: "FEATURES             Location/Qualifiers",
		lines: []string{"      source          1..8959",
			"                     /organism=\"Homo sapiens\"",
			"                     /db_xref=\"taxon:9606\"",
			"                     /mol_type=\"genomic DNA\"",
			"      gene            212..8668",
			"                      /gene=\"NF1\"",
			"      CDS             212..8668",
			"                      /gene=\"NF1\"",
			"                      /note=\"putative\"",
			"                      /codon_start=1",
			"                      /product=\"GAP-related protein\"",
			"                      /protein_id=\"AAA59924.1\"",
			"                      /translation=\"MAAHRPVEWVQAVVSRFDEQLPIKTGQQNTHTKVSTE",
			"                      MAAHRPVEWVQAVVSRFDEQLPIKTGQQNTHTKVSTE\"",
		},
	}

	feats := parseGenbankFEATURES(gf)
	if !reflect.DeepEqual(feats, []GenbankFeature{
		GenbankFeature{
			Feature:  "source",
			Location: Location{Representation: "1..8959"},
			Info: map[string]string{
				"organism": "Homo sapiens",
				"db_xref":  "taxon:9606",
				"mol_type": "genomic DNA",
			},
		},
		GenbankFeature{
			Feature:  "gene",
			Location: Location{Representation: "212..8668"},
			Info: map[string]string{
				"gene": "NF1",
			},
		},
		GenbankFeature{
			Feature:  "CDS",
			Location: Location{Representation: "212..8668"},
			Info: map[string]string{
				"gene":        "NF1",
				"note":        "putative",
				"codon_start": "1",
				"product":     "GAP-related protein",
				"protein_id":  "AAA59924.1",
				"translation": "MAAHRPVEWVQAVVSRFDEQLPIKTGQQNTHTKVSTEMAAHRPVEWVQAVVSRFDEQLPIKTGQQNTHTKVSTE",
			},
		},
	}) {
		t.Errorf("Problem in TestParseGenbankFEATURES()")
	}
}

func TestParseGenbankORIGIN(t *testing.T) {
	gf := genbankField{
		header: "ORIGIN                                ",
		lines: []string{
			"        1 cgttatttaa ggtgttacat agttctatgg aaatagggtc tatacctttc gccttacaat",
			"       61 gtaatttctt",
		},
	}

	origin := parseGenbankORIGIN(gf)

	if !reflect.DeepEqual(origin, []byte("cgttatttaaggtgttacatagttctatggaaatagggtctatacctttcgccttacaatgtaatttctt")) {
		t.Errorf("Problem in TestParseGenbankORIGIN()")
	}
}

func TestReadGenbank(t *testing.T) {
	data := []byte(`FEATURES             Location/Qualifiers
      source          1..8959
                     /organism="Homo sapiens"
                     /db_xref="taxon:9606"
                     /mol_type="genomic DNA"
      gene            212..8668
                      /gene="NF1"
      CDS             212..8668
                      /gene="NF1"
                      /note="putative"
                      /codon_start=1
                      /product="GAP-related protein"
                      /protein_id="AAA59924.1"
                      /translation="MAAHRPVEWVQAVVSRFDEQLPIKTGQQNTHTKVSTE
                      MAAHRPVEWVQAVVSRFDEQLPIKTGQQNTHTKVSTE"
ORIGIN                                
        1 cgttatttaa ggtgttacat agttctatgg aaatagggtc tatacctttc gccttacaat
       61 gtaatttctt
`)

	reader := bytes.NewReader(data)

	gb, err := ReadGenBank(reader)
	if err != nil {
		t.Error(err)
	}

	desiredResult := Genbank{}
	desiredResult.FEATURES = []GenbankFeature{
		GenbankFeature{
			Feature:  "source",
			Location: Location{Representation: "1..8959"},
			Info: map[string]string{
				"organism": "Homo sapiens",
				"db_xref":  "taxon:9606",
				"mol_type": "genomic DNA",
			},
		},
		GenbankFeature{
			Feature:  "gene",
			Location: Location{Representation: "212..8668"},
			Info: map[string]string{
				"gene": "NF1",
			},
		},
		GenbankFeature{
			Feature:  "CDS",
			Location: Location{Representation: "212..8668"},
			Info: map[string]string{
				"gene":        "NF1",
				"note":        "putative",
				"codon_start": "1",
				"product":     "GAP-related protein",
				"protein_id":  "AAA59924.1",
				"translation": "MAAHRPVEWVQAVVSRFDEQLPIKTGQQNTHTKVSTEMAAHRPVEWVQAVVSRFDEQLPIKTGQQNTHTKVSTE",
			},
		},
	}
	desiredResult.ORIGIN = []byte("cgttatttaaggtgttacatagttctatggaaatagggtctatacctttcgccttacaatgtaatttctt")

	if !reflect.DeepEqual(gb, desiredResult) {
		t.Errorf("Problem in TestReadGenbank()")
	}
}
