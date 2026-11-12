package variants

import (
	"strings"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/genbank"
	"github.com/virus-evolution/gofasta/pkg/gff"
)

// Malformed annotations give an error, never a panic.
func TestMalformedGenbankCDSGivesError(t *testing.T) {
	cds := func(loc, quals string) string {
		return "LOCUS X\nFEATURES             Location/Qualifiers\n     CDS             " + loc + "\n" + quals + "ORIGIN\n        1 atgaaataa\n//\n"
	}
	good := "                     /gene=\"a\"\n                     /codon_start=1\n                     /translation=\"MK\"\n"
	for name, in := range map[string]string{
		"past the reference":    cds("1..30", good),
		"codon_start too big":   cds("1..9", strings.Replace(good, "codon_start=1", "codon_start=100", 1)),
		"codon_start zero":      cds("1..9", strings.Replace(good, "codon_start=1", "codon_start=0", 1)),
		"codon_start not int":   cds("1..9", strings.Replace(good, "codon_start=1", "codon_start=x", 1)),
		"no translation":        cds("1..9", "                     /gene=\"a\"\n                     /codon_start=1\n"),
		"translation too short": cds("1..9", strings.Replace(good, "MK", "", 1)),
	} {
		gb, err := genbank.ReadGenBank(strings.NewReader(in))
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		if _, _, err := RegionsFromGenbank(gb, 9); err == nil {
			t.Errorf("%s: expected an error", name)
		}
	}
	gb, _ := genbank.ReadGenBank(strings.NewReader(cds("1..9", good)))
	if _, _, err := RegionsFromGenbank(gb, 9); err != nil {
		t.Errorf("well-formed CDS: %v", err)
	}
	if _, _, err := RegionsFromGenbank(gb, -1); err == nil {
		t.Errorf("negative reference length: expected an error")
	}
}

func TestMalformedGFFCDSGivesError(t *testing.T) {
	for name, line := range map[string]string{
		"past the reference": "s\t.\tCDS\t1\t99\t.\t+\t0\tName=a",
		"end before start":   "s\t.\tCDS\t5\t1\t.\t+\t0\tName=a",
		"start zero":         "s\t.\tCDS\t0\t9\t.\t+\t0\tName=a",
		"reverse, past":      "s\t.\tCDS\t3\t30\t.\t-\t0\tName=a",
		"unknown strand":     "s\t.\tCDS\t1\t9\t.\t?\t0\tName=a",
		"unstranded":         "s\t.\tCDS\t1\t9\t.\t.\t0\tName=a",
		"? then +":           "s\t.\tCDS\t1\t3\t.\t?\t0\tID=c;Name=a\ns\t.\tCDS\t4\t9\t.\t+\t0\tID=c;Name=a",
		"+ then ?":           "s\t.\tCDS\t1\t3\t.\t+\t0\tID=c;Name=a\ns\t.\tCDS\t4\t9\t.\t?\t0\tID=c;Name=a",
	} {
		g, err := gff.ReadGFF(strings.NewReader("##gff-version 3\n" + line + "\n"))
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		if _, _, err := RegionsFromGFF(g, "ATGAAATAA"); err == nil {
			t.Errorf("%s: expected an error", name)
		}
	}
}

// A CDS is reported, under its Name, only if it has a Name attribute; gene= is not used.
func TestGFFCDSNaming(t *testing.T) {
	for attrs, want := range map[string]string{
		"ID=cds-YP_1.1;Name=YP_1.1;gene=OPG001": "YP_1.1",
		"Name=S":                                "S",
		"ID=x;gene=N":                           "",
		"ID=x":                                  "",
	} {
		g, err := gff.ReadGFF(strings.NewReader("##gff-version 3\ns\t.\tCDS\t1\t9\t.\t+\t0\t" + attrs + "\n"))
		if err != nil {
			t.Fatal(err)
		}
		cds, _, err := RegionsFromGFF(g, "atgaaataa")
		if err != nil {
			t.Fatalf("%s: %v", attrs, err)
		}
		got := ""
		if len(cds) == 1 {
			got = cds[0].Name
		}
		if len(cds) > 1 || got != want {
			t.Errorf("%s: got %d regions named %q, want %q", attrs, len(cds), got, want)
		}
	}
}
