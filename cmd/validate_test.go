package cmd

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// resetFlags sets every flag of c back to its default
func resetFlags(t *testing.T, c *cobra.Command) {
	c.InheritedFlags() // merges the parent's persistent flags into c.Flags()
	c.Flags().VisitAll(func(f *pflag.Flag) {
		if err := f.Value.Set(f.DefValue); err != nil {
			t.Fatalf("resetting --%s: %v", f.Name, err)
		}
		f.Changed = false
	})
}

// preRun resets every flag of c to its default, sets the given ones, and runs c's PreRunE. The flags are
// reset again when the test ends, so later tests see the defaults.
func preRun(t *testing.T, c *cobra.Command, set map[string]string) error {
	t.Helper()
	resetFlags(t, c)
	t.Cleanup(func() { resetFlags(t, c) })
	for k, v := range set {
		if err := c.Flags().Set(k, v); err != nil {
			t.Fatalf("setting --%s: %v", k, err)
		}
	}
	return c.PreRunE(c, nil)
}

func TestFlagValidation(t *testing.T) {
	dir := t.TempDir()
	f := filepath.Join(dir, "in.fa")
	if err := os.WriteFile(f, []byte(">a\nA\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	other := filepath.Join(dir, "out.csv")

	type tc struct {
		c   *cobra.Command
		set map[string]string
	}

	bad := []tc{
		{closestCmd, map[string]string{"threads": "-1"}},
		{closestCmd, map[string]string{"number": "-3"}},
		{closestCmd, map[string]string{"query": f, "outfile": f}},
		{closestCmd, map[string]string{"max-dist": "-1"}},
		{closestCmd, map[string]string{"max-dist": "-0.5"}},
		{closestCmd, map[string]string{"max-dist": "nan"}},
		{closestCmd, map[string]string{"max-dist": "x"}},
		{snpCmd, map[string]string{"threshold": "-0.5"}},
		{snpCmd, map[string]string{"reference": f, "outfile": f}},
		{variantsCmd, map[string]string{"threads": "0"}},
		{variantsCmd, map[string]string{"threshold": "-0.2"}},
		{variantsCmd, map[string]string{"msa": f, "outfile": f}},
		{variantsCmd, map[string]string{"start": "9", "end": "5"}},
		{samVariantsCmd, map[string]string{"start": "9", "end": "5"}},
		{samVariantsCmd, map[string]string{"threads": "0"}},
		{samVariantsCmd, map[string]string{"threshold": "-0.5"}},
		{samVariantsCmd, map[string]string{"reference": f, "outfile": f}},
		{toMultiAlignCmd, map[string]string{"samfile": f, "fasta-out": f}},
		{toPairAlignCmd, map[string]string{"threads": "-1"}},
		{toprankingCmd, map[string]string{"size-up": "-2"}},
		{toprankingCmd, map[string]string{"dist-push": "-1"}},
		{toprankingCmd, map[string]string{"threshold-pair": "-0.1"}},
		{toprankingCmd, map[string]string{"ignore": f, "outfile": f}},
		{updownListCmd, map[string]string{"query": f, "outfile": f}},
	}
	for _, b := range bad {
		if err := preRun(t, b.c, b.set); err == nil {
			t.Errorf("%s %v: expected an error", b.c.CommandPath(), b.set)
		}
	}

	good := []tc{
		{closestCmd, map[string]string{}},
		{closestCmd, map[string]string{"threads": "0"}},
		{closestCmd, map[string]string{"number": "5", "table": "true"}},
		{closestCmd, map[string]string{"table": "true"}}, // ignored without -n/-d
		{closestCmd, map[string]string{"max-dist": "0"}},
		{closestCmd, map[string]string{"max-dist": "inf"}},
		{closestCmd, map[string]string{"query": f, "outfile": other}},
		{snpCmd, map[string]string{"threshold": "0.5"}},
		{variantsCmd, map[string]string{"start": "5", "end": "5", "threshold": "0.1"}},
		{variantsCmd, map[string]string{"start": "5"}},
		{variantsCmd, map[string]string{"end": "5"}},
		{samVariantsCmd, map[string]string{"start": "5"}},
		{variantsCmd, map[string]string{"reference": f, "outfile": f}}, // --reference is a record ID here
		{samVariantsCmd, map[string]string{"threads": "4", "reference": f, "outfile": other}},
		{toMultiAlignCmd, map[string]string{}},
		{toprankingCmd, map[string]string{"size-up": "2", "threshold-pair": "0.2"}},
		{updownListCmd, map[string]string{}},
	}
	for _, g := range good {
		if err := preRun(t, g.c, g.set); err != nil {
			t.Errorf("%s %v: unexpected error %v", g.c.CommandPath(), g.set, err)
		}
	}
}

func TestAnnotationFormat(t *testing.T) {
	for path, want := range map[string]string{
		"a.gb": "gb", "a.gbk": "gb", "a.GBFF": "gb", "dir.x/a.gff": "gff", "a.gff3": "gff", "A.GFF3": "gff",
	} {
		if got, err := annotationFormat(path); err != nil || got != want {
			t.Errorf("annotationFormat(%q) = %q, %v; want %q", path, got, err, want)
		}
	}
	for _, path := range []string{"a.gtf", "a", "a.gff.gz", ""} {
		if _, err := annotationFormat(path); err == nil {
			t.Errorf("annotationFormat(%q): want an error", path)
		}
	}
}
