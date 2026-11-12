package variants

import (
	"bytes"
	"fmt"
	"os"
	"runtime/debug"
	"strings"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/genbank"
	"github.com/virus-evolution/gofasta/pkg/gff"
)

// FuzzRegionsFromGenbank: any annotation the parser accepts must turn into regions (or an error), not a panic.
func FuzzRegionsFromGenbank(f *testing.F) {
	if b, err := os.ReadFile("../../resources/MN908947.gb"); err == nil {
		f.Add(b[:min(len(b), 6000)], 29903)
	}
	f.Add([]byte("LOCUS       X 10 bp\nFEATURES             Location/Qualifiers\n     CDS             1..9\n                     /gene=\"a\"\n                     /codon_start=1\n                     /translation=\"MK\"\nORIGIN\n        1 atgaaataa\n//\n"), 9)
	f.Fuzz(func(t *testing.T, data []byte, refLen int) {
		if refLen > 10000000 || refLen < -10000000 {
			return
		}
		collect(data, func() {
			gb, err := genbank.ReadGenBank(bytes.NewReader(data))
			if err != nil {
				return
			}
			_, _, _ = RegionsFromGenbank(gb, refLen)
		})
	})
}

// FuzzRegionsFromGFF: likewise for GFF plus a reference sequence.
func FuzzRegionsFromGFF(f *testing.F) {
	f.Add([]byte("##gff-version 3\ns\t.\tCDS\t1\t9\t.\t+\t0\tName=a;ID=b\n"), "ATGAAATAA")
	f.Fuzz(func(t *testing.T, data []byte, ref string) {
		collect(data, func() {
			g, err := gff.ReadGFF(bytes.NewReader(data))
			if err != nil {
				return
			}
			_, _, _ = RegionsFromGFF(g, ref)
		})
	})
}

// collect runs fn on data. Normally a panic fails the fuzz test. With GOFASTA_FUZZ_COLLECT=DIR it is
// recovered and the first input seen for each distinct panic site is written to DIR instead, so one run
// can survey every distinct crash rather than stopping at the first.
func collect(data []byte, fn func()) {
	dir := os.Getenv("GOFASTA_FUZZ_COLLECT")
	if dir == "" {
		fn()
		return
	}
	defer func() {
		if r := recover(); r != nil {
			site := "?"
			for _, l := range strings.Split(string(debug.Stack()), "\n") {
				if i := strings.Index(l, "/gofasta/pkg/"); i >= 0 && strings.HasPrefix(l, "\t") && !strings.Contains(l, "_test.go") {
					site = strings.Fields(l[i+len("/gofasta/"):])[0]
					break
				}
			}
			name := strings.NewReplacer("/", "_", ":", "_", ".", "_").Replace(fmt.Sprint(r, " ", site))
			if len(name) > 120 {
				name = name[:120]
			}
			_ = os.MkdirAll(dir, 0o755)
			_ = os.WriteFile(dir+"/"+name, data, 0o644)
		}
	}()
	fn()
}
