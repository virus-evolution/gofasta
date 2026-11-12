package genbank

import (
	"bytes"
	"fmt"
	"os"
	"runtime/debug"
	"strings"
	"testing"
)

// FuzzReadGenBank: parsing arbitrary text, and asking every feature for its positions and strand.
func FuzzReadGenBank(f *testing.F) {
	if b, err := os.ReadFile("../../resources/MN908947.gb"); err == nil {
		f.Add(b[:min(len(b), 6000)])
	}
	f.Add([]byte("LOCUS       X 10 bp\nFEATURES             Location/Qualifiers\n     CDS             1..9\n                     /gene=\"a\"\nORIGIN\n        1 atgaaataa\n//\n"))
	f.Fuzz(func(t *testing.T, data []byte) {
		collect(data, func() {
			gb, err := ReadGenBank(bytes.NewReader(data))
			if err != nil {
				return
			}
			for _, feat := range gb.FEATURES {
				_, _ = feat.Location.GetPositions()
				_, _ = feat.Location.IsReverse()
			}
		})
	})
}

// FuzzLocation: location strings on their own.
func FuzzLocation(f *testing.F) {
	for _, s := range []string{
		"1..9", "complement(1..9)", "join(1..3,5..9)", "complement(join(1..3,5..9))",
		"order(1..3,5..9)", "<1..9",
	} {
		f.Add(s)
	}
	f.Fuzz(func(t *testing.T, s string) {
		collect([]byte(s), func() {
			l := Location{Representation: s}
			_, _ = l.GetPositions()
			_, _ = l.IsReverse()
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
