package gff

import (
	"bytes"
	"fmt"
	"os"
	"runtime/debug"
	"strings"
	"testing"
)

// FuzzReadGFF: parsing arbitrary text.
func FuzzReadGFF(f *testing.F) {
	if b, err := os.ReadFile("../../resources/sarscov2-reduced.gff"); err == nil {
		f.Add(b[:min(len(b), 4000)])
	}
	f.Add([]byte("##gff-version 3\n##sequence-region s 1 10\ns\t.\tCDS\t1\t9\t.\t+\t0\tName=a;ID=b\n##FASTA\n>s\nACGTACGTAC\n"))
	f.Add([]byte("##gff-version 3\n\ns\t.\tCDS\t1\t9\t.\t?\t.\tName=a;\n"))
	f.Add([]byte("##sequence-region\n##gff-version\n"))
	f.Fuzz(func(t *testing.T, data []byte) {
		collect(data, func() { _, _ = ReadGFF(bytes.NewReader(data)) })
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
