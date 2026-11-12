package fasta

import (
	"bytes"
	"fmt"
	"os"
	"runtime/debug"
	"strings"
	"testing"
)

// FuzzReader: reading arbitrary bytes must neither panic nor loop forever, and must end.
func FuzzReader(f *testing.F) {
	for _, s := range []string{
		">a\nACGT\n>b\nAC-T\n", ">a\nAC\n>b", "", "\n", "AC\n>a\n",
		">a\r\nAC\r\n", ">a\n\nAC\n", ">a\nAC", ">a b c\nAC\nGT\n",
	} {
		f.Add([]byte(s))
	}
	f.Fuzz(func(t *testing.T, data []byte) {
		collect(data, func() {
			r := NewReader(bytes.NewReader(data))
			// each successful Read consumes at least one byte, so this bounds the loop
			for i := 0; i <= len(data)+1; i++ {
				_, err := r.Read()
				if err != nil {
					return
				}
			}
			t.Fatalf("Reader did not reach the end of %d bytes", len(data))
		})
	})
}

// FuzzStreamAlignment: the alignment reader on arbitrary bytes.
func FuzzGetAlignmentDims(f *testing.F) {
	for _, s := range []string{">a\nACGT\n>b\nAC-T\n", ">a\nAC\n>b\nACG\n", ">a\n", ""} {
		f.Add([]byte(s))
	}
	f.Fuzz(func(t *testing.T, data []byte) {
		collect(data, func() { _, _, _ = getAlignmentDims(bytes.NewReader(data)) })
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
