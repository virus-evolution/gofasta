package sam

import (
	"bytes"
	"io"
	"os"
	"strconv"
	"strings"
	"testing"
	"time"
)

const fuzzRef = "ATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAATAA"

const fuzzGenbank = `LOCUS       X                 45 bp    RNA     linear   VRL 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             1..45
                     /gene="a"
                     /codon_start=1
                     /translation="MKPGFKPGFKPGFK"
ORIGIN
        1 atgaaacccg ggtttaaacc cgggtttaaa cccgggttta aataa
//
`

var fuzzSeeds = []string{
	"@HD\tVN:1.6\n@SQ\tSN:r\tLN:45\nq1\t0\tr\t1\t60\t45M\t*\t0\t0\t" + fuzzRef + "\t*\n",
	"@SQ\tSN:r\tLN:45\nq1\t0\tr\t4\t60\t3M2I5M3D10M\t*\t0\t0\tAAACCCGGGTTTAAACCCGG\t*\nq1\t2048\tr\t20\t60\t5S10M2H\t*\t0\t0\tGGGTTTAAACCCGGG\t*\n",
	"@SQ\tSN:r\tLN:45\nq1\t0\tr\t1\t60\t10M\t*\t0\t0\tATGAAACCCG\t*\nq2\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\nq3\t256\tr\t1\t60\t10M\t*\t0\t0\tATGAAACCCG\t*\n",
	"@SQ\tSN:r\tLN:45\nq1\t0\tr\t1\t60\t3=1X5M2N3D4I10S\t*\t0\t0\tATGAAACCCGGGTTTAAACCCGG\t*\n",
	"@SQ\tSN:r\tLN:45\n@SQ\tSN:s\tLN:90\nq1\t16\ts\t50\t60\t10M\t*\t0\t0\tATGAAACCCG\t*\n",
	"@SQ\tSN:r\tLN:45\nq1\t0\tr\t40\t60\t10M\t*\t0\t0\tATGAAACCCG\t*\n",
	"@SQ\tSN:r\tLN:45\nq1\t0\tr\t1\t60\t10M\t*\t0\t0\t*\t*\n",
	"@SQ\tSN:r\tLN:45\nq1\t0\tr\t1\t60\t*\t*\t0\t0\t*\t*\n",
	"@SQ\tSN:r\tLN:45\nq1\t0\tr\t1\t60\t5H10M5P\t*\t0\t0\tATGAAACCCG\tIIIIIIIIII\n",
	"",
}

// runWithTimeout fails the test if fn does not return soon (a hang) and discards stdout while it runs.
// The limit is 10 s, or GOFASTA_FUZZ_HANG_TIMEOUT seconds if that is set.
func runWithTimeout(t *testing.T, fn func()) {
	limit := 10 * time.Second
	if s, err := strconv.Atoi(os.Getenv("GOFASTA_FUZZ_HANG_TIMEOUT")); err == nil && s > 0 {
		limit = time.Duration(s) * time.Second
	}
	done := make(chan struct{})
	go func() { fn(); close(done) }()
	select {
	case <-done:
	case <-time.After(limit):
		t.Fatalf("hang: did not return within %v", limit)
	}
}

func discardStdout(t *testing.T) {
	devnull, err := os.OpenFile(os.DevNull, os.O_WRONLY, 0)
	if err != nil {
		t.Fatal(err)
	}
	old := os.Stdout
	os.Stdout = devnull
	t.Cleanup(func() { os.Stdout = old; devnull.Close() })
}

func quietStderr(t *testing.T) {
	devnull, err := os.OpenFile(os.DevNull, os.O_WRONLY, 0)
	if err != nil {
		t.Fatal(err)
	}
	old := os.Stderr
	os.Stderr = devnull
	t.Cleanup(func() { os.Stderr = old; devnull.Close() })
}

// FuzzToMultiAlign: arbitrary bytes as a sam file must give an error or an alignment, never a panic or a hang.
func FuzzToMultiAlign(f *testing.F) {
	for _, s := range fuzzSeeds {
		f.Add([]byte(s))
	}
	f.Fuzz(func(t *testing.T, data []byte) {
		quietStderr(t)
		runWithTimeout(t, func() {
			_ = ToMultiAlign(bytes.NewReader(data), io.Discard, 60, -1, -1, false, 2, false)
		})
	})
}

func FuzzToPairAlign(f *testing.F) {
	for _, s := range fuzzSeeds {
		f.Add([]byte(s))
	}
	f.Fuzz(func(t *testing.T, data []byte) {
		discardStdout(t)
		quietStderr(t)
		runWithTimeout(t, func() {
			_ = ToPairAlign(bytes.NewReader(data), strings.NewReader(">r\n"+fuzzRef+"\n"), "stdout", 60, -1, -1, false, false, 2, false)
		})
	})
}

func FuzzVariants(f *testing.F) {
	for _, s := range fuzzSeeds {
		f.Add([]byte(s))
	}
	f.Fuzz(func(t *testing.T, data []byte) {
		quietStderr(t)
		runWithTimeout(t, func() {
			_ = Variants(bytes.NewReader(data), nil, false, strings.NewReader(fuzzGenbank), "gb", io.Discard, -1, -1, false, 0, false, 2, false)
		})
	})
}
