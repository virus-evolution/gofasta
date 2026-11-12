package cmd

import (
	"bytes"
	"context"
	"errors"
	"flag"
	"fmt"
	"io/fs"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"
	"time"
)

// End-to-end tests: each case runs gofasta (this test binary, re-executed as the command line
// program) in a fresh directory holding copies of testdata/e2e/inputs, and compares its exit
// code, stdout, stderr and any files it writes with testdata/e2e/golden/<case>.txt (only the
// exit code and stderr when it fails).
// Regenerate the golden files with: go test ./cmd -run TestEndToEnd -update
// The inputs are written by testdata/e2e/make_inputs.py (needs Biopython), except query.csv and
// target.csv, which are the output of gofasta updown list for query.fasta and target.fasta, and
// target_bad.csv, which is target.csv with two malformed lines added.

var updateGolden = flag.Bool("update", false, "rewrite the end-to-end golden files")

const e2eExecEnv = "GOFASTA_E2E_EXEC"

func TestMain(m *testing.M) {
	if os.Getenv(e2eExecEnv) == "1" {
		Execute()
		os.Exit(0)
	}
	os.Exit(m.Run())
}

type e2eCase struct {
	name  string
	args  []string
	stdin string // input file to use as stdin, if any
	// threadsFlag, if set, names a threads flag: the case is also run with it set to 4 and
	// must give exactly the same result
	threadsFlag string
}

var e2eCases = []e2eCase{
	// snps
	{name: "snps", args: []string{"snps", "-r", "ref.fasta", "-q", "msa.fasta"}},
	{name: "snps_stdin", args: []string{"snps", "-r", "ref.fasta"}, stdin: "msa.fasta"},
	{name: "snps_hardgaps", args: []string{"snps", "-r", "ref.fasta", "-q", "msa.fasta", "--hard-gaps"}},
	{name: "snps_aggregate", args: []string{"snps", "-r", "ref.fasta", "-q", "msa.fasta", "--aggregate"}},
	{name: "snps_aggregate_threshold", args: []string{"snps", "-r", "ref.fasta", "-q", "msa.fasta", "--aggregate", "--threshold", "0.3"}},
	{name: "snps_outfile", args: []string{"snps", "-r", "ref.fasta", "-q", "msa.fasta", "-o", "out.csv"}},
	{name: "snps_badnuc", args: []string{"snps", "-r", "ref.fasta", "-q", "msa_badnuc.fasta"}},
	{name: "snps_badrow", args: []string{"snps", "-r", "ref.fasta", "-q", "msa_badrow.fasta"}},
	{name: "snps_badnuc_strict", args: []string{"snps", "-r", "ref.fasta", "-q", "msa_badnuc.fasta", "--strict"}},
	{name: "snps_badrow_strict", args: []string{"snps", "-r", "ref.fasta", "-q", "msa_badrow.fasta", "--strict"}},
	{name: "snps_aggregate_badrow", args: []string{"snps", "-r", "ref.fasta", "-q", "msa_badrow.fasta", "--aggregate"}},
	{name: "snps_all_bad", args: []string{"snps", "-r", "ref_short.fasta", "-q", "msa.fasta"}},
	{name: "snps_empty", args: []string{"snps", "-r", "ref.fasta", "-q", "empty.fasta"}},
	{name: "snps_missing_file", args: []string{"snps", "-r", "ref.fasta", "-q", "nosuchfile.fasta"}},

	// variants
	{name: "variants_gb", args: []string{"variants", "--msa", "msa.fasta", "-r", "ref", "-a", "ref.gb"}, threadsFlag: "-t"},
	{name: "variants_gff", args: []string{"variants", "--msa", "msa.fasta", "-r", "ref", "-a", "ref.gff"}, threadsFlag: "-t"},
	{name: "variants_stdin", args: []string{"variants", "-r", "ref", "-a", "ref.gb"}, stdin: "msa.fasta"},
	{name: "variants_append_snps", args: []string{"variants", "--msa", "msa.fasta", "-r", "ref", "-a", "ref.gb", "--append-snps"}},
	{name: "variants_aggregate", args: []string{"variants", "--msa", "msa.fasta", "-r", "ref", "-a", "ref.gb", "--aggregate"}, threadsFlag: "-t"},
	{name: "variants_aggregate_threshold", args: []string{"variants", "--msa", "msa.fasta", "-r", "ref", "-a", "ref.gb", "--aggregate", "--threshold", "0.2", "--append-snps"}},
	{name: "variants_start_end", args: []string{"variants", "--msa", "msa.fasta", "-r", "ref", "-a", "ref.gb", "--start", "10", "--end", "90"}},
	{name: "variants_insertions", args: []string{"variants", "--msa", "msa_ins.fasta", "-r", "ref", "-a", "ref.gb"}, threadsFlag: "-t"},
	{name: "variants_outfile", args: []string{"variants", "--msa", "msa.fasta", "-r", "ref", "-a", "ref.gb", "-o", "out.csv"}},
	{name: "variants_ref_not_found", args: []string{"variants", "--msa", "msa.fasta", "-r", "nosuchref", "-a", "ref.gb"}},
	{name: "variants_badrow", args: []string{"variants", "--msa", "msa_badrow.fasta", "-r", "ref", "-a", "ref.gb"}, threadsFlag: "-t"},
	{name: "variants_badnuc", args: []string{"variants", "--msa", "msa_badnuc.fasta", "-r", "ref", "-a", "ref.gb"}, threadsFlag: "-t"},
	{name: "variants_badnuc_strict", args: []string{"variants", "--msa", "msa_badnuc.fasta", "-r", "ref", "-a", "ref.gb", "--strict"}},
	{name: "variants_badrow_strict", args: []string{"variants", "--msa", "msa_badrow.fasta", "-r", "ref", "-a", "ref.gb", "--strict"}},
	{name: "variants_aggregate_badrow", args: []string{"variants", "--msa", "msa_badrow.fasta", "-r", "ref", "-a", "ref.gb", "--aggregate"}, threadsFlag: "-t"},
	{name: "variants_ref_after_bad", args: []string{"variants", "--msa", "msa_badnuc.fasta", "-r", "q8", "-a", "ref.gb"}},
	{name: "variants_stdin_badnuc", args: []string{"variants", "-r", "ref", "-a", "ref.gb"}, stdin: "msa_badnuc.fasta"},
	{name: "variants_empty", args: []string{"variants", "--msa", "empty.fasta", "-r", "ref", "-a", "ref.gb"}},
	{name: "variants_bad_suffix", args: []string{"variants", "--msa", "msa.fasta", "-r", "ref", "-a", "ignore.txt"}},

	// closest
	{name: "closest_n", args: []string{"closest", "--query", "query.fasta", "--target", "target.fasta", "-n", "2"}, threadsFlag: "-t"},
	{name: "closest_default", args: []string{"closest", "--query", "query.fasta", "--target", "target.fasta"}, threadsFlag: "-t"},
	{name: "closest_d", args: []string{"closest", "--query", "query.fasta", "--target", "target.fasta", "-d", "0.02"}, threadsFlag: "-t"},
	{name: "closest_d_inf", args: []string{"closest", "--query", "query.fasta", "--target", "target.fasta", "-d", "inf"}},
	{name: "closest_n_table", args: []string{"closest", "--query", "query.fasta", "--target", "target.fasta", "-n", "2", "--table"}},
	{name: "closest_snp", args: []string{"closest", "--query", "query.fasta", "--target", "target.fasta", "-m", "snp", "-n", "3"}, threadsFlag: "-t"},
	{name: "closest_tn93", args: []string{"closest", "--query", "query.fasta", "--target", "target.fasta", "-m", "tn93", "-n", "3"}, threadsFlag: "-t"},
	{name: "closest_badnuc_target", args: []string{"closest", "--query", "query.fasta", "--target", "msa_badnuc.fasta", "-n", "1"}},
	{name: "closest_badnuc_target_strict", args: []string{"closest", "--query", "query.fasta", "--target", "msa_badnuc.fasta", "-n", "1", "--strict"}},
	{name: "closest_badrow_target", args: []string{"closest", "--query", "query.fasta", "--target", "msa_badrow.fasta", "-n", "2"}, threadsFlag: "-t"},
	{name: "closest_bad_queries", args: []string{"closest", "--query", "msa_badrow.fasta", "--target", "target.fasta"}, threadsFlag: "-t"},
	{name: "closest_bad_queries_strict", args: []string{"closest", "--query", "msa_badnuc.fasta", "--target", "target.fasta", "--strict"}},
	{name: "closest_empty_query", args: []string{"closest", "--query", "empty.fasta", "--target", "target.fasta", "-n", "1"}},

	// updown
	{name: "updown_list", args: []string{"updown", "list", "-r", "ref.fasta", "-q", "msa.fasta"}},
	{name: "updown_list_stdin", args: []string{"updown", "list", "-r", "ref.fasta"}, stdin: "msa.fasta"},
	{name: "updown_list_badnuc", args: []string{"updown", "list", "-r", "ref.fasta", "-q", "msa_badnuc.fasta"}},
	{name: "updown_list_badnuc_strict", args: []string{"updown", "list", "-r", "ref.fasta", "-q", "msa_badnuc.fasta", "--strict"}},
	{name: "updown_topranking_bad_target", args: []string{"updown", "topranking", "-r", "ref.fasta", "-q", "query.fasta", "-t", "msa_badrow.fasta", "--size-total", "4"}},
	{name: "updown_topranking_bad_csv_target", args: []string{"updown", "topranking", "-q", "query.csv", "-t", "target_bad.csv", "--size-total", "4"}},
	{name: "updown_topranking_bad_csv_query", args: []string{"updown", "topranking", "-q", "target_bad.csv", "-t", "target.csv", "--size-total", "4"}},
	{name: "updown_topranking_bad_csv_strict", args: []string{"updown", "topranking", "-q", "query.csv", "-t", "target_bad.csv", "--size-total", "4", "--strict"}},
	{name: "updown_topranking_fasta", args: []string{"updown", "topranking", "-r", "ref.fasta", "-q", "query.fasta", "-t", "target.fasta", "--size-total", "4"}},
	{name: "updown_topranking_table", args: []string{"updown", "topranking", "-r", "ref.fasta", "-q", "query.fasta", "-t", "target.fasta", "--size-total", "4", "--table"}},
	{name: "updown_topranking_ignore", args: []string{"updown", "topranking", "-r", "ref.fasta", "-q", "query.fasta", "-t", "target.fasta", "--size-total", "4", "--ignore", "ignore.txt"}},
	{name: "updown_topranking_csv", args: []string{"updown", "topranking", "-q", "query.csv", "-t", "target.csv", "--size-total", "4"}},
	{name: "updown_topranking_csv_query", args: []string{"updown", "topranking", "-r", "ref.fasta", "-q", "query.csv", "-t", "target.fasta", "--size-total", "4"}},
	{name: "updown_topranking_csv_target", args: []string{"updown", "topranking", "-r", "ref.fasta", "-q", "query.fasta", "-t", "target.csv", "--size-total", "4", "--table"}},
	{name: "updown_topranking_bins", args: []string{"updown", "topranking", "-r", "ref.fasta", "-q", "query.fasta", "-t", "target.fasta", "--size-up", "1", "--size-down", "1", "--size-side", "1", "--size-same", "1", "--dist-all", "3"}},

	// sam toMultiAlign
	{name: "sam_toma", args: []string{"sam", "toMultiAlign", "-s", "aln.sam"}, threadsFlag: "-t"},
	{name: "sam_toma_stdin", args: []string{"sam", "toMultiAlign"}, stdin: "aln.sam"},
	{name: "sam_toma_pad", args: []string{"sam", "toMultiAlign", "-s", "aln.sam", "--pad"}},
	{name: "sam_toma_start_end_pad", args: []string{"sam", "toMultiAlign", "-s", "aln.sam", "--start", "20", "--end", "100", "--pad"}},
	{name: "sam_toma_start_end", args: []string{"sam", "toMultiAlign", "-s", "aln.sam", "--start", "20", "--end", "100"}},
	{name: "sam_toma_wrap", args: []string{"sam", "toMultiAlign", "-s", "aln.sam", "-w", "50"}},
	{name: "sam_toma_outfile", args: []string{"sam", "toMultiAlign", "-s", "aln.sam", "-o", "out.fasta"}},
	{name: "sam_toma_badrecord", args: []string{"sam", "toMultiAlign", "-s", "aln_badrecord.sam"}, threadsFlag: "-t"},
	{name: "sam_toma_badrecord_strict", args: []string{"sam", "toMultiAlign", "-s", "aln_badrecord.sam", "--strict"}},
	{name: "sam_toma_badquery", args: []string{"sam", "toMultiAlign", "-s", "aln_badquery.sam"}, threadsFlag: "-t"},
	{name: "sam_toma_badquery_strict", args: []string{"sam", "toMultiAlign", "-s", "aln_badquery.sam", "--strict"}},
	{name: "sam_toma_untidy", args: []string{"sam", "toMultiAlign", "-s", "aln_untidy.sam"}},
	{name: "sam_toma_twosq", args: []string{"sam", "toMultiAlign", "-s", "aln_twosq.sam"}},
	{name: "sam_toma_headeronly", args: []string{"sam", "toMultiAlign", "-s", "aln_headeronly.sam"}},
	{name: "sam_toma_empty", args: []string{"sam", "toMultiAlign", "-s", "empty.sam"}},

	// sam toPairAlign
	{name: "sam_topa_stdout", args: []string{"sam", "toPairAlign", "-s", "aln.sam", "-r", "ref.fasta", "-o", "stdout"}, threadsFlag: "-t"},
	{name: "sam_topa_omitref", args: []string{"sam", "toPairAlign", "-s", "aln.sam", "-r", "ref.fasta", "-o", "stdout", "--omit-reference"}},
	{name: "sam_topa_skipins", args: []string{"sam", "toPairAlign", "-s", "aln.sam", "-r", "ref.fasta", "-o", "stdout", "--skip-insertions"}},
	{name: "sam_topa_start_end", args: []string{"sam", "toPairAlign", "-s", "aln.sam", "-r", "ref.fasta", "-o", "stdout", "--start", "20", "--end", "100"}},
	{name: "sam_topa_dir", args: []string{"sam", "toPairAlign", "-s", "aln.sam", "-r", "ref.fasta", "-o", "pairs", "-w", "60"}, threadsFlag: "-t"},
	{name: "sam_topa_badrecord", args: []string{"sam", "toPairAlign", "-s", "aln_badrecord.sam", "-r", "ref.fasta", "-o", "stdout"}},
	{name: "sam_topa_empty", args: []string{"sam", "toPairAlign", "-s", "empty.sam", "-r", "ref.fasta", "-o", "stdout"}},

	// sam variants
	{name: "sam_variants_gb", args: []string{"sam", "variants", "-s", "aln.sam", "-r", "ref.fasta", "-a", "ref.gb"}, threadsFlag: "-t"},
	{name: "sam_variants_gff", args: []string{"sam", "variants", "-s", "aln.sam", "-r", "ref.fasta", "-a", "ref.gff"}, threadsFlag: "-t"},
	{name: "sam_variants_aggregate", args: []string{"sam", "variants", "-s", "aln.sam", "-r", "ref.fasta", "-a", "ref.gb", "--aggregate", "--append-snps"}, threadsFlag: "-t"},
	{name: "sam_variants_start_end", args: []string{"sam", "variants", "-s", "aln.sam", "-r", "ref.fasta", "-a", "ref.gb", "--start", "20", "--end", "90"}},
	{name: "sam_variants_badrecord", args: []string{"sam", "variants", "-s", "aln_badrecord.sam", "-r", "ref.fasta", "-a", "ref.gb"}, threadsFlag: "-t"},
	{name: "sam_variants_badquery_aggregate", args: []string{"sam", "variants", "-s", "aln_badquery.sam", "-r", "ref.fasta", "-a", "ref.gb", "--aggregate"}, threadsFlag: "-t"},
	{name: "sam_topa_badquery_dir", args: []string{"sam", "toPairAlign", "-s", "aln_badquery.sam", "-r", "ref.fasta", "-o", "pairs"}},
	{name: "sam_variants_twosq", args: []string{"sam", "variants", "-s", "aln_twosq.sam", "-r", "ref.fasta", "-a", "ref.gb"}},
	{name: "sam_variants_empty", args: []string{"sam", "variants", "-s", "empty.sam", "-r", "ref.fasta", "-a", "ref.gb"}},
}

// runE2E runs one command line in a fresh copy of the inputs and renders the result.
func runE2E(t *testing.T, args []string, stdin string) string {
	t.Helper()
	inputs, err := filepath.Abs(filepath.Join("testdata", "e2e", "inputs"))
	if err != nil {
		t.Fatal(err)
	}
	dir := t.TempDir()
	before := map[string]bool{}
	entries, err := os.ReadDir(inputs)
	if err != nil {
		t.Fatal(err)
	}
	for _, e := range entries {
		b, err := os.ReadFile(filepath.Join(inputs, e.Name()))
		if err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(filepath.Join(dir, e.Name()), b, 0o644); err != nil {
			t.Fatal(err)
		}
		before[e.Name()] = true
	}

	self, err := os.Executable()
	if err != nil {
		t.Fatal(err)
	}
	ctx, cancel := context.WithTimeout(context.Background(), 60*time.Second)
	defer cancel()
	c := exec.CommandContext(ctx, self, args...)
	c.Dir = dir
	c.Env = append(os.Environ(), e2eExecEnv+"=1")
	var stdout, stderr bytes.Buffer
	c.Stdout, c.Stderr = &stdout, &stderr
	if stdin != "" {
		f, err := os.Open(filepath.Join(dir, stdin))
		if err != nil {
			t.Fatal(err)
		}
		defer f.Close()
		c.Stdin = f
	}
	err = c.Run()
	if ctx.Err() != nil {
		t.Fatalf("%v: timed out", args)
	}
	code := 0
	var exitErr *exec.ExitError
	if errors.As(err, &exitErr) {
		code = exitErr.ExitCode()
	} else if err != nil {
		t.Fatal(err)
	}

	var b strings.Builder
	fmt.Fprintf(&b, "args: %s\n", strings.Join(args, " "))
	if stdin != "" {
		fmt.Fprintf(&b, "stdin: %s\n", stdin)
	}
	fmt.Fprintf(&b, "exit: %d\n", code)
	if code != 0 {
		// how much output is written before a failure depends on goroutine timing, so only the
		// exit code and stderr are compared
		fmt.Fprintf(&b, "--- stderr\n%s", stderr.String())
		return b.String()
	}
	fmt.Fprintf(&b, "--- stdout\n%s", stdout.String())
	fmt.Fprintf(&b, "--- stderr\n%s", stderr.String())

	// any file that was not an input is an output
	var outputs []string
	err = filepath.WalkDir(dir, func(path string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		rel, _ := filepath.Rel(dir, path)
		if d.IsDir() || before[rel] {
			return nil
		}
		outputs = append(outputs, rel)
		return nil
	})
	if err != nil {
		t.Fatal(err)
	}
	sort.Strings(outputs)
	for _, o := range outputs {
		content, err := os.ReadFile(filepath.Join(dir, o))
		if err != nil {
			t.Fatal(err)
		}
		fmt.Fprintf(&b, "--- file %s\n%s", filepath.ToSlash(o), content)
	}
	return b.String()
}

func TestEndToEnd(t *testing.T) {
	seen := map[string]bool{}
	for _, tc := range e2eCases {
		if seen[tc.name] {
			t.Fatalf("duplicate case name %s", tc.name)
		}
		seen[tc.name] = true
		t.Run(tc.name, func(t *testing.T) {
			t.Parallel()
			got := runE2E(t, tc.args, tc.stdin)
			golden := filepath.Join("testdata", "e2e", "golden", tc.name+".txt")
			if *updateGolden {
				if err := os.MkdirAll(filepath.Dir(golden), 0o755); err != nil {
					t.Fatal(err)
				}
				if err := os.WriteFile(golden, []byte(got), 0o644); err != nil {
					t.Fatal(err)
				}
			} else {
				want, err := os.ReadFile(golden)
				if err != nil {
					t.Fatalf("%v (run with -update to create it)", err)
				}
				if got != string(want) {
					t.Errorf("output differs from %s:\n%s", golden, lineDiff(string(want), got))
				}
			}
			if tc.threadsFlag != "" {
				args := append(append([]string{}, tc.args...), tc.threadsFlag, "4")
				threaded := runE2E(t, args, tc.stdin)
				// the recorded command line differs; everything after it must not
				_, a, _ := strings.Cut(got, "\n")
				_, b, _ := strings.Cut(threaded, "\n")
				if a != b {
					t.Errorf("%s 4 gives a different result:\n%s", tc.threadsFlag, lineDiff(a, b))
				}
			}
		})
	}
}

// lineDiff shows the first differing line of want and got, with a little context.
func lineDiff(want, got string) string {
	w, g := strings.Split(want, "\n"), strings.Split(got, "\n")
	for i := 0; i < len(w) || i < len(g); i++ {
		var wl, gl string
		if i < len(w) {
			wl = w[i]
		}
		if i < len(g) {
			gl = g[i]
		}
		if wl != gl {
			return fmt.Sprintf("line %d:\n  want: %q\n  got:  %q", i+1, wl, gl)
		}
	}
	return "(no line differs; check line endings)"
}
