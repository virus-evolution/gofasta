package cmd

import (
	"errors"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/spf13/cobra"
)

// The helpers below check the values of flags that the user set (defaults are not checked). Each command calls
// the ones it needs from its PreRunE, which cobra runs before RunE opens any input or creates any output.

// firstError returns the first non-nil error
func firstError(errs ...error) error {
	for _, err := range errs {
		if err != nil {
			return err
		}
	}
	return nil
}

// intAtLeast checks that each named int flag, if set, is >= min
func intAtLeast(cmd *cobra.Command, min int, names ...string) error {
	for _, name := range names {
		if !cmd.Flags().Changed(name) {
			continue
		}
		v, err := cmd.Flags().GetInt(name)
		if err != nil {
			return err
		}
		if v < min {
			return fmt.Errorf("--%s must be at least %d (got %d)", name, min, v)
		}
	}
	return nil
}

// floatAtLeast checks that each named float64 or float32 flag, if set, is a finite number >= min
func floatAtLeast(cmd *cobra.Command, min float64, names ...string) error {
	for _, name := range names {
		f := cmd.Flags().Lookup(name)
		if f == nil || !f.Changed {
			continue
		}
		v, err := strconv.ParseFloat(f.Value.String(), 64)
		if err != nil {
			return fmt.Errorf("couldn't parse --%s as a number: %v", name, err)
		}
		if math.IsNaN(v) || math.IsInf(v, 0) || v < min {
			return fmt.Errorf("--%s must be a finite number >= %v (got %s)", name, min, f.Value.String())
		}
	}
	return nil
}

// startNotAfterEnd checks that, if both --start and --end are given, start <= end
func startNotAfterEnd(cmd *cobra.Command) error {
	if !cmd.Flags().Changed("start") || !cmd.Flags().Changed("end") {
		return nil
	}
	start, err := cmd.Flags().GetInt("start")
	if err != nil {
		return err
	}
	end, err := cmd.Flags().GetInt("end")
	if err != nil {
		return err
	}
	if start > end {
		return fmt.Errorf("--start must not be after --end (got %d and %d)", start, end)
	}
	return nil
}

// limitAtLeast checks that each named flag, if set, is a number >= min; unlike floatAtLeast it accepts +Inf,
// which means "no limit"
func limitAtLeast(cmd *cobra.Command, min float64, names ...string) error {
	for _, name := range names {
		f := cmd.Flags().Lookup(name)
		if f == nil || !f.Changed {
			continue
		}
		v, err := strconv.ParseFloat(f.Value.String(), 64)
		if err != nil {
			return fmt.Errorf("couldn't parse --%s as a number: %v", name, err)
		}
		if math.IsNaN(v) || v < min {
			return fmt.Errorf("--%s must be a number >= %v, or inf for no limit (got %s)", name, min, f.Value.String())
		}
	}
	return nil
}

// outputNotInput checks that the file named by the output flag is none of the files named by the input flags:
// creating the output would truncate that input before it is read
func outputNotInput(cmd *cobra.Command, output string, inputs ...string) error {
	o := cmd.Flags().Lookup(output)
	if o == nil {
		return nil
	}
	for _, input := range inputs {
		i := cmd.Flags().Lookup(input)
		if i != nil && sameFile(o.Value.String(), i.Value.String()) {
			return fmt.Errorf("--%s and --%s are the same file; writing the output would overwrite the input", output, input)
		}
	}
	return nil
}

// sameFile is true if both are names of the same existing file (stdin/stdout/empty names are never files)
func sameFile(a, b string) bool {
	for _, n := range []string{a, b} {
		if n == "" || n == "stdin" || n == "stdout" {
			return false
		}
	}
	fa, err := os.Stat(a)
	if err != nil {
		return false
	}
	fb, err := os.Stat(b)
	if err != nil {
		return false
	}
	return os.SameFile(fa, fb)
}

// annotationFormat tells a GenBank file ("gb") from a GFF3 file ("gff") by its suffix, ignoring case.
func annotationFormat(path string) (string, error) {
	switch strings.ToLower(filepath.Ext(path)) {
	case ".gb", ".gbk", ".gbff":
		return "gb", nil
	case ".gff", ".gff3":
		return "gff", nil
	}
	return "", errors.New("couldn't tell if --annotation was a GenBank (.gb, .gbk, .gbff) or a GFF3 (.gff, .gff3) file")
}
