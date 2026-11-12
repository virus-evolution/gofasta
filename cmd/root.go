package cmd

import (
	"os"

	"github.com/spf13/cobra"
)

var (
	rootCmd = &cobra.Command{
		Use:   "gofasta",
		Short: "Command-line utilities for genomic epidemiology research",
		Long: `Command-line utilities for genomic epidemiology research
		
If you use gofasta in your work, please cite:

Jackson B (2022). gofasta: command-line utilities for genomic epidemiology research. Bioinformatics 38 (16), 4033-4035
https://doi.org/10.1093/bioinformatics/btac424`,
		Version: "1.2.3",
	}
)

// strict makes a record that cannot be processed stop the run with an error, instead of a warning
var strict bool

func init() {
	rootCmd.PersistentFlags().BoolVar(&strict, "strict", false, "Stop with an error at the first record that can't be processed, instead of warning and skipping it")
}

// Execute executes the root command.
func Execute() {
	// cobra prints the error (and usage) to stderr itself; stdout is reserved for data
	if err := rootCmd.Execute(); err != nil {
		os.Exit(1)
	}
}
