package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

var (
	rootCmd = &cobra.Command{
		Use:     "gofasta",
		Short:   "genomic epidemiology utilities for short genome alignments",
		Long:    `genomic epidemiology utilities for short genome alignments`,
		Version: "1.0.0",
	}
)

// Execute executes the root command.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}
