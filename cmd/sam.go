package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

var samThreads int
var samFile string
var samReference string

func init() {
	rootCmd.AddCommand(samCmd)

	samCmd.PersistentFlags().IntVarP(&samThreads, "threads", "t", 1, "Number of threads to use")
	samCmd.PersistentFlags().StringVarP(&samFile, "samfile", "s", "stdin", "Samfile to read. If none is specified, will read from stdin")
	samCmd.PersistentFlags().StringVarP(&samReference, "reference", "r", "", "Reference fasta file used to generate the sam file")
}

var samCmd = &cobra.Command{
	Use:   "sam",
	Short: "Do things with sam files",
	Long:  `Do things with sam files`,
	Args:  cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		// Args is MinimumNArgs(1), so there is always an argument here, and it is not a subcommand
		// (cobra would have run that instead)
		return fmt.Errorf("unknown command %q for %q", args[0], cmd.CommandPath())
	},
}
