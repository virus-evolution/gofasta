package cmd

import (
	"github.com/spf13/cobra"
)

var samThreads int
var samFile string
var samReference string

func init() {
	rootCmd.AddCommand(samCmd)

	samCmd.PersistentFlags().IntVarP(&samThreads, "threads", "t", 1, "Number of threads to use")
	samCmd.PersistentFlags().StringVarP(&samFile, "samfile", "s", "", "samfile to read")
	samCmd.PersistentFlags().StringVarP(&samReference, "reference", "r", "", "reference fasta file used to generate the sam file")
}

var samCmd = &cobra.Command{
	Use:   "sam",
	Short: "Do things with sam files",
	Long:  `Do things with sam files`,
	Args:  cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {

		return nil
	},
}
