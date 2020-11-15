package cmd

import (
	"github.com/spf13/cobra"

	"github.com/cov-ert/gofasta/pkg/sam"
)

var toMultiAlignOutfile string
var toMultiAlignTrim bool
var toMultiAlignPad bool
var toMultiAlignTrimStart int
var toMultiAlignTrimEnd int

func init() {
	samCmd.AddCommand(toMultiAlignCmd)

	toMultiAlignCmd.Flags().StringVarP(&toMultiAlignOutfile, "fasta-out", "o", "stdout", "fasta file to write")
	toMultiAlignCmd.Flags().BoolVarP(&toMultiAlignTrim, "trim", "", false, "trim the alignment")
	toMultiAlignCmd.Flags().BoolVarP(&toMultiAlignPad, "pad", "", false, "if trim, pad the alignment with Ns")
	toMultiAlignCmd.Flags().IntVarP(&toMultiAlignTrimStart, "trimstart", "", -1, "start coordinate for trimming")
	toMultiAlignCmd.Flags().IntVarP(&toMultiAlignTrimEnd, "trimend", "", -1, "end coordinate for trimming")

	toMultiAlignCmd.Flags().SortFlags = false
}

var toMultiAlignCmd = &cobra.Command{
	Use:   "toMultiAlign",
	Short: "convert a SAM file to a multiple alignment in fasta format",
	Long:  `convert a SAM file to a multiple alignment in fasta format

		insertions relative to the reference are omitted, so all sequences
		in the output are the same ( = reference) length`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		err = sam.ToMultiAlign(samFile, toMultiAlignOutfile, toMultiAlignTrim, toMultiAlignPad, toMultiAlignTrimStart, toMultiAlignTrimEnd, samThreads)

		return
	},
}
