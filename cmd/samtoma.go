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

	toMultiAlignCmd.Flags().StringVarP(&toMultiAlignOutfile, "fasta-out", "o", "stdout", "Where to write the alignment")
	toMultiAlignCmd.Flags().BoolVarP(&toMultiAlignTrim, "trim", "", false, "Trim the alignment")
	toMultiAlignCmd.Flags().BoolVarP(&toMultiAlignPad, "pad", "", false, "If trim, replace the trimmed regions with Ns")
	toMultiAlignCmd.Flags().IntVarP(&toMultiAlignTrimStart, "trimstart", "", -1, "Start coordinate for trimming (0-based, half open)")
	toMultiAlignCmd.Flags().IntVarP(&toMultiAlignTrimEnd, "trimend", "", -1, "End coordinate for trimming (0-based, half open)")

	toMultiAlignCmd.Flags().SortFlags = false
}

var toMultiAlignCmd = &cobra.Command{
	Use:     "toMultiAlign",
	Aliases: []string{"tomultialign"},
	Short:   "Convert a SAM file to a multiple alignment in fasta format",
	Long: `Convert a SAM file to a multiple alignment in fasta format

Insertions relative to the reference are omitted, so all sequences
in the output are the same ( = reference) length.

Example usage:
	gofasta sam toMultiAlign -s aligned.sam -o aligned.fasta

If you want, you can trim (and optionally pad) the output alignment to coordinates of your choosing:
	gofasta sam toMultiAlign -s aligned.sam --trim --trimstart 250 --trimend 29000 --pad -o aligned.fasta

If input and output files are not specified, the behaviour is to read the sam file from stdin and write
the fasta file to stdout, e.g.:
	minimap2 -a -x asm5 reference.fasta unaligned.fasta | gofasta sam toMultiAlign > aligned.fasta`,

	RunE: func(cmd *cobra.Command, args []string) (err error) {

		err = sam.ToMultiAlign(samFile, samReference, toMultiAlignOutfile, toMultiAlignTrim, toMultiAlignPad, toMultiAlignTrimStart, toMultiAlignTrimEnd, samThreads)

		return
	},
}
