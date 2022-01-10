package cmd

import (
	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/sam"
)

var toPairAlignOutpath string
var toPairAlignOmitReference bool
var toPairAlignSkipInsertions bool
var toPairAlignTrim bool
var toPairAlignTrimStart int
var toPairAlignTrimEnd int

func init() {
	samCmd.AddCommand(toPairAlignCmd)

	toPairAlignCmd.Flags().StringVarP(&toPairAlignOutpath, "outpath", "o", "", "Output path where fasta files will be written")
	toPairAlignCmd.Flags().BoolVarP(&toPairAlignOmitReference, "omit-reference", "", false, "Omit the reference sequences from the output alignments")
	toPairAlignCmd.Flags().BoolVarP(&toPairAlignSkipInsertions, "skip-insertions", "", false, "Skip insertions relative to the reference from the output alignments")
	toPairAlignCmd.Flags().BoolVarP(&toPairAlignTrim, "trim", "", false, "Trim the alignment (to reference coordinates")
	toPairAlignCmd.Flags().IntVarP(&toPairAlignTrimStart, "trimstart", "", -1, "Start coordinate for trimming (0-based, half open)")
	toPairAlignCmd.Flags().IntVarP(&toPairAlignTrimEnd, "trimend", "", -1, "End coordinate for trimming (0-based, half open)")

	toPairAlignCmd.Flags().Lookup("omit-reference").NoOptDefVal = "true"
	toPairAlignCmd.Flags().Lookup("skip-insertions").NoOptDefVal = "true"

	toPairAlignCmd.Flags().SortFlags = false
}

var toPairAlignCmd = &cobra.Command{
	Use:     "toPairAlign",
	Aliases: []string{"topairalign", "topa"},
	Short:   "convert a SAM file to pairwise alignments in fasta format",
	Long:    `convert a SAM file to pairwise alignments in fasta format`,

	RunE: func(cmd *cobra.Command, args []string) (err error) {

		samIn, err := gfio.OpenIn(*cmd.Flag("samfile"))
		if err != nil {
			return err
		}
		defer samIn.Close()

		ref, err := gfio.OpenIn(*cmd.Flag("reference"))
		if err != nil {
			return err
		}
		defer ref.Close()

		err = sam.ToPairAlign(samIn, ref, toPairAlignOutpath, toPairAlignTrim, toPairAlignTrimStart, toPairAlignTrimEnd, toPairAlignOmitReference, toPairAlignSkipInsertions, samThreads)

		return err
	},
}
