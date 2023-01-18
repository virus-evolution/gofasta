package cmd

import (
	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/sam"
)

var toPairAlignOutpath string
var toPairAlignOmitReference bool
var toPairAlignSkipInsertions bool
var toPairAlignStart int
var toPairAlignEnd int
var toPairAlignWrap int

func init() {
	samCmd.AddCommand(toPairAlignCmd)

	toPairAlignCmd.Flags().StringVarP(&toPairAlignOutpath, "outpath", "o", "", "Output path where fasta files will be written")
	toPairAlignCmd.Flags().BoolVarP(&toPairAlignOmitReference, "omit-reference", "", false, "Omit the reference sequences from the output alignments")
	toPairAlignCmd.Flags().BoolVarP(&toPairAlignSkipInsertions, "skip-insertions", "", false, "Skip insertions relative to the reference from the output alignments")
	toPairAlignCmd.Flags().IntVarP(&toPairAlignStart, "start", "", -1, "1-based first nucleotide position (in reference coordinates) to retain in the output. Bases before this position are omitted")
	toPairAlignCmd.Flags().IntVarP(&toPairAlignEnd, "end", "", -1, "1-based last nucleotide position (in reference coordinates) to retain in the output. Bases after this position are omitted")
	toPairAlignCmd.Flags().IntVarP(&toPairAlignWrap, "wrap", "w", -1, "Wrap the output alignment to this number of nucleotides wide. Omit this option not to wrap the output.")

	toPairAlignCmd.Flags().Lookup("wrap").NoOptDefVal = "80"

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

		err = sam.ToPairAlign(samIn, ref, toPairAlignOutpath, toPairAlignWrap, toPairAlignStart, toPairAlignEnd, toPairAlignOmitReference, toPairAlignSkipInsertions, samThreads)

		return err
	},
}
