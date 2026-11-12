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

	toPairAlignCmd.Flags().StringVarP(&toPairAlignOutpath, "outpath", "o", "", "Directory to write one fasta file per query to, or stdout to write all the alignments to stdout (required)")
	toPairAlignCmd.Flags().BoolVarP(&toPairAlignOmitReference, "omit-reference", "", false, "Omit the reference sequences from the output alignments")
	toPairAlignCmd.Flags().BoolVarP(&toPairAlignSkipInsertions, "skip-insertions", "", false, "Skip insertions relative to the reference from the output alignments")
	toPairAlignCmd.Flags().IntVarP(&toPairAlignStart, "start", "", -1, "1-based first nucleotide position (in reference coordinates) to retain in the output. Bases before this position are omitted")
	toPairAlignCmd.Flags().IntVarP(&toPairAlignEnd, "end", "", -1, "1-based last nucleotide position (in reference coordinates) to retain in the output. Bases after this position are omitted")
	toPairAlignCmd.Flags().IntVarP(&toPairAlignWrap, "wrap", "w", -1, "Wrap the output alignment to this number of nucleotides wide. Omit this option not to wrap the output.")

	toPairAlignCmd.Flags().Lookup("omit-reference").NoOptDefVal = "true"
	toPairAlignCmd.Flags().Lookup("skip-insertions").NoOptDefVal = "true"

	toPairAlignCmd.Flags().BoolVar(getStrictBoolVarArgs())
	toPairAlignCmd.MarkFlagRequired("outpath")
	toPairAlignCmd.Flags().SortFlags = false
}

var toPairAlignCmd = &cobra.Command{
	Use:     "toPairAlign",
	Aliases: []string{"topairalign", "topa"},
	Short:   "Convert a SAM file to pairwise alignments in fasta format",
	Long: `Convert a SAM file to pairwise alignments in fasta format

Each query is written as a pairwise alignment with the reference: the reference sequence first, then the query.
Unlike toMultiAlign, insertions relative to the reference are kept (unless you use --skip-insertions), so the output
alignment can be longer than the reference. External alignment gaps are always written as N.

--reference must be the fasta file that was used to generate the sam file, and its sequence must be the same
length as the @SQ line in the sam header.

Example usage:
	gofasta sam toPairAlign -s aligned.sam -r reference.fasta -o alignments

writes one file per query to the directory alignments (which is created if it doesn't exist). Each file is
named after its query, with any "/" replaced by "_", plus ".fasta". To write all the alignments to stdout
instead, use -o stdout, e.g.:
	minimap2 -a -x asm20 --score-N=0 reference.fasta unaligned.fasta | gofasta sam toPairAlign -r reference.fasta -o stdout`,

	PreRunE: func(cmd *cobra.Command, args []string) error {
		return intAtLeast(cmd, 1, "threads")
	},
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

		err = sam.ToPairAlign(samIn, ref, toPairAlignOutpath, toPairAlignWrap, toPairAlignStart, toPairAlignEnd, toPairAlignOmitReference, toPairAlignSkipInsertions, samThreads, strict)

		return err
	},
}
