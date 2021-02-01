package cmd

import (
	"github.com/spf13/cobra"

	"github.com/cov-ert/gofasta/pkg/sam"
)

var toPairAlignGenbankFile string
var toPairAlignGenbankFeature string
var toPairAlignOutpath string
var toPairAlignOmitReference bool
var toPairAlignSkipInsertions bool

func init() {
	samCmd.AddCommand(toPairAlignCmd)

	toPairAlignCmd.Flags().StringVarP(&toPairAlignGenbankFile, "genbank", "g", "", "Genbank format annotation of a sequence in the same coordinates as the alignment")
	toPairAlignCmd.Flags().StringVarP(&toPairAlignGenbankFeature, "feature", "", "", "Feature to output (choose one of: gene, CDS). If none is specified, will output the entire alignment")
	toPairAlignCmd.Flags().StringVarP(&toPairAlignOutpath, "outpath", "o", "", "Output path where fasta files will be written")
	toPairAlignCmd.Flags().BoolVarP(&toPairAlignOmitReference, "omit-reference", "", false, "Omit the reference sequences from the output alignments")
	toPairAlignCmd.Flags().BoolVarP(&toPairAlignSkipInsertions, "skip-insertions", "", false, "Skip insertions relative to the reference")

	toPairAlignCmd.Flags().Lookup("omit-reference").NoOptDefVal = "true"
	toPairAlignCmd.Flags().Lookup("skip-insertions").NoOptDefVal = "true"

	toPairAlignCmd.Flags().SortFlags = false
}

var toPairAlignCmd = &cobra.Command{
	Use:   "toPairAlign",
	Short: "convert a SAM file to pairwise alignments in fasta format",
	Long:  `convert a SAM file to pairwise alignments in fasta format`,

	RunE: func(cmd *cobra.Command, args []string) (err error) {

		err = sam.ToPairAlign(samFile, samReference, toPairAlignGenbankFile, toPairAlignGenbankFeature, toPairAlignOutpath, toPairAlignOmitReference, toPairAlignSkipInsertions, samThreads)

		return err
	},
}
