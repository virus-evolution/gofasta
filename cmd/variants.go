package cmd

import (
	"github.com/spf13/cobra"

	"github.com/cov-ert/gofasta/pkg/sam"
)

var variantGenbankFile string
var variantOutfile string
// var variantSkipInsertions bool

func init() {
	samCmd.AddCommand(variantCmd)

	variantCmd.Flags().StringVarP(&variantGenbankFile, "genbank", "g", "", "Genbank format annotation of a sequence in the same coordinates as the alignment")
	variantCmd.Flags().StringVarP(&variantOutfile, "outfile", "o", "stdout", "Where to write the variants")
	// variantCmd.Flags().BoolVarP(&variantSkipInsertions, "skip-insertions", "", false, "skip insertions relative to the reference")

	// variantCmd.Flags().Lookup("skip-insertions").NoOptDefVal = "true"

	variantCmd.Flags().SortFlags = false
}

var variantCmd = &cobra.Command{
	Use:   "variants",
	Short: "output variants between ref and query from a SAM file",
	Long:  `output variants between ref and query from a SAM file`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		err = sam.Variants(samFile, samReference, variantGenbankFile, variantOutfile, samThreads)
		// err = sam.Variants(samFile, samReference, variantGenbankFile, variantOutfile, toPairAlignSkipInsertions, samThreads)

		return err
	},
}
