package cmd

import (
	"github.com/spf13/cobra"

	"github.com/cov-ert/gofasta/pkg/gfio"
	"github.com/cov-ert/gofasta/pkg/variants"
)

var variantsMSA string
var variantsReference string
var variantsGenbankFile string
var variantsOutfile string

func init() {
	rootCmd.AddCommand(variantsCmd)

	variantsCmd.Flags().StringVarP(&variantsMSA, "msa", "", "stdin", "multiple sequence alignment in fasta format")
	variantsCmd.Flags().StringVarP(&variantsReference, "reference", "r", "", "the name of the reference sequence in the msa")
	variantsCmd.Flags().StringVarP(&variantsGenbankFile, "genbank", "", "", "genbank format annotation")
	variantsCmd.Flags().StringVarP(&variantsOutfile, "outfile", "o", "stdout", "name of the file of variants to write")

	variantsCmd.Flags().SortFlags = false
}

var variantsCmd = &cobra.Command{
	Use:   "variants",
	Short: "find variants relative to a reference in a multiple sequence alignment",
	Long: `find variants relative to a reference in a multiple sequence alignment

Example usage:

./gofasta variants --msa alignment.fasta --genbank MN908947.gb --reference MN908947.3 > variants.csv
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		msa, err := gfio.OpenIn(variantsMSA)
		if err != nil {
			return err
		}
		defer msa.Close()

		genbank, err := gfio.OpenIn(variantsGenbankFile)
		if err != nil {
			return err
		}
		defer genbank.Close()

		out, err := gfio.OpenOut(variantsOutfile)
		if err != nil {
			return err
		}
		defer out.Close()

		err = variants.Variants(msa, variantsReference, genbank, out)

		return
	},
}
