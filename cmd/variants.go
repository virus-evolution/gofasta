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
	Short: "find mutations relative to a reference in a multiple sequence alignment",
	Long: `find mutations relative to a reference in a multiple sequence alignment

Example usage:

	./gofasta variants --msa alignment.fasta --genbank MN908947.gb --reference MN908947.3 > variants.csv

Mutations are annotated with ins (insertion), del (deletion), aa (amino acid change) or nuc (a nucleotide change that
isn't in a codon that is represented by an amino acid change). The formats are:

ins:2028:3 - a 3-base insertion immediately after (1-based) position 2028 in reference coordinates
del:11288:9 - a 9-base deletion whose first missing nucleotide is at (1-based) position 11288 in reference coordinates
aa:s:D614G - the amino acid at (1-based) residue 614 in the S gene is a D in the reference and a G in this sequence
nuc:C3037T - the nucleotide at (1-based) position 3037 is a C in the reference and a T in this sequence
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		msa, err := gfio.OpenIn(*cmd.Flag("msa"))
		if err != nil {
			return err
		}
		defer msa.Close()

		genbank, err := gfio.OpenIn(*cmd.Flag("genbank"))
		if err != nil {
			return err
		}
		defer genbank.Close()

		out, err := gfio.OpenOut(*cmd.Flag("outfile"))
		if err != nil {
			return err
		}
		defer out.Close()

		err = variants.Variants(msa, variantsReference, genbank, out)

		return
	},
}
