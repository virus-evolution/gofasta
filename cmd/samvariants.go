package cmd

import (
	"github.com/spf13/cobra"

	"github.com/cov-ert/gofasta/pkg/gfio"
	"github.com/cov-ert/gofasta/pkg/sam"
)

var variantGenbankFile string
var variantOutfile string

func init() {
	samCmd.AddCommand(variantCmd)

	variantCmd.Flags().StringVarP(&variantGenbankFile, "genbank", "g", "", "Genbank format annotation of a sequence in the same coordinates as the alignment")
	variantCmd.Flags().StringVarP(&variantOutfile, "outfile", "o", "stdout", "Where to write the variants")

	variantCmd.Flags().SortFlags = false
}

var variantCmd = &cobra.Command{
	Use:   "variants",
	Short: "Call variants between ref and query from a SAM file",
	Long: `Call variants between a reference sequence and query sequences aligned in sam format

Example usage:
	gofasta sam variants -s aligned.sam -r reference.fasta -g annotation.gb -o variants.csv

The output is a csv-format file with one line per query sequence, and two columns: 'query' and
'variants', the second of which is a "|"-delimited list of amino acid changes and synonymous SNPs
in that query relative to the reference sequence specified using --reference/-r.

If input sam and output csv files are not specified, the behaviour is to read the sam from stdin and write
the variants to stdout.`,

	RunE: func(cmd *cobra.Command, args []string) (err error) {

		samIn, err := gfio.OpenIn(samFile)
		if err != nil {
			return err
		}
		defer samIn.Close()

		ref, err := gfio.OpenIn(samReference)
		if err != nil {
			return err
		}
		defer ref.Close()

		genbank, err := gfio.OpenIn(variantGenbankFile)
		if err != nil {
			return err
		}
		defer genbank.Close()

		out, err := gfio.OpenOut(variantOutfile)
		if err != nil {
			return err
		}
		defer out.Close()

		err = sam.Variants(samIn, ref, genbank, out, samThreads)

		return err
	},
}
