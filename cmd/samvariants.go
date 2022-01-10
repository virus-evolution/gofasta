package cmd

import (
	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/sam"
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
	Short: "Find mutations relative to a reference from an alignment in sam format",
	Long: `Find mutations relative to a reference from an alignment in sam format

Example usage:
	gofasta sam variants -s aligned.sam -r reference.fasta -g annotation.gb -o variants.csv

Mutations are annotated with ins (insertion), del (deletion), aa (amino acid change) or nuc (a nucleotide change that
isn't in a codon that is represented by an amino acid change). The formats are:

ins:2028:3 - a 3-base insertion immediately after (1-based) position 2028 in reference coordinates
del:11288:9 - a 9-base deletion whose first missing nucleotide is at (1-based) position 11288 in reference coordinates
aa:s:D614G - the amino acid at (1-based) residue 614 in the S gene is a D in the reference and a G in this sequence
nuc:C3037T - the nucleotide at (1-based) position 3037 is a C in the reference and a T in this sequence

If input sam and output csv files are not specified, the behaviour is to read the sam from stdin and write
the variants to stdout.`,

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

		err = sam.Variants(samIn, ref, genbank, out, samThreads)

		return err
	},
}
