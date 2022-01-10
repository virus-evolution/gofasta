package cmd

import (
	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/sam"
)

var samVariantsGenbankFile string
var samVariantsOutfile string
var samVariantsAggregate bool
var samVariantsThreshold float64
var samVariantsAppendSNP bool

func init() {
	samCmd.AddCommand(samVariantsCmd)

	samVariantsCmd.Flags().StringVarP(&samVariantsGenbankFile, "genbank", "g", "", "Genbank format annotation of a sequence in the same coordinates as the alignment")
	samVariantsCmd.Flags().StringVarP(&samVariantsOutfile, "outfile", "o", "stdout", "Where to write the variants")
	samVariantsCmd.Flags().BoolVarP(&samVariantsAggregate, "aggregate", "", false, "Report the proportions of each change")
	samVariantsCmd.Flags().Float64VarP(&samVariantsThreshold, "threshold", "", 0.0, "If --aggregate, only report changes with a freq greater than or equal to this value")
	samVariantsCmd.Flags().BoolVarP(&samVariantsAppendSNP, "append-snps", "", false, "Report the codon's SNPs in parenthesis after each amino acid mutation")

	samVariantsCmd.Flags().Lookup("aggregate").NoOptDefVal = "true"
	samVariantsCmd.Flags().Lookup("append-snps").NoOptDefVal = "true"

	samVariantsCmd.Flags().SortFlags = false
}

var samVariantsCmd = &cobra.Command{
	Use:   "variants",
	Short: "Find mutations relative to a reference from an alignment in sam format",
	Long: `Find mutations relative to a reference from an alignment in sam format

Example usage:
	gofasta sam variants -s aligned.sam -r reference.fasta -g annotation.gb -o variants.csv

If input sam and output csv files are not specified, the behaviour is to read the sam from stdin and write
the variants to stdout.

Mutations are annotated with ins (insertion), del (deletion), aa (amino acid change) or nuc (a nucleotide change that
isn't in a codon that is represented by an amino acid change). The formats are:

ins:2028:3 - a 3-base insertion immediately after (1-based) position 2028 in reference coordinates
del:11288:9 - a 9-base deletion whose first missing nucleotide is at (1-based) position 11288 in reference coordinates
aa:s:D614G - the amino acid at (1-based) residue 614 in the S gene is a D in the reference and a G in this sequence
nuc:C3037T - the nucleotide at (1-based) position 3037 is a C in the reference and a T in this sequence

Frame-shifting mutations in coding sequence are reported as indels but are ignored for subsequent amino-acids in the alignment.
`,

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

		err = sam.Variants(samIn, ref, genbank, out, samVariantsAggregate, samVariantsThreshold, samVariantsAppendSNP, samThreads)

		return err
	},
}
