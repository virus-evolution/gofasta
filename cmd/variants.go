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
var variantsThreads int
var variantsAggregate bool
var variantsThreshold float64
var variantsAppendSNP bool

func init() {
	rootCmd.AddCommand(variantsCmd)

	variantsCmd.Flags().StringVarP(&variantsMSA, "msa", "", "stdin", "Multiple sequence alignment in fasta format")
	variantsCmd.Flags().StringVarP(&variantsReference, "reference", "r", "", "The name of the reference sequence in the msa")
	variantsCmd.Flags().StringVarP(&variantsGenbankFile, "genbank", "", "", "Genbank format annotation")
	variantsCmd.Flags().StringVarP(&variantsOutfile, "outfile", "o", "stdout", "Name of the file of variants to write")
	variantsCmd.Flags().BoolVarP(&variantsAggregate, "aggregate", "", false, "Report the proportions of each change")
	variantsCmd.Flags().Float64VarP(&variantsThreshold, "threshold", "", 0.0, "If --aggregate, only report changes with a freq greater than or equal to this value")
	variantsCmd.Flags().BoolVarP(&variantsAppendSNP, "append-snps", "", false, "Report the codon's SNPs in parenthesis after each amino acid mutation")
	variantsCmd.Flags().IntVarP(&variantsThreads, "threads", "t", 1, "Number of threads to use")

	variantsCmd.Flags().Lookup("aggregate").NoOptDefVal = "true"
	variantsCmd.Flags().Lookup("append-snps").NoOptDefVal = "true"

	variantsCmd.Flags().SortFlags = false
}

var variantsCmd = &cobra.Command{
	Use:   "variants",
	Short: "Find mutations relative to a reference from a multiple sequence alignment in fasta format",
	Long: `Find mutations relative to a reference from a multiple sequence alignment in fasta format

Example usage:

	./gofasta variants --msa alignment.fasta --genbank MN908947.gb --reference MN908947.3 > variants.csv

If input --msa and output csv files are not specified, the behaviour is to read the alignment from stdin and write
the variants to stdout.

If you are reading the --msa from stdin, you must either:
	1) provide a value to the --reference argument and this sequence name must be first in the alignment
							OR
	2) NOT provide a value to the --reference argument, in which case the genbank SOURCE will be used BUT 
	   the --msa must then be in reference coordinates (i.e. there can be no insertions relative to the reference 
	   in the alignment)

Mutations are annotated with ins (insertion), del (deletion), aa (amino acid change) or nuc (a nucleotide change that
isn't in a codon that is represented by an amino acid change). The formats are:

ins:2028:3 - a 3-base insertion immediately after (1-based) position 2028 in reference coordinates
del:11288:9 - a 9-base deletion whose first missing nucleotide is at (1-based) position 11288 in reference coordinates
aa:s:D614G - the amino acid at (1-based) residue 614 in the S gene is a D in the reference and a G in this sequence
nuc:C3037T - the nucleotide at (1-based) position 3037 is a C in the reference and a T in this sequence

Frame-shifting mutations in coding sequence are reported as indels but are ignored for subsequent amino-acids in the alignment.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		msa, err := gfio.OpenIn(*cmd.Flag("msa"))
		if err != nil {
			return err
		}
		defer msa.Close()

		stdin := false
		if variantsMSA == "stdin" {
			stdin = true
		}

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

		err = variants.Variants(msa, stdin, variantsReference, genbank, out, variantsAggregate, variantsThreshold, variantsAppendSNP, variantsThreads)

		return
	},
}
