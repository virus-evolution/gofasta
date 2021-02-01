package cmd

import (
	"github.com/spf13/cobra"

	"github.com/cov-ert/gofasta/pkg/snps"
)

var snpsReference string
var snpsQuery string
var snpsOutfile string

func init() {
	rootCmd.AddCommand(snpCmd)

	snpCmd.Flags().StringVarP(&snpsReference, "reference", "r", "", "Reference sequence, in fasta format")
	snpCmd.Flags().StringVarP(&snpsQuery, "query", "q", "stdin", "Alignment of sequences to find snps in, in fasta format")
	snpCmd.Flags().StringVarP(&snpsOutfile, "outfile", "o", "stdout", "Output to write")
}

var snpCmd = &cobra.Command{
	Use:   "snps",
	Short: "Find snps relative to a reference",
	Long:  `Find snps relative to a reference.

Example usage:
	gofasta snps -r reference.fasta -q alignment.fasta -o snps.csv

The output is a csv-format file with one line per query sequence, and two columns:
'query' and 'SNPs', the second of which is a "|"-delimited list of snps in that query.

If input and output files are not specified, the behaviour is to read the query alignment
from stdin and write the snps file to stdout, e.g. you could do this:
	cat alignment.fasta | gofasta snps -r reference.fasta > snps.csv`,

	RunE: func(cmd *cobra.Command, args []string) (err error) {

		err = snps.SNPs(snpsReference, snpsQuery, snpsOutfile)

		return
	},
}
