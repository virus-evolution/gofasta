package cmd

import (
	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/snps"
)

var snpsReference string
var snpsQuery string
var snpsOutfile string
var hardGaps bool
var aggregate bool
var thresh float64

func init() {
	rootCmd.AddCommand(snpCmd)

	snpCmd.Flags().StringVarP(&snpsReference, "reference", "r", "", "Reference sequence, in fasta format")
	snpCmd.Flags().StringVarP(&snpsQuery, "query", "q", "stdin", "Alignment of sequences to find snps in, in fasta format")
	snpCmd.Flags().StringVarP(&snpsOutfile, "outfile", "o", "stdout", "Output to write")
	snpCmd.Flags().BoolVarP(&hardGaps, "hard-gaps", "", false, "Don't treat alignment gaps as missing data")
	snpCmd.Flags().BoolVarP(&aggregate, "aggregate", "", false, "Report the proportions of each change")
	snpCmd.Flags().Float64VarP(&thresh, "threshold", "", 0.0, "If --aggregate, only report snps with a freq greater than or equal to this value")

	snpCmd.Flags().Lookup("hard-gaps").NoOptDefVal = "true"
	snpCmd.Flags().Lookup("aggregate").NoOptDefVal = "true"

	snpCmd.Flags().SortFlags = false
}

var snpCmd = &cobra.Command{
	Use:   "snps",
	Short: "Find snps relative to a reference",
	Long: `Find snps relative to a reference.

Example usage:
	gofasta snps -r reference.fasta -q alignment.fasta -o snps.csv

reference.fasta and alignment.fasta must be the same length.

With the default settings the output is a csv-format file with one line per query sequence, and two columns:
'query' and 'SNPs', the second of which is a "|"-delimited list of snps in that query.

If you set --aggregate (and optionally a --threshold) it will return the SNPs present in the entire sample
(whose frequency is equal to/above --threshold) and their frequencies.

Setting --hard-gaps treats alignment gaps as different from {ATGC}.

If query and outfile are not specified, the behaviour is to read the query alignment
from stdin and write the snps file to stdout, e.g. you could do this:
	cat alignment.fasta | gofasta snps -r reference.fasta > snps.csv`,

	RunE: func(cmd *cobra.Command, args []string) (err error) {

		query, err := gfio.OpenIn(*cmd.Flag("query"))
		if err != nil {
			return err
		}
		defer query.Close()

		ref, err := gfio.OpenIn(*cmd.Flag("reference"))
		if err != nil {
			return err
		}
		defer ref.Close()

		out, err := gfio.OpenIn(*cmd.Flag("outfile"))
		if err != nil {
			return err
		}
		defer out.Close()

		err = snps.SNPs(ref, query, hardGaps, aggregate, thresh, out)

		return
	},
}
