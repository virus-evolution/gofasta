package cmd

import (
	"github.com/spf13/cobra"

	"github.com/cov-ert/gofasta/pkg/gfio"
	"github.com/cov-ert/gofasta/pkg/updown"
)

var UDListReference string
var UDListQuery string
var UDListOutfile string

func init() {
	updownCmd.AddCommand(updownListCmd)

	updownListCmd.Flags().StringVarP(&UDListQuery, "query", "q", "stdin", "Alignment of sequences to parse, in fasta format")
	updownListCmd.Flags().StringVarP(&UDListOutfile, "outfile", "o", "stdout", "Output to write")

	updownListCmd.Flags().SortFlags = false
}

var updownListCmd = &cobra.Command{
	Use:   "list",
	Short: "Generate input CSV files for gofasta updown topranking",
	Long: `Generate input CSV files for gofasta updown topranking

Example usage:

	gofasta updown list -r reference.fasta -q alignment.fasta -o mutationlist.csv

Non-ATGC nucleotides are not recommended in the --reference, and --reference and --query must
be aligned to the same thing.

--outfile is a CSV-format file with the columns: query,SNPs,ambiguities,SNPcount,ambcount. There is one row
for each sequence in --query. SNPs is a "|"-delimited list of SNPs relative to --reference. ambiguities is
a "|"-delimited list of ranges (1-based, inclusive) of tracts of ambiguities (anything that isn't ATGC).
`,

	RunE: func(cmd *cobra.Command, args []string) (err error) {

		ref, err := gfio.OpenIn(*cmd.Flag("reference"))
		if err != nil {
			return err
		}
		defer ref.Close()

		query, err := gfio.OpenIn(*cmd.Flag("query"))
		if err != nil {
			return err
		}
		defer query.Close()

		out, err := gfio.OpenOut(*cmd.Flag("outfile"))
		if err != nil {
			return err
		}
		defer out.Close()

		err = updown.List(ref, query, out)

		return
	},
}
