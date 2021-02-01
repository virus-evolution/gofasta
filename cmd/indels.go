package cmd

import (
	"github.com/spf13/cobra"

	"github.com/cov-ert/gofasta/pkg/sam"
)

var indelsInsOut string
var indelsDelOut string
var indelsThreshold int

func init() {
	samCmd.AddCommand(indelCmd)

	indelCmd.Flags().StringVarP(&indelsInsOut, "insertions-out", "", "insertions.txt", "Where to write the insertions")
	indelCmd.Flags().StringVarP(&indelsDelOut, "deletions-out", "", "deletions.txt", "Where to write the deletions")
	indelCmd.Flags().IntVarP(&indelsThreshold, "threshold", "", 2, "Minimum count for an indel to be included in the output")

	indelCmd.Flags().SortFlags = false
}

var indelCmd = &cobra.Command{
	Use:   "indels",
	Short: "Parse a SAM file for raw indel information",
	Long:  `Parse a SAM file for raw indel information

Parse a sam file for raw insertion and deletion information stored in the CIGAR. No attempt
is made to consolidate indels within one query sequence's sam lines, so there may be some conflict.

Two files (default: insertions.txt and deletions.txt) are written, which contain a list of insertions
or deletions represented in more than a threshold (default: 2) number of queries.

The format of insertions.txt is a three-column, tab-separated file with the headers: ref_start	insertion	samples
The format of deletions.txt is a three-column, tab-separated file with the headers: ref_start	length	samples

the 'samples' column is a "|"-separated list of the queries with the insertion/deletion described by the first two columns.

Example usage:
	gofasta sam indels -s aligned.sam --threshold 2 --insertions-out insertions.txt --deletions-out deletions.txt
`,

	RunE: func(cmd *cobra.Command, args []string) (err error) {

		err = sam.Indels(samFile, indelsInsOut, indelsDelOut, indelsThreshold)

		return
	},
}
