package cmd

import (
	"github.com/spf13/cobra"

	"github.com/cov-ert/gofasta/pkg/closest"
)

var closestThreads int
var closestQuery string
var closestTarget string
var closestOutfile string

func init() {
	rootCmd.AddCommand(closestCmd)

	closestCmd.Flags().IntVarP(&closestThreads, "threads", "t", 1, "Number of threads to use")
	closestCmd.Flags().StringVarP(&closestQuery, "query", "", "", "Alignment of sequences to find matches for, in fasta format")
	closestCmd.Flags().StringVarP(&closestTarget, "target", "", "", "Alignment of sequences to search for matches in, in fasta format")
	closestCmd.Flags().StringVarP(&closestOutfile, "outfile", "o", "snps.csv", "The output file to write")
}

var closestCmd = &cobra.Command{
	Use:   "closest",
	Short: "Find the closest sequence to a query by SNP-distance",
	Long:  `Find the closest sequence to a query by SNP-distance`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		err = closest.Closest(closestQuery, closestTarget, closestOutfile, closestThreads)

		return
	},
}
