package cmd

import (
	"github.com/spf13/cobra"

	"github.com/cov-ert/gofasta/pkg/closest"
)

var closestThreads int
var closestQuery string
var closestTarget string
var closestOutfile string
var closestN int

func init() {
	rootCmd.AddCommand(closestCmd)

	closestCmd.Flags().IntVarP(&closestThreads, "threads", "t", 0, "Number of CPUs to use (Default: all available CPUs)")
	closestCmd.Flags().StringVarP(&closestQuery, "query", "", "", "Alignment of sequences to find neighbours for, in fasta format")
	closestCmd.Flags().StringVarP(&closestTarget, "target", "", "", "Alignment of sequences to search for neighbours in, in fasta format")
	closestCmd.Flags().IntVarP(&closestN, "number", "n", 0, "(Optional) the closest n sequences to each query will be returned")
	closestCmd.Flags().StringVarP(&closestOutfile, "outfile", "o", "stdout", "The output file to write")
}

var closestCmd = &cobra.Command{
	Use:   "closest",
	Short: "Find the closest sequence(s) to a query by raw distance",
	Long:  `Find the closest sequence(s) to a query by raw distance

The usecase is imagined to be a small number of query sequences, whose nearest neighbours
need to be found amongst a large number of target sequences. The query alignment is read
into memory and the target alignment is read from disk and iterated over once.

Closest neighbours are those with the lowest raw distance per site to the query sequence,
and ties for this score are broken by how unambiguous the target genomes are.

You can find the single closest neighbour like:

	gofasta closest -t 2 --query query.fasta --target target.fasta -o closest.csv

and the output will be a CSV format file with the headers query, closest, SNPdistance, SNPs.

Or you can find the nearest n neighbours:

	gofasta closest -t 2 -n 1000 --query query.fasta --target target.fasta -o closest.n1000.csv

and the output will be a CSV format file with just the headers query, closest. The 'closest' column
is a ";"-delimited list of neighbours, closest first.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		if closestN > 0 {
			err = closest.ClosestN(closestN, closestQuery, closestTarget, closestOutfile, closestThreads)
		} else {
			err = closest.Closest(closestQuery, closestTarget, closestOutfile, closestThreads)
		}

		return err
	},
}
