package cmd

import (
	"errors"
	"strconv"
	"strings"

	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/closest"
	"github.com/virus-evolution/gofasta/pkg/gfio"
)

var closestThreads int
var closestQuery string
var closestTarget string
var closestOutfile string
var closestN int
var closestDist string
var closestMeasure string
var closestTable bool

func init() {
	rootCmd.AddCommand(closestCmd)

	closestCmd.Flags().IntVarP(&closestThreads, "threads", "t", 0, "Number of CPUs to use (Default: all available CPUs)")
	closestCmd.Flags().StringVarP(&closestQuery, "query", "", "", "Alignment of sequences to find neighbours for, in fasta format")
	closestCmd.Flags().StringVarP(&closestTarget, "target", "", "", "Alignment of sequences to search for neighbours in, in fasta format")
	closestCmd.Flags().StringVarP(&closestMeasure, "measure", "m", "raw", "which distance measure to use (raw, snp or tn93)")
	closestCmd.Flags().IntVarP(&closestN, "number", "n", 0, "(Optional) the closest n sequences to each query will be returned")
	closestCmd.Flags().StringVarP(&closestDist, "max-dist", "d", "", "(Optional) return all sequences less than or equal to this distance away")
	closestCmd.Flags().StringVarP(&closestOutfile, "outfile", "o", "stdout", "The output file to write")
	closestCmd.Flags().BoolVarP(&closestTable, "table", "", false, "write a long-form table of the output")

	closestCmd.Flags().SortFlags = false
}

var closestCmd = &cobra.Command{
	Use:   "closest",
	Short: "Find the closest sequence(s) to a query by genetic distance",
	Long: `Find the closest sequence(s) to a query by genetic distance

The usecase is imagined to be a small number of query sequences, whose nearest neighbours
by genetic distance need to be found amongst a large number of target sequences. The query 
alignment is read into memory and the target alignment is streamed from disk and iterated 
over once, so it can be arbitrarily large.

You can find the single closest neighbour like:

	gofasta closest -t 2 --query query.fasta --target target.fasta -o closest.csv

and the output will be a CSV format file with the headers query, closest, distance, SNPs.

Or you can find the nearest n neighbours:

	gofasta closest -t 2 -n 1000 --query query.fasta --target target.fasta -o closest.n1000.csv

or you can find all the neighbours under a -d distance away:

	gofasta closest -t 2 --measure snp -d 5 --query query.fasta --target target.fasta -o closest.d5snps.csv

and the default output will be a CSV format file with the headers query, closest.  The 'closest' column
is a ";"-delimited list of neighbours, closest first. Ties for distance are broken by genome completeness.

You can combine -d with -n to find the nearest n neighbours less than or equal to a distance, d.

Possible measures of distance are raw number of nucleotide changes per site (the default, raw), raw number
of nucleotide changes in total (snp), or Tamura and Nei's 1993 evolutionary distance (tn93).

Use --table in combination with the -n and/or -d flags to write a long-form output including the distance
between every pair.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		queryIn, err := gfio.OpenIn(*cmd.Flag("query"))
		if err != nil {
			return err
		}
		defer queryIn.Close()

		targetIn, err := gfio.OpenIn(*cmd.Flag("target"))
		if err != nil {
			return err
		}
		defer targetIn.Close()

		var measure string
		switch strings.ToLower(closestMeasure) {
		case "raw":
			measure = "raw"
		case "snp":
			measure = "snp"
		case "tn93":
			measure = "tn93"
		default:
			return errors.New("Couldn't tell which distance --measure / -m to use (choose one of \"raw\", \"snp\" or \"tn93\")")
		}

		dist := -1.0
		if closestDist != "" {
			dist, err = strconv.ParseFloat(closestDist, 64)
			if err != nil {
				return err
			}
		}

		closestOut, err := gfio.OpenOut(*cmd.Flag("outfile"))
		if err != nil {
			return err
		}
		defer closestOut.Close()

		if closestN > 0 || dist != -1.0 {
			err = closest.ClosestN(closestN, dist, queryIn, targetIn, measure, closestOut, closestTable, closestThreads)
		} else {
			err = closest.Closest(queryIn, targetIn, measure, closestOut, closestThreads)
		}

		return err
	},
}
