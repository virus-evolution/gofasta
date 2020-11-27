package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

var (
	rootCmd = &cobra.Command{
		Use:   "gofasta",
		Short: "some functions for working with alignments",
		Long:  `some functions for working with alignments`,
	}
)

// Execute executes the root command.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}

// func main() {
//   closestCmd := flag.NewFlagSet("closest", flag.ExitOnError)
//
//   closestThreads := closestCmd.Int("t", 1, "number of threads")
//   closestQuery := closestCmd.String("query", "query.fasta", "query alignment")
//   closestTarget := closestCmd.String("target", "target.fasta", "target alignment")
//   closestOutfile := closestCmd.String("outfile", "snps.csv", "snps file to write")
//
//   if len(os.Args) < 2 {
//     fmt.Println("expected 'closest' subcommand(s)")
//     os.Exit(1)
//   }
//
//   if os.Args[1] != "closest" {
//     fmt.Println("expected 'closest' subcommand(s)")
//     os.Exit(1)
//   }
//
//   switch os.Args[1] {
//   case "closest":
//     closestCmd.Parse(os.Args[2:])
//     fmt.Println("subcommand 'closest'")
//     fmt.Println("  query:", *closestQuery)
//     fmt.Println("  target:", *closestTarget)
//     fmt.Println("  outfile:", *closestOutfile)
//     fmt.Println("  threads:", *closestThreads)
//
//     closest.Closest(*closestQuery, *closestTarget, *closestOutfile, *closestThreads)
//   }
//
// }
