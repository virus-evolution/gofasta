package cmd

import (
	"github.com/spf13/cobra"

	"github.com/cov-ert/gofasta/pkg/sam"
)

var tofastaOutfile string
var tofastaTrim bool
var tofastaPad bool
var tofastaTrimStart int
var tofastaTrimEnd int

func init() {
	samCmd.AddCommand(tofasta)

	tofasta.Flags().StringVarP(&tofastaOutfile, "fasta-out", "o", "stdout", "fasta file to write")
	tofasta.Flags().BoolVarP(&tofastaTrim, "trim", "", false, "trim the alignment")
	tofasta.Flags().BoolVarP(&tofastaPad, "pad", "", false, "if trim, pad the alignment with Ns")
	tofasta.Flags().IntVarP(&tofastaTrimStart, "trimstart", "", -1, "start coordinate for trimming")
	tofasta.Flags().IntVarP(&tofastaTrimEnd, "trimend", "", -1, "end coordinate for trimming")

	tofasta.Flags().SortFlags = false
}

var tofasta = &cobra.Command{
	Use:   "tofasta",
	Short: "convert a SAM file to an alignment in fasta format",
	Long:  `convert a SAM file to an alignment in fasta format`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		err = sam.ToFasta(samInfile, tofastaOutfile, tofastaTrim, tofastaPad, tofastaTrimStart, tofastaTrimEnd, samThreads)

		return
	},
}
