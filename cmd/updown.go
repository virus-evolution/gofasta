package cmd

import (
	"github.com/spf13/cobra"
)

var udReference string


func init() {
	rootCmd.AddCommand(updownCmd)

	updownCmd.PersistentFlags().StringVarP(&udReference, "reference", "r", "", "Reference sequence, in fasta format")
}

var updownCmd = &cobra.Command{
	Use:   "updown",
	Short: "get pseudo-tree-aware catchments for query sequences from alignments",
	Long:  `get pseudo-tree-aware catchments for query sequences from alignments`,
	Args:  cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {

		return nil
	},
}
