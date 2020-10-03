package cmd

import (
	"github.com/spf13/cobra"
)

var samThreads int
var samInfile string

func init() {
	rootCmd.AddCommand(samCmd)

	tofasta.PersistentFlags().IntVarP(&samThreads, "threads", "t", 1, "Number of threads to use")
	tofasta.PersistentFlags().StringVarP(&samInfile, "samfile", "s", "", "samfile to read")
}

var samCmd = &cobra.Command{
	Use:   "sam",
	Short: "Do things with SAM files",
	Long:  `Do things with SAM files`,
	Args:  cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {

		return nil
	},
}
