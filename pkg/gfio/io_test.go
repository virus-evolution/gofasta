package gfio

import (
	"errors"
	"testing"

	"github.com/spf13/cobra"
)

func TestOpenIn(t *testing.T) {

	var (
		Cmd = &cobra.Command{
			Use:     "test",
			Short:   "test",
			Long:    `test`,
			Version: "1.0",
		}
	)

	var reference string
	Cmd.PersistentFlags().StringVarP(&reference, "reference", "r", "", "Reference fasta file")
	Cmd.PersistentFlags().Set("reference", "not/a/file.whatever")

	_, err := OpenIn(*Cmd.Flag("reference"))
	if err.Error() != errors.New("open"+" "+"--reference / -r"+" "+"not/a/file.whatever"+": "+"no such file or directory").Error() {
		t.Error(err)
	}
}
