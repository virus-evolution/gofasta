package gfio

import (
	"errors"
	"os"
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
	if err.Error() != errors.New("open -r / --reference not/a/file.whatever: no such file or directory").Error() {
		t.Error(err)
	}
}

func TestCloseOutReportsCloseError(t *testing.T) {
	f, err := os.CreateTemp(t.TempDir(), "x")
	if err != nil {
		t.Fatal(err)
	}
	f.Close() // already closed, so the next Close fails

	var got error
	CloseOut(f, &got)
	if got == nil {
		t.Errorf("expected the close error")
	}

	first := errors.New("first")
	got = first
	CloseOut(f, &got)
	if got != first {
		t.Errorf("an existing error must not be overwritten")
	}
}
