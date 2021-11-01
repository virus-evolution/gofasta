package gfio

import (
	"errors"
	"io/fs"
	"os"

	"github.com/spf13/pflag"
)

// func (e *fs.PathError) Error() string { return e.Op + " " + e.Path + ": " + e.Err.Error() }

func parseInErr(err error, flagString string) error {
	switch x := err.(type) {
	case *fs.PathError:
		return errors.New(x.Op + " " + flagString + " " + x.Path + ": " + x.Err.Error())
	default:
		return err
	}
}

func OpenIn(flag pflag.Flag) (*os.File, error) {
	var err error
	var f *os.File

	inFile := flag.Value.String()
	flagString := "--" + flag.Name + " / -" + flag.Shorthand

	if inFile != "stdin" {
		if f, err = os.Open(inFile); err != nil {
			err = parseInErr(err, flagString)
			return f, err
		}
	} else {
		f = os.Stdin
	}

	return f, nil
}

func OpenOut(flag pflag.Flag) (*os.File, error) {
	var err error
	var f *os.File

	outFile := flag.Value.String()

	if outFile != "stdout" {
		f, err = os.Create(outFile)
		if err != nil {
			return f, err
		}
	} else {
		f = os.Stdout
	}

	return f, nil
}
