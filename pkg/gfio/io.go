package gfio

import (
	"os"
)

func OpenIn(inFile string) (*os.File, error) {
	var err error
	var f *os.File

	if inFile != "stdin" {
		f, err = os.Open(inFile)
		if err != nil {
			return f, err
		}
	} else {
		f = os.Stdin
	}

	return f, nil
}

func OpenOut(outFile string) (*os.File, error) {
	var err error
	var f *os.File

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
