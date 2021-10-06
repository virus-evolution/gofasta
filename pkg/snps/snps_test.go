package snps

import (
	"bytes"
	"testing"
)

func TestSNPs(t *testing.T) {
	refData := []byte(`>ref
ATGATG
`)
	queryData := []byte(
		`>Query1
ATGATG
>Query2
ATGATC
>Query3
ATTTTW
`)

	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)

	out := new(bytes.Buffer)

	err := SNPs(ref, query, out)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != `query,SNPs
Query1,
Query2,G6C
Query3,G3T|A4T|G6W
` {
		t.Errorf("problem in TestSNPs()")
	}
}
