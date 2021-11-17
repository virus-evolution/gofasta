package updown

import (
	"bytes"
	"testing"
)

func TestList(t *testing.T) {
	refData := []byte(`>ref
ATGATG
`)
	queryData := []byte(
		`>Target1
ATGATG
>Target2
ATGATC
>Target3
ATTTTW
`)

	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)

	out := new(bytes.Buffer)

	err := List(ref, query, out)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != `query,SNPs,ambiguities,SNPcount,ambcount
Target1,,,0,0
Target2,G6C,,1,0
Target3,G3T|A4T,6,2,1
` {
		t.Errorf("problem in TestList()")
	}
}
