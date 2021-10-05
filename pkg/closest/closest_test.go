package closest

import (
	"bytes"
	"testing"
)

func TestClosest(t *testing.T) {
	targetData := []byte(
		`>Target1
ATGATC
>Target2
ATGATG
>Target3
ATTTTC
`)

	queryData := []byte(
		`>Query1
ATGATG
>Query2
ATGATC
>Query3
ATTTTG
`)

	target := bytes.NewReader(targetData)

	query := bytes.NewReader(queryData)

	out := new(bytes.Buffer)

	err := Closest(query, target, out, 2)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != `query,closest,SNPdistance,SNPs
Query1,Target2,0,
Query2,Target1,0,
Query3,Target3,1,6GC
` {
		t.Errorf("problem in closest test")
	}
}
