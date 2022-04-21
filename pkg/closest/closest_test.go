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
WTGATG
>Target3
WTTTTC
>Target4
ATGATG
>Target5
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

	err := Closest(query, target, "snp", out, 2)
	if err != nil {
		t.Error(err)
	}

	// fmt.Println(string(out.Bytes()))

	if string(out.Bytes()) != `query,closest,distance,SNPs
Query1,Target4,0,
Query2,Target1,0,
Query3,Target5,1,6GC
` {
		t.Errorf("problem in closest test")
	}
}
