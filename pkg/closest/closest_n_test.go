package closest

import (
	"bytes"
	"testing"
)

func TestClosestN(t *testing.T) {
	targetData := []byte(
		`>Target1
ATGATC
>Target2
ATGATG
>Target3
ATTAGG
>Target4
ATTATG
>Target5
ATTATT
`)

	queryData := []byte(
		`>Query1
ATGATG
>Query2
ATGATC
>Query3
ATTATT
`)

	target := bytes.NewReader(targetData)

	query := bytes.NewReader(queryData)

	buf := new(bytes.Buffer)

	err := ClosestN(2, -1.0, query, target, "raw", buf, false, 2)
	if err != nil {
		t.Error(err)
	}

	if string(buf.Bytes()) != `query,closest
Query1,Target2;Target1
Query2,Target1;Target2
Query3,Target5;Target4
` {
		t.Errorf("problem in closest test")
	}
}
