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
`)

	queryData := []byte(
		`>Query1
ATGATG
>Query2
ATGATC
`)

	target := bytes.NewReader(targetData)

	query := bytes.NewReader(queryData)

	buf := new(bytes.Buffer)

	err := ClosestN(2, query, target, buf, 2)
	if err != nil {
		t.Error(err)
	}

	if string(buf.Bytes()) != `query,closest
Query1,Target2;Target1
Query2,Target1;Target2
` {
		t.Errorf("problem in closest test")
	}
}
