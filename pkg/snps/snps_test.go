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

	err := SNPs(ref, query, false, false, 0.0, out)
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

func TestSNPsHardGaps(t *testing.T) {
	refData := []byte(`>ref
ATGATG
`)
	queryData := []byte(
		`>Query1
--GATG
>Query2
ATGATC
>Query3
ATTTTW
`)

	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)

	out := new(bytes.Buffer)

	err := SNPs(ref, query, true, false, 0.0, out)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != `query,SNPs
Query1,A1-|T2-
Query2,G6C
Query3,G3T|A4T|G6W
` {
		t.Errorf("problem in TestSNPsHardGaps()")
	}
}

func TestSNPsAggregate(t *testing.T) {
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
>Query4
ATTTTG
`)

	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)

	out := new(bytes.Buffer)

	err := SNPs(ref, query, false, true, 0.0, out)
	if err != nil {
		t.Error(err)
	}

	// strconv.FormatFloat(propMap[snp]/counter, 'f', 4, 64) + "\n")

	if string(out.Bytes()) != `SNP,frequency
G3T,0.500000000
A4T,0.500000000
G6C,0.250000000
G6W,0.250000000
` {
		t.Errorf("problem in TestSNPsAggregate()")
	}
}

func TestSNPsAggregateThresh(t *testing.T) {
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
>Query4
ATTTTG
`)

	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)

	out := new(bytes.Buffer)

	err := SNPs(ref, query, false, true, 0.26, out)
	if err != nil {
		t.Error(err)
	}

	// strconv.FormatFloat(propMap[snp]/counter, 'f', 4, 64) + "\n")

	if string(out.Bytes()) != `SNP,frequency
G3T,0.500000000
A4T,0.500000000
` {
		t.Errorf("problem in TestSNPsAggregateThresh()")
	}
}
