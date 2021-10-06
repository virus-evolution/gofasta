package updown

import (
	"bytes"
	"testing"
)

func TestTopRanking(t *testing.T) {
	refData := []byte(`>ref
ATGATG
`)
	queryData := []byte(
		`>Query1
ATTATT
`)

	targetData := []byte(`>TargetUp1
ATGATG
>TargetSame1
ATTATT
>TargetDown1
ATTACT
>TargetUp2
ATGATT
>TargetSide1
ATGCTT
`)
	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)
	target := bytes.NewReader(targetData)
	out := new(bytes.Buffer)

	qtype := "fasta"
	ttype := "fasta"
	ignoreArray := make([]string, 0)
	TRsizetotal := 5
	TRsizeup := 0
	TRsizedown := 0
	TRsizeside := 0
	TRsizesame := 0
	TRdistall := 0
	TRdistup := 0
	TRdistdown := 0
	TRdistside := 0
	TRthresholdpair := float32(0.1)
	TRthresholdtarget := 10000
	TRnofill := false
	TRdistpush := 0

	err := TopRanking(query, target, ref, out,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != `query,closestsame,closestup,closestdown,closestside
Query1,TargetSame1,TargetUp1;TargetUp2,TargetDown1,TargetSide1
` {
		t.Errorf("problem in TestTopRanking()")
	}
}
