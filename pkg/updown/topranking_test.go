package updown

import (
	"bytes"
	"fmt"
	"testing"
)

func TestTopRanking1(t *testing.T) {
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
CCCCCC
>TargetSide2
ATGCTT
`)
	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)
	target := bytes.NewReader(targetData)
	out := new(bytes.Buffer)
	table := false
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

	err := TopRanking(query, target, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	desiredResult := `query,closestsame,closestup,closestdown,closestside
Query1,TargetSame1,TargetUp2;TargetUp1,TargetDown1,TargetSide2
`
	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRanking1(fasta)")
		fmt.Println(string(out.Bytes()))
	}

	ref = bytes.NewReader(refData)
	query = bytes.NewReader(queryData)
	queryList := new(bytes.Buffer)
	err = List(ref, query, queryList)
	if err != nil {
		t.Error(err)
	}

	ref = bytes.NewReader(refData)
	target = bytes.NewReader(targetData)
	targetList := new(bytes.Buffer)
	err = List(ref, target, targetList)
	if err != nil {
		t.Error(err)
	}

	qtype = "csv"
	ttype = "csv"
	out = new(bytes.Buffer)

	err = TopRanking(queryList, targetList, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRanking1(csv)")
	}
}

func TestTopRanking2(t *testing.T) {
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
CCCCCC
>TargetSide2
ATGCTT
`)
	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)
	target := bytes.NewReader(targetData)
	out := new(bytes.Buffer)
	table := false
	qtype := "fasta"
	ttype := "fasta"
	ignoreArray := make([]string, 0)
	TRsizetotal := 0
	TRsizeup := 1
	TRsizedown := 1
	TRsizeside := 1
	TRsizesame := 1
	TRdistall := 0
	TRdistup := 0
	TRdistdown := 0
	TRdistside := 0
	TRthresholdpair := float32(0.1)
	TRthresholdtarget := 10000
	TRnofill := false
	TRdistpush := 0

	err := TopRanking(query, target, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	desiredResult := `query,closestsame,closestup,closestdown,closestside
Query1,TargetSame1,TargetUp2,TargetDown1,TargetSide2
`

	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRanking2()")
	}

	ref = bytes.NewReader(refData)
	query = bytes.NewReader(queryData)
	queryList := new(bytes.Buffer)
	err = List(ref, query, queryList)
	if err != nil {
		t.Error(err)
	}

	ref = bytes.NewReader(refData)
	target = bytes.NewReader(targetData)
	targetList := new(bytes.Buffer)
	err = List(ref, target, targetList)
	if err != nil {
		t.Error(err)
	}

	qtype = "csv"
	ttype = "csv"
	out = new(bytes.Buffer)

	err = TopRanking(queryList, targetList, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRanking2(csv)")
	}
}

func TestTopRanking3(t *testing.T) {
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
CCCCCC
>TargetSide2
ATGCTT
`)
	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)
	target := bytes.NewReader(targetData)
	out := new(bytes.Buffer)
	table := false
	qtype := "fasta"
	ttype := "fasta"
	ignoreArray := make([]string, 0)
	TRsizetotal := 0
	TRsizeup := 0
	TRsizedown := 0
	TRsizeside := 0
	TRsizesame := 0
	TRdistall := 1
	TRdistup := 0
	TRdistdown := 0
	TRdistside := 0
	TRthresholdpair := float32(0.1)
	TRthresholdtarget := 10000
	TRnofill := false
	TRdistpush := 2

	err := TopRanking(query, target, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	desiredResult := `query,closestsame,closestup,closestdown,closestside
Query1,TargetSame1,TargetUp2;TargetUp1,TargetDown1,TargetSide2;TargetSide1
`

	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRanking3(fasta)")
	}

	ref = bytes.NewReader(refData)
	query = bytes.NewReader(queryData)
	queryList := new(bytes.Buffer)
	err = List(ref, query, queryList)
	if err != nil {
		t.Error(err)
	}

	ref = bytes.NewReader(refData)
	target = bytes.NewReader(targetData)
	targetList := new(bytes.Buffer)
	err = List(ref, target, targetList)
	if err != nil {
		t.Error(err)
	}

	qtype = "csv"
	ttype = "csv"
	out = new(bytes.Buffer)

	err = TopRanking(queryList, targetList, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRanking3(csv)")
	}
}

func TestTopRanking4(t *testing.T) {
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
CCCCCC
>TargetSide2
ATNCTN
>TargetSide3
ATNATC
`)
	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)
	target := bytes.NewReader(targetData)
	out := new(bytes.Buffer)
	table := false
	qtype := "fasta"
	ttype := "fasta"
	ignoreArray := make([]string, 0)
	TRsizetotal := 0
	TRsizeup := 2
	TRsizedown := 2
	TRsizeside := 2
	TRsizesame := 2
	TRdistall := 0
	TRdistup := 0
	TRdistdown := 0
	TRdistside := 0
	TRthresholdpair := float32(0.5)
	TRthresholdtarget := 10000
	TRnofill := false
	TRdistpush := 0

	err := TopRanking(query, target, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	desiredResult := `query,closestsame,closestup,closestdown,closestside
Query1,TargetSame1,TargetUp2;TargetUp1,TargetDown1,TargetSide3;TargetSide1
`

	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRanking4(fasta)")
	}

	ref = bytes.NewReader(refData)
	query = bytes.NewReader(queryData)
	queryList := new(bytes.Buffer)
	err = List(ref, query, queryList)
	if err != nil {
		t.Error(err)
	}

	ref = bytes.NewReader(refData)
	target = bytes.NewReader(targetData)
	targetList := new(bytes.Buffer)
	err = List(ref, target, targetList)
	if err != nil {
		t.Error(err)
	}

	qtype = "csv"
	ttype = "csv"
	out = new(bytes.Buffer)

	err = TopRanking(queryList, targetList, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRanking4(csv)")
	}
}

func TestTopRanking5(t *testing.T) {
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
CCCCCC
>TargetSide2
ATNCTN
>TargetSide3
ATNATC
`)
	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)
	target := bytes.NewReader(targetData)
	out := new(bytes.Buffer)
	table := false
	qtype := "fasta"
	ttype := "fasta"
	ignoreArray := make([]string, 0)
	TRsizetotal := 100
	TRsizeup := 0
	TRsizedown := 0
	TRsizeside := 0
	TRsizesame := 0
	TRdistall := 0
	TRdistup := 0
	TRdistdown := 0
	TRdistside := 0
	TRthresholdpair := float32(1.0)
	TRthresholdtarget := 0
	TRnofill := false
	TRdistpush := 0

	err := TopRanking(query, target, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	desiredResult := `query,closestsame,closestup,closestdown,closestside
Query1,TargetSame1,TargetUp2;TargetUp1,TargetDown1,TargetSide1
`

	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRanking5(fasta)")
	}

	ref = bytes.NewReader(refData)
	query = bytes.NewReader(queryData)
	queryList := new(bytes.Buffer)
	err = List(ref, query, queryList)
	if err != nil {
		t.Error(err)
	}

	ref = bytes.NewReader(refData)
	target = bytes.NewReader(targetData)
	targetList := new(bytes.Buffer)
	err = List(ref, target, targetList)
	if err != nil {
		t.Error(err)
	}

	qtype = "csv"
	ttype = "csv"
	out = new(bytes.Buffer)

	err = TopRanking(queryList, targetList, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRanking5(csv)")
	}
}

func TestTopRankingDist2(t *testing.T) {
	refData := []byte(`>ref
ATGATG
`)
	queryData := []byte(
		`>Query1
ATTATT
`)

	targetData := []byte(`>TargetUp1_dist2
ATGATG
>TargetSame1_dist0
ATTATT
>TargetDown1_dist1
ATTACT
>TargetUp2_dist1
ATGATT
>TargetSide1_dist6
CCCCCC
>TargetSide2_2
ATGCTT
`)
	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)
	target := bytes.NewReader(targetData)
	out := new(bytes.Buffer)
	table := false
	qtype := "fasta"
	ttype := "fasta"
	ignoreArray := make([]string, 0)
	TRsizetotal := 0
	TRsizeup := 0
	TRsizedown := 0
	TRsizeside := 0
	TRsizesame := 0
	TRdistall := 0
	TRdistup := 1
	TRdistdown := 1
	TRdistside := 1
	TRthresholdpair := float32(0.1)
	TRthresholdtarget := 10000
	TRnofill := false
	TRdistpush := 0

	err := TopRanking(query, target, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	desiredResult := `query,closestsame,closestup,closestdown,closestside
Query1,TargetSame1_dist0,TargetUp2_dist1,TargetDown1_dist1,
`
	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRankingDist2(fasta)")
	}

	ref = bytes.NewReader(refData)
	query = bytes.NewReader(queryData)
	queryList := new(bytes.Buffer)
	err = List(ref, query, queryList)
	if err != nil {
		t.Error(err)
	}

	ref = bytes.NewReader(refData)
	target = bytes.NewReader(targetData)
	targetList := new(bytes.Buffer)
	err = List(ref, target, targetList)
	if err != nil {
		t.Error(err)
	}

	qtype = "csv"
	ttype = "csv"
	out = new(bytes.Buffer)

	err = TopRanking(queryList, targetList, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRankingDist2(csv)")
	}
}

func TestTopRankingDist3(t *testing.T) {
	refData := []byte(`>ref
ATGATG
`)
	queryData := []byte(
		`>Query1
ATTATT
`)

	targetData := []byte(`>TargetUp1_dist2
ATGATG
>TargetSame1_dist0
ATTATT
>TargetDown1_dist1
ATTACT
>TargetUp2_dist1
ATGATT
>TargetSide1_dist6
CCCCCC
>TargetSide3_dist3
GTGCTT
>TargetSide2_dist2
ATGCTT
`)
	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)
	target := bytes.NewReader(targetData)
	out := new(bytes.Buffer)
	table := false
	qtype := "fasta"
	ttype := "fasta"
	ignoreArray := make([]string, 0)
	TRsizetotal := 0
	TRsizeup := 0
	TRsizedown := 0
	TRsizeside := 0
	TRsizesame := 0
	TRdistall := 0
	TRdistup := 2
	TRdistdown := 1
	TRdistside := 1
	TRthresholdpair := float32(0.1)
	TRthresholdtarget := 10000
	TRnofill := false
	TRdistpush := 2

	err := TopRanking(query, target, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	desiredResult := `query,closestsame,closestup,closestdown,closestside
Query1,TargetSame1_dist0,TargetUp2_dist1;TargetUp1_dist2,TargetDown1_dist1,TargetSide2_dist2;TargetSide3_dist3
`
	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRankingDist3(fasta)")
	}

	ref = bytes.NewReader(refData)
	query = bytes.NewReader(queryData)
	queryList := new(bytes.Buffer)
	err = List(ref, query, queryList)
	if err != nil {
		t.Error(err)
	}

	ref = bytes.NewReader(refData)
	target = bytes.NewReader(targetData)
	targetList := new(bytes.Buffer)
	err = List(ref, target, targetList)
	if err != nil {
		t.Error(err)
	}

	qtype = "csv"
	ttype = "csv"
	out = new(bytes.Buffer)

	err = TopRanking(queryList, targetList, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRankingDist3(csv)")
	}
}

func TestTopRankingTable1(t *testing.T) {
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
CCCCCC
>TargetSide2
ATGCTT
`)
	ref := bytes.NewReader(refData)
	query := bytes.NewReader(queryData)
	target := bytes.NewReader(targetData)
	out := new(bytes.Buffer)
	table := true
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

	err := TopRanking(query, target, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	desiredResult := `query,direction,distance,target
Query1,same,0,TargetSame1
Query1,up,1,TargetUp2
Query1,up,2,TargetUp1
Query1,down,1,TargetDown1
Query1,side,2,TargetSide2
`
	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRankingTable1(fasta)")
	}

	ref = bytes.NewReader(refData)
	query = bytes.NewReader(queryData)
	queryList := new(bytes.Buffer)
	err = List(ref, query, queryList)
	if err != nil {
		t.Error(err)
	}

	ref = bytes.NewReader(refData)
	target = bytes.NewReader(targetData)
	targetList := new(bytes.Buffer)
	err = List(ref, target, targetList)
	if err != nil {
		t.Error(err)
	}

	qtype = "csv"
	ttype = "csv"
	out = new(bytes.Buffer)

	err = TopRanking(queryList, targetList, ref, out, table,
		qtype, ttype, ignoreArray,
		TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
		TRdistall, TRdistup, TRdistdown, TRdistside,
		TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)
	if err != nil {
		t.Error(err)
	}

	if string(out.Bytes()) != desiredResult {
		t.Errorf("problem in TestTopRankingTable1(csv)")
	}
}
