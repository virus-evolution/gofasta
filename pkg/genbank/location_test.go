package genbank

import (
	"fmt"
	"reflect"
	"testing"
)

func TestGetPositions(t *testing.T) {
	l := Location{Representation: "1..5"}
	p, err := l.GetPositions()
	if err != nil {
		t.Error(err)
	}
	if !reflect.DeepEqual(p, []int{1, 2, 3, 4, 5}) {
		t.Errorf("Problem in TestGetPositions()")
		fmt.Println(p)
	}

	l = Location{Representation: "join(12..14,18..20)"}
	p, err = l.GetPositions()
	if err != nil {
		t.Error(err)
	}
	if !reflect.DeepEqual(p, []int{12, 13, 14, 18, 19, 20}) {
		t.Errorf("Problem in TestGetPositions()")
		fmt.Println(p)
	}

	l = Location{Representation: "complement(1..5)"}
	p, err = l.GetPositions()
	if err != nil {
		t.Error(err)
	}
	if !reflect.DeepEqual(p, []int{5, 4, 3, 2, 1}) {
		t.Errorf("Problem in TestGetPositions()")
		fmt.Println(p)
	}

	l = Location{Representation: "complement(join(1..5,8..10))"}
	p, err = l.GetPositions()
	if err != nil {
		t.Error(err)
	}
	if !reflect.DeepEqual(p, []int{10, 9, 8, 5, 4, 3, 2, 1}) {
		t.Errorf("Problem in TestGetPositions()")
		fmt.Println(p)
	}

	l = Location{Representation: "join(complement(8..10),complement(1..5))"}
	p, err = l.GetPositions()
	if err != nil {
		t.Error(err)
	}
	if !reflect.DeepEqual(p, []int{10, 9, 8, 5, 4, 3, 2, 1}) {
		t.Errorf("Problem in TestGetPositions()")
		fmt.Println(p)
	}

	l = Location{Representation: "join(complement(join(11..12,15..17)),complement(1..3))"}
	p, err = l.GetPositions()
	if err != nil {
		t.Error(err)
	}
	if !reflect.DeepEqual(p, []int{17, 16, 15, 12, 11, 3, 2, 1}) {
		t.Errorf("Problem in TestGetPositions()")
		fmt.Println(p)
	}
}

func TestIsReverse(t *testing.T) {
	l := Location{Representation: "1..5"}
	rev, err := l.IsReverse()
	if err != nil {
		t.Error(err)
	}
	if rev {
		t.Errorf("Problem in TestIsReverse()")
	}

	l = Location{Representation: "complement(join(1..5,8..10))"}
	rev, err = l.IsReverse()
	if err != nil {
		t.Error(err)
	}
	if !rev {
		t.Errorf("Problem in TestIsReverse()")
	}

	l = Location{Representation: "join(complement(join(11..12,15..17)),complement(1..3))"}
	rev, err = l.IsReverse()
	if err != nil {
		t.Error(err)
	}
	if !rev {
		t.Errorf("Problem in TestIsReverse()")
	}
}
