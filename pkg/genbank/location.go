package genbank

import (
	"strconv"
	"strings"
)

type Location struct {
	Representation string
}

func (l Location) String() string {
	return l.Representation
}

func posFromComp(s string) ([]int, error) {
	s = strings.TrimRight(strings.TrimLeft(s, "complement("), ")")
	f := strings.Split(s, "..")
	fivep, err := strconv.Atoi(f[0])
	if err != nil {
		return []int{}, err
	}
	threep, err := strconv.Atoi(f[1])
	if err != nil {
		return []int{}, err
	}
	ps := make([]int, 0)
	for i := threep; i >= fivep; i-- {
		ps = append(ps, i)
	}
	return ps, nil
}

func posFromJoin(s string) ([]int, error) {
	s = strings.TrimRight(strings.TrimLeft(s, "join("), ")")
	j := strings.Split(s, ",")
	ps := make([]int, 0)
	for k := range j {
		temp := make([]int, 0)
		f := strings.Split(j[k], "..")
		fivep, err := strconv.Atoi(f[0])
		if err != nil {
			return []int{}, err
		}
		threep, err := strconv.Atoi(f[1])
		if err != nil {
			return []int{}, err
		}
		for i := fivep; i <= threep; i++ {
			temp = append(temp, i)
		}
		ps = append(ps, temp...)
	}
	return ps, nil
}

func isNested(s string) bool {
	open_count := 0
	for _, r := range s {
		if r == '(' {
			open_count++
		} else if r == ')' {
			open_count--
		}
		if open_count > 1 {
			return true
		}
	}
	return false
}

func SplitOnOuterCommas(s string) []string {
	commas := make([]int, 0)
	open_count := 0
	for i, r := range s {
		if r == '(' {
			open_count++
		} else if r == ')' {
			open_count--
		} else if r == ',' {
			if open_count == 0 {
				commas = append(commas, i)
			}
		}
	}

	starts := []int{0}
	stops := make([]int, 0)
	for _, x := range commas {
		starts = append(starts, x+1)
		stops = append(stops, x)
	}
	stops = append(stops, len(s))

	fields := make([]string, 0)
	for i := 0; i < len(starts); i++ {
		fields = append(fields, s[starts[i]:stops[i]])
	}
	return fields
}

func joinFromRanges(as [][]int) []int {
	x := make([]int, 0)
	for i := range as {
		x = append(x, as[i]...)
	}
	return x
}

func compFromRange(as []int) []int {
	for i, j := 0, len(as)-1; i < j; i, j = i+1, j-1 {
		as[i], as[j] = as[j], as[i]
	}
	return as
}

type positionBuilder struct {
	finished   []int
	inprogress [][]int
}

func (pb *positionBuilder) join() {
	pb.finished = append(pb.finished, joinFromRanges(pb.inprogress)...)
	pb.inprogress = make([][]int, 0)
}

func (pb *positionBuilder) complement() {
	if len(pb.inprogress) != 1 {
		panic("arghh!!")
	}
	pb.finished = append(pb.finished, compFromRange(pb.inprogress[0])...)
	pb.inprogress = make([][]int, 0)
}

func unNestRecur(s string, pb *positionBuilder) {

	fields := SplitOnOuterCommas(s)

	for _, f := range fields {

		// base state - if s is not nested, then we can record the range?
		if !isNested(f) {
			switch f[0:4] {
			case "join":
				pos, _ := posFromJoin(f)
				pb.inprogress = append(pb.inprogress, pos)
				continue
			case "comp":
				pos, _ := posFromComp(f)
				pb.inprogress = append(pb.inprogress, pos)
				continue
			}
		}

		open, closed, open_idx, closed_idx := 0, 0, 0, 0
		for i, r := range f {
			if open == 1 && open_idx == 0 {
				open_idx = i
			}
			if r == '(' {
				open++
			} else if r == ')' {
				closed++
				if open == closed {
					closed_idx = i
				}
			}
		}
		outer := f[0:open_idx] + f[closed_idx:]
		inner := f[open_idx:closed_idx]

		unNestRecur(inner, pb)

		switch outer {
		case "join()":
			pb.join()
		case "complement()":
			pb.complement()
		}
	}
}

func UnNest(l Location, pb *positionBuilder) []int {
	s := l.Representation
	var pos []int
	switch isNested(s) {
	case true:
		unNestRecur(s, pb)
		pos = pb.finished
	case false:
		switch s[0:4] {
		case "join":
			pos, _ = posFromJoin(s)
		case "comp":
			pos, _ = posFromComp(s)
		}
	}
	return pos
}

/*
join(12..78,134..202)
complement(34..126)
complement(join(2691..4571,4918..5163))
join(complement(4918..5163),complement(2691..4571))
join(complement(join(4918..5163,1..2)),complement(2691..4571))
*/
