package genbank

import (
	"errors"
	"strconv"
	"strings"
)

var (
	locationErr = errors.New("Error parsing Genbank location")
)

// Genbank location field
type Location struct {
	Representation string
}

// Return the Location's representation as a string
func (l Location) String() string {
	return l.Representation
}

func isNumeric(s string) bool {
	_, err := strconv.ParseFloat(s, 64)
	return err == nil
}

// Convert a Genbank location to the list of each genomic position that it spans
// Positions are 1-based inclusive
func (l Location) GetPositions() ([]int, error) {
	s := l.Representation
	if isNested(s) {
		temp, err := unNestRecur(s)
		if len(temp) != 1 || err != nil {
			return []int{}, locationErr
		}
		return temp[0], nil
	} else if isNumeric(string(s[0])) {
		if isRange(s) {
			pos, err := posFromRange(s)
			if err != nil {
				return []int{}, err
			}
			return pos, nil
		} else {
			return []int{}, locationErr
		}
	} else {
		switch s[0:4] {
		case "join":
			pos, err := posFromJoin(s)
			if err != nil {
				return []int{}, err
			}
			return pos, nil
		case "comp":
			pos, err := posFromComp(s)
			if err != nil {
				return []int{}, err
			}
			return pos, nil
		default:
			return []int{}, locationErr
		}
	}
}

// True/false this location is on the reverse strand
func (l Location) IsReverse() (bool, error) {
	pos, err := l.GetPositions()
	if err != nil {
		return true, err
	}
	if pos[0] > pos[len(pos)-1] {
		return true, nil
	}
	return false, nil
}

// // 5'-most position relative to the forward strand
// func (l Location) FivePrimeMost() (int, error) {
// 	pos, err := l.GetPositions()
// 	if err != nil {
// 		return 0, err
// 	}
// 	min := pos[0]
// 	for _, p := range pos {
// 		if p < min {
// 			min = p
// 		}
// 	}
// 	return min, nil
// }

// // 3'-most position relative to the forward strand
// func (l Location) ThreePrimeMost() (int, error) {
// 	pos, err := l.GetPositions()
// 	if err != nil {
// 		return 0, err
// 	}
// 	max := pos[0]
// 	for _, p := range pos {
// 		if p > max {
// 			max = p
// 		}
// 	}
// 	return max, nil
// }

// Recursively unpacks the component regions of a Genbank Location
func unNestRecur(s string) ([][]int, error) {

	fields := splitOnOuterCommas(s)

	this_result := make([][]int, 0)

	for _, f := range fields {
		// base state - if f is not nested, then we only need to calculate the range
		if !isNested(f) {
			switch f[0:4] {
			case "join":
				pos, _ := posFromJoin(f)
				this_result = append(this_result, pos)
				continue
			case "comp":
				pos, _ := posFromComp(f)
				this_result = append(this_result, pos)
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

		inner_result, err := unNestRecur(inner)
		if err != nil {
			return make([][]int, 0), locationErr
		}
		var pos []int

		switch outer {
		case "join()":
			pos = joinFromRanges(inner_result)
		case "complement()":
			if len(inner_result) != 1 {
				return make([][]int, 0), locationErr
			}
			pos = compFromRange(inner_result[0])
		}

		this_result = append(this_result, pos)
	}

	return this_result, nil
}

// Converts a simple "complement(5prime..3prime)" location feature to a list of the genomic
// positions it represents
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

// Converts a simple "join(5prime..3prime, 5prime..3prime, ...)" location feature to a list of the genomic
// positions it represents
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

func isRange(s string) bool {
	return strings.Contains(s, "..")
}

func posFromRange(s string) ([]int, error) {
	fields := strings.Split(s, "..")
	if len(fields) != 2 {
		return []int{}, locationErr
	}
	start, err := strconv.Atoi(fields[0])
	if err != nil {
		return []int{}, err
	}
	stop, err := strconv.Atoi(fields[1])
	if err != nil {
		return []int{}, err
	}
	pos := make([]int, 0)
	for i := start; i <= stop; i++ {
		pos = append(pos, i)
	}
	return pos, nil
}

// Flattens a slice of slices of integers
func joinFromRanges(as [][]int) []int {
	x := make([]int, 0)
	for i := range as {
		x = append(x, as[i]...)
	}
	return x
}

// Reverses a slice of integers
func compFromRange(as []int) []int {
	for i, j := 0, len(as)-1; i < j; i, j = i+1, j-1 {
		as[i], as[j] = as[j], as[i]
	}
	return as
}

// True/false this Genbank Location contains nested features
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

// Divides a Genbank Location feature into top level fields separated by commas.
// Comma-separated fields in subfeatures remain un-split.
func splitOnOuterCommas(s string) []string {
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
