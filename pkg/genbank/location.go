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

// maxLocationSpan is the largest number of positions one range may cover. It stops a typo such as
// "1..122222222222222222" from allocating all of memory.
const maxLocationSpan = 100_000_000

// Convert a Genbank location to the list of each genomic position that it spans
// Positions are 1-based inclusive. Supported are ranges ("5..9"), join(...) and complement(...), nested to
// a bounded depth. A complemented location lists its positions in reverse order (3' to 5' on the forward strand).
func (l Location) GetPositions() ([]int, error) {
	pos, err := parseLocation(strings.TrimSpace(l.Representation), 0)
	if err != nil || len(pos) == 0 {
		return []int{}, locationErr
	}
	return pos, nil
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

// maxLocationDepth bounds the nesting of join/complement
const maxLocationDepth = 32

func parseLocation(s string, depth int) ([]int, error) {
	if depth > maxLocationDepth {
		return nil, locationErr
	}
	s = strings.TrimSpace(s)
	switch {
	case strings.HasPrefix(s, "complement(") && strings.HasSuffix(s, ")"):
		pos, err := parseLocation(s[len("complement("):len(s)-1], depth+1)
		if err != nil {
			return nil, err
		}
		for i, j := 0, len(pos)-1; i < j; i, j = i+1, j-1 {
			pos[i], pos[j] = pos[j], pos[i]
		}
		return pos, nil
	case strings.HasPrefix(s, "join(") && strings.HasSuffix(s, ")"):
		fields, err := splitOnOuterCommas(s[len("join(") : len(s)-1])
		if err != nil {
			return nil, err
		}
		pos := make([]int, 0)
		for _, f := range fields {
			p, err := parseLocation(f, depth+1)
			if err != nil {
				return nil, err
			}
			if len(pos)+len(p) > maxLocationSpan {
				return nil, locationErr
			}
			pos = append(pos, p...)
		}
		return pos, nil
	default:
		return posFromRange(s)
	}
}

// posFromRange lists the positions of "start..stop" (1-based, inclusive, start <= stop)
func posFromRange(s string) ([]int, error) {
	fields := strings.Split(s, "..")
	if len(fields) != 2 {
		return nil, locationErr
	}
	start, err := strconv.Atoi(strings.TrimSpace(fields[0]))
	if err != nil {
		return nil, err
	}
	stop, err := strconv.Atoi(strings.TrimSpace(fields[1]))
	if err != nil {
		return nil, err
	}
	if start < 1 || stop < start || stop-start >= maxLocationSpan {
		return nil, locationErr
	}
	pos := make([]int, 0, stop-start+1)
	for i := start; i <= stop; i++ {
		pos = append(pos, i)
	}
	return pos, nil
}

// Divides the inside of a join(...) into top level fields separated by commas.
// Comma-separated fields in subfeatures remain un-split. Unbalanced brackets are an error.
func splitOnOuterCommas(s string) ([]string, error) {
	fields := make([]string, 0)
	open, start := 0, 0
	for i, r := range s {
		switch r {
		case '(':
			open++
		case ')':
			open--
			if open < 0 {
				return nil, locationErr
			}
		case ',':
			if open == 0 {
				fields = append(fields, s[start:i])
				start = i + 1
			}
		}
	}
	if open != 0 {
		return nil, locationErr
	}
	return append(fields, s[start:]), nil
}
