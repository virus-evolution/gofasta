/*
TO DO
*/
package genbank

import (
	"bufio"
	"io"
	"strconv"
	"strings"
	"unicode"
	"unicode/utf8"
)

// Genbank is a master struct containing information from a single genbank record
type Genbank struct {
	LOCUS struct {
		Name     string
		Length   int
		Type     string
		Division string
		Date     string
	} // NOT implemented
	DEFINITION string // NOT implemented
	ACCESSION  string // NOT implemented
	VERSION    string // NOT implemented
	KEYWORDS   string // NOT implemented
	SOURCE     struct {
		Source   string
		Organism string
	} // NOT implemented
	REFERENCE struct {
		Authors string
		Title   string
		Journal string
		Pubmed  string
		Remark  string
	} // NOT implemented
	COMMENT  string           // NOT implemented
	FEATURES []GenbankFeature // implemented
	ORIGIN   []byte           // implemented
}

// genbankField is a utility struct for moving main toplevel genbank FIELDS +
// their associated lines around through channels, etc.
type genbankField struct {
	header string
	lines  []string
}

// GenbankFeature is contains information about one feature
// from a genbank record's FEATURES section
type GenbankFeature struct {
	Feature  string
	Location Location
	Info     map[string]string
}

func (F *GenbankFeature) HasAttribute(tag string) bool {
	if _, ok := F.Info[tag]; ok {
		return true
	}
	return false
}

// isFeatureLine returns true/false does this line of the file code a new FEATURE (CDS, gene, 5'UTR etc)
func isFeatureLine(line string, quoteClosed bool) bool {

	lineFields := strings.Fields(line)

	if quoteClosed {
		if len(lineFields) == 2 {
			if lineFields[0][0] != '/' {
				return true
			}
		}
	}

	return false
}

// parseGenbankFEATURES gets the FEATURES info from a genbank file
func parseGenbankFEATURES(field genbankField) []GenbankFeature {

	rawLines := field.lines

	features := make([]GenbankFeature, 0)

	quoteClosed := true
	var feature string
	var pos string
	var gb GenbankFeature
	var keyBuffer []rune
	var valueBuffer []rune
	var isKey bool

	for linecounter, line := range rawLines {

		newFeature := isFeatureLine(line, quoteClosed)

		if newFeature && linecounter == 0 {

			lineFields := strings.Fields(line)

			feature = lineFields[0]
			pos = lineFields[1]

			gb = GenbankFeature{}
			gb.Feature = feature
			gb.Location = Location{Representation: pos}
			gb.Info = make(map[string]string)

			keyBuffer = make([]rune, 0)
			valueBuffer = make([]rune, 0)

		} else if strings.TrimSpace(line)[0] == '/' && len(keyBuffer) == 0 {

			keyBuffer = make([]rune, 0)
			valueBuffer = make([]rune, 0)

			isKey = true

			quoteClosed = true

			for _, character := range strings.TrimSpace(line)[1:] {

				if character == '=' {
					isKey = false
					continue
				}

				if isKey == true {
					keyBuffer = append(keyBuffer, character)
				} else {
					if character == '"' {
						quoteClosed = !quoteClosed
						continue
					}
					valueBuffer = append(valueBuffer, character)
				}
			}

		} else if !quoteClosed {

			for _, character := range strings.TrimSpace(line) {
				if character == '"' {
					quoteClosed = !quoteClosed
					continue
				}

				valueBuffer = append(valueBuffer, character)
			}

		} else if strings.TrimSpace(line)[0] == '/' && len(keyBuffer) != 0 {

			quoteClosed = true

			gb.Info[string(keyBuffer)] = string(valueBuffer)

			keyBuffer = make([]rune, 0)
			valueBuffer = make([]rune, 0)

			isKey = true

			for _, character := range strings.TrimSpace(line)[1:] {

				if character == '=' {
					isKey = false
					continue
				}

				if isKey {
					keyBuffer = append(keyBuffer, character)
				} else {
					if character == '"' {
						quoteClosed = !quoteClosed
						continue
					}
					valueBuffer = append(valueBuffer, character)
				}
			}

		} else if newFeature && linecounter != 0 {

			quoteClosed = true

			gb.Info[string(keyBuffer)] = string(valueBuffer)
			features = append(features, gb)

			lineFields := strings.Fields(line)
			feature = lineFields[0]
			pos = lineFields[1]

			gb = GenbankFeature{}
			gb.Feature = feature
			gb.Location = Location{Representation: pos}
			gb.Info = make(map[string]string)

			keyBuffer = make([]rune, 0)
			valueBuffer = make([]rune, 0)
		}
	}

	if len(keyBuffer) > 0 && len(valueBuffer) > 0 {
		gb.Info[string(keyBuffer)] = string(valueBuffer)
	}

	features = append(features, gb)

	return features
}

// TO DO - delete this
// ParsePositions returns the genomic positions of a feature, handling regions that are join()ed together
func ParsePositions(position string) ([]int, error) {
	var A []int
	if position[0:4] == "join" {
		A = make([]int, 0)
		position = strings.TrimLeft(position, "join(")
		position = strings.TrimRight(position, ")")
		ranges := strings.Split(position, ",")
		for _, x := range ranges {
			y := strings.Split(x, "..")
			for _, z := range y {
				temp, err := strconv.Atoi(z)
				if err != nil {
					return []int{}, err
				}
				A = append(A, temp)
			}
		}
	} else {
		A = make([]int, 0)
		y := strings.Split(position, "..")
		for _, z := range y {
			temp, err := strconv.Atoi(z)
			if err != nil {
				return []int{}, err
			}
			A = append(A, temp)
		}
	}

	return A, nil
}

// parseGenbankORIGIN gets the ORIGIN info from a genbank file
func parseGenbankORIGIN(field genbankField) []byte {

	rawLines := field.lines

	seq := make([]byte, 0)

	for _, line := range rawLines {
		for _, character := range line {
			if unicode.IsLetter(character) {
				seq = append(seq, []byte(string(character))...)
			}
		}
	}

	return seq
}

// ReadGenBank reads a genbank annotation file and returns a struct that contains
// parsed versions of the fields it contains. Not all fields are currently parsed.
func ReadGenBank(r io.Reader) (Genbank, error) {

	gb := Genbank{}

	s := bufio.NewScanner(r)

	first := true
	var header string
	var lines []string
	var field genbankField

	for s.Scan() {
		line := s.Text()

		if len(line) == 0 {
			continue
		}

		r, _ := utf8.DecodeRune([]byte{line[0]})

		if unicode.IsUpper(r) {
			if first {
				header = strings.Fields(line)[0]
				first = false
				continue
			}

			switch {
			case header == "FEATURES":
				field = genbankField{header: header, lines: lines}
				gb.FEATURES = parseGenbankFEATURES(field)
			case header == "ORIGIN":
				field = genbankField{header: header, lines: lines}
				gb.ORIGIN = parseGenbankORIGIN(field)
			}

			header = strings.Fields(line)[0]
			lines = make([]string, 0)

			continue
		}

		lines = append(lines, line)
	}

	switch {
	case header == "FEATURES":
		field = genbankField{header: header, lines: lines}
		gb.FEATURES = parseGenbankFEATURES(field)
	case header == "ORIGIN":
		field = genbankField{header: header, lines: lines}
		gb.ORIGIN = parseGenbankORIGIN(field)
	}

	return gb, nil
}
