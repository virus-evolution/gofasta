package encoding

import (
	"strings"
	"testing"
)

func intersectionStringArrays(A []string, B []string) []string {
	intersection := make([]string, 0)
	for i := 0; i < len(A); i++ {
		for j := 0; j < len(B); j++ {
			test := A[i] == B[j]
			if test {
				intersection = append(intersection, A[i])
			}
		}
	}
	return intersection
}

func makeLookupChar() map[byte][]string {
	lookupChar := make(map[byte][]string)

	lookupChar['A'] = []string{"A"}
	lookupChar['a'] = []string{"A"}
	lookupChar['C'] = []string{"C"}
	lookupChar['c'] = []string{"C"}
	lookupChar['G'] = []string{"G"}
	lookupChar['g'] = []string{"G"}
	lookupChar['T'] = []string{"T"}
	lookupChar['t'] = []string{"T"}
	lookupChar['R'] = []string{"A", "G"}
	lookupChar['r'] = []string{"A", "G"}
	lookupChar['Y'] = []string{"C", "T"}
	lookupChar['y'] = []string{"C", "T"}
	lookupChar['S'] = []string{"G", "C"}
	lookupChar['s'] = []string{"G", "C"}
	lookupChar['W'] = []string{"A", "T"}
	lookupChar['w'] = []string{"A", "T"}
	lookupChar['K'] = []string{"G", "T"}
	lookupChar['k'] = []string{"G", "T"}
	lookupChar['M'] = []string{"A", "C"}
	lookupChar['m'] = []string{"A", "C"}
	lookupChar['B'] = []string{"C", "G", "T"}
	lookupChar['b'] = []string{"C", "G", "T"}
	lookupChar['D'] = []string{"A", "G", "T"}
	lookupChar['d'] = []string{"A", "G", "T"}
	lookupChar['H'] = []string{"A", "C", "T"}
	lookupChar['h'] = []string{"A", "C", "T"}
	lookupChar['V'] = []string{"A", "C", "G"}
	lookupChar['v'] = []string{"A", "C", "G"}
	lookupChar['N'] = []string{"A", "C", "G", "T"}
	lookupChar['n'] = []string{"A", "C", "G", "T"}
	lookupChar['?'] = []string{"A", "C", "G", "T"}
	lookupChar['-'] = []string{"A", "C", "G", "T"}

	return lookupChar
}

// tests new encoding using arrays
func TestEncoding(t *testing.T) {

	nucs := []byte{'A', 'G', 'C', 'T', 'R', 'M', 'W', 'S', 'K', 'Y', 'V', 'H', 'D', 'B', 'N', '-', '?',
		'a', 'g', 'c', 't', 'r', 'm', 'w', 's', 'k', 'y', 'v', 'h', 'd', 'b', 'n'}

	lookupChar := makeLookupChar()

	lookupByte := MakeEncodingArray()

	for i := 0; i < len(nucs); i++ {
		for j := 0; j < len(nucs); j++ {
			nuc1 := nucs[i]
			nuc2 := nucs[j]

			nuc1Chars := lookupChar[nuc1]
			nuc2Chars := lookupChar[nuc2]

			byte1 := lookupByte[nuc1]
			byte2 := lookupByte[nuc2]

			byteDifferent := (byte1 & byte2) < 16
			byteSame := (byte1&8 == 8) && byte1 == byte2

			nucDifferent := len(intersectionStringArrays(nuc1Chars, nuc2Chars)) == 0
			nucSame := len(intersectionStringArrays([]string{strings.ToUpper(string(nuc1))}, []string{"A", "C", "G", "T"})) == 1 && strings.ToUpper(string(nuc1)) == strings.ToUpper(string(nuc2))

			test := byteDifferent == nucDifferent && byteSame == nucSame

			if !test {
				t.Errorf("problem in encoding test: %s %s", string(nuc1), string(nuc2))
			}
		}
	}
}

func makeLookupCharHardGaps() map[byte][]string {
	lookupChar := make(map[byte][]string)

	lookupChar['A'] = []string{"A"}
	lookupChar['a'] = []string{"A"}
	lookupChar['C'] = []string{"C"}
	lookupChar['c'] = []string{"C"}
	lookupChar['G'] = []string{"G"}
	lookupChar['g'] = []string{"G"}
	lookupChar['T'] = []string{"T"}
	lookupChar['t'] = []string{"T"}
	lookupChar['R'] = []string{"A", "G"}
	lookupChar['r'] = []string{"A", "G"}
	lookupChar['Y'] = []string{"C", "T"}
	lookupChar['y'] = []string{"C", "T"}
	lookupChar['S'] = []string{"G", "C"}
	lookupChar['s'] = []string{"G", "C"}
	lookupChar['W'] = []string{"A", "T"}
	lookupChar['w'] = []string{"A", "T"}
	lookupChar['K'] = []string{"G", "T"}
	lookupChar['k'] = []string{"G", "T"}
	lookupChar['M'] = []string{"A", "C"}
	lookupChar['m'] = []string{"A", "C"}
	lookupChar['B'] = []string{"C", "G", "T"}
	lookupChar['b'] = []string{"C", "G", "T"}
	lookupChar['D'] = []string{"A", "G", "T"}
	lookupChar['d'] = []string{"A", "G", "T"}
	lookupChar['H'] = []string{"A", "C", "T"}
	lookupChar['h'] = []string{"A", "C", "T"}
	lookupChar['V'] = []string{"A", "C", "G"}
	lookupChar['v'] = []string{"A", "C", "G"}
	lookupChar['N'] = []string{"A", "C", "G", "T"}
	lookupChar['n'] = []string{"A", "C", "G", "T"}
	lookupChar['?'] = []string{"A", "C", "G", "T"}
	lookupChar['-'] = []string{}

	return lookupChar
}

// tests new encoding using arrays
func TestEncodingHardGaps(t *testing.T) {

	nucs := []byte{'A', 'G', 'C', 'T', 'R', 'M', 'W', 'S', 'K', 'Y', 'V', 'H', 'D', 'B', 'N', '-', '?',
		'a', 'g', 'c', 't', 'r', 'm', 'w', 's', 'k', 'y', 'v', 'h', 'd', 'b', 'n'}

	lookupChar := makeLookupCharHardGaps()

	lookupByte := MakeEncodingArrayHardGaps()

	for i := 0; i < len(nucs); i++ {
		for j := 0; j < len(nucs); j++ {
			nuc1 := nucs[i]
			nuc2 := nucs[j]

			nuc1Chars := lookupChar[nuc1]
			nuc2Chars := lookupChar[nuc2]

			byte1 := lookupByte[nuc1]
			byte2 := lookupByte[nuc2]

			byteDifferent := (byte1 & byte2) < 16
			byteSame := (byte1&8 == 8) && byte1 == byte2

			nucDifferent := len(intersectionStringArrays(nuc1Chars, nuc2Chars)) == 0
			nucSame := len(intersectionStringArrays([]string{strings.ToUpper(string(nuc1))}, []string{"A", "C", "G", "T"})) == 1 && strings.ToUpper(string(nuc1)) == strings.ToUpper(string(nuc2))

			test := byteDifferent == nucDifferent && byteSame == nucSame

			if !test {
				t.Errorf("problem in encoding test: %s %s", string(nuc1), string(nuc2))
			}
		}
	}
}

// tests new decoding using maps
func TestDecoding(t *testing.T) {
	nucs := []byte{'A', 'G', 'C', 'T', 'R', 'M', 'W', 'S', 'K', 'Y', 'V', 'H', 'D', 'B', 'N', '-', '?',
		'a', 'g', 'c', 't', 'r', 'm', 'w', 's', 'k', 'y', 'v', 'h', 'd', 'b', 'n'}

	EA := MakeEncodingArray()
	DA := MakeDecodingArray()

	for _, nuc := range nucs {
		a := EA[nuc]
		b := DA[a]
		if strings.ToUpper(string(nuc)) != b {
			t.Errorf("problem in decoding test: %s", string(nuc))
		}
	}

	EA = MakeEncodingArrayHardGaps()
	for _, nuc := range nucs {
		a := EA[nuc]
		b := DA[a]
		if strings.ToUpper(string(nuc)) != b {
			t.Errorf("problem in decoding test: %s", string(nuc))
		}
	}
}

func TestScoring(t *testing.T) {
	SA := MakeScoreArray()
	EA := MakeEncodingArray()
	ESA := MakeEncodedScoreArray()

	lookupChar := makeLookupChar()

	for iupac1, nucs1 := range lookupChar {
		for iupac2, nucs2 := range lookupChar {
			score1 := SA[iupac1]
			score2 := SA[iupac2]
			if float64(score1/score2) != float64(len(nucs2)/len(nucs1)) {
				t.Errorf("problem in scoring test")
			}
		}
	}

	for iupac1, nucs1 := range lookupChar {
		iupac1 = EA[iupac1]
		for iupac2, nucs2 := range lookupChar {
			iupac2 = EA[iupac2]
			score1 := ESA[iupac1]
			score2 := ESA[iupac2]
			if float64(score1/score2) != float64(len(nucs2)/len(nucs1)) {
				t.Errorf("problem in scoring test")
			}
		}
	}
}
