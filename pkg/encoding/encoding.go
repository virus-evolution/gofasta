/*
Package encoding provides mappings between character representations of nucleotides
and bitwise encodings, for fast sequence comparison.

The coding scheme is Emmanual Paradis' design, which is described here:
http://ape-package.ird.fr/misc/BitLevelCodingScheme.html
*/
package encoding

func DecodeToString(bs []byte) string {
	da := MakeDecodingArray()
	s := ""
	for _, n := range bs {
		s = s + da[n]
	}
	return s
}

// MakeEncodingArray returns an array whose indices are the byte representations
// of IUPAC codes and whose contents are Emmanual Paradis encodings.
// Lower case nucleotides are mapped to their upper case nucleotides's encoding
func MakeEncodingArray() [256]byte {
	var byteArray [256]byte

	byteArray['A'] = 136
	byteArray['a'] = 136
	byteArray['G'] = 72
	byteArray['g'] = 72
	byteArray['C'] = 40
	byteArray['c'] = 40
	byteArray['T'] = 24
	byteArray['t'] = 24
	byteArray['R'] = 192
	byteArray['r'] = 192
	byteArray['M'] = 160
	byteArray['m'] = 160
	byteArray['W'] = 144
	byteArray['w'] = 144
	byteArray['S'] = 96
	byteArray['s'] = 96
	byteArray['K'] = 80
	byteArray['k'] = 80
	byteArray['Y'] = 48
	byteArray['y'] = 48
	byteArray['V'] = 224
	byteArray['v'] = 224
	byteArray['H'] = 176
	byteArray['h'] = 176
	byteArray['D'] = 208
	byteArray['d'] = 208
	byteArray['B'] = 112
	byteArray['b'] = 112
	byteArray['N'] = 240
	byteArray['n'] = 240
	byteArray['-'] = 244
	byteArray['?'] = 242

	return byteArray
}

// MakeEncodingArrayHardGaps is as MakeEncodingArray but with '-'
// set to EP's original code (4; character is not completely unknown,
// character represents NONE of ACTG)
func MakeEncodingArrayHardGaps() [256]byte {
	var byteArray [256]byte

	byteArray['A'] = 136
	byteArray['a'] = 136
	byteArray['G'] = 72
	byteArray['g'] = 72
	byteArray['C'] = 40
	byteArray['c'] = 40
	byteArray['T'] = 24
	byteArray['t'] = 24
	byteArray['R'] = 192
	byteArray['r'] = 192
	byteArray['M'] = 160
	byteArray['m'] = 160
	byteArray['W'] = 144
	byteArray['w'] = 144
	byteArray['S'] = 96
	byteArray['s'] = 96
	byteArray['K'] = 80
	byteArray['k'] = 80
	byteArray['Y'] = 48
	byteArray['y'] = 48
	byteArray['V'] = 224
	byteArray['v'] = 224
	byteArray['H'] = 176
	byteArray['h'] = 176
	byteArray['D'] = 208
	byteArray['d'] = 208
	byteArray['B'] = 112
	byteArray['b'] = 112
	byteArray['N'] = 240
	byteArray['n'] = 240
	byteArray['-'] = 4
	byteArray['?'] = 242

	return byteArray
}

// MakeScoreArray makes an array that maps (the byte value of) IUPAC nucleotide
// characters (which are the arrays indices) to a score for how ambiguous that
// nucleotide is. The score is caculated as:
// 12 * 1/possible real nucleotides.
// E.g. an 'A' scores 12, but an 'N' (A, C, G or T) scores 3
func MakeScoreArray() [256]int64 {
	var byteArray [256]int64

	byteArray['A'] = 12
	byteArray['a'] = 12
	byteArray['G'] = 12
	byteArray['g'] = 12
	byteArray['C'] = 12
	byteArray['c'] = 12
	byteArray['T'] = 12
	byteArray['t'] = 12
	byteArray['R'] = 6
	byteArray['r'] = 6
	byteArray['M'] = 6
	byteArray['m'] = 6
	byteArray['W'] = 6
	byteArray['w'] = 6
	byteArray['S'] = 6
	byteArray['s'] = 6
	byteArray['K'] = 6
	byteArray['k'] = 6
	byteArray['Y'] = 6
	byteArray['y'] = 6
	byteArray['V'] = 4
	byteArray['v'] = 4
	byteArray['H'] = 4
	byteArray['h'] = 4
	byteArray['D'] = 4
	byteArray['d'] = 4
	byteArray['B'] = 4
	byteArray['b'] = 4
	byteArray['N'] = 3
	byteArray['n'] = 3
	byteArray['-'] = 3
	byteArray['?'] = 3

	return byteArray
}

// MakeEncodedScoreArray is like MakeScoreArray() but the indices are EP bitwise representations of
// nucleotides
func MakeEncodedScoreArray() [256]int64 {

	var byteArray [256]int64

	byteArray[136] = 12
	byteArray[72] = 12
	byteArray[40] = 12
	byteArray[24] = 12
	byteArray[192] = 6
	byteArray[160] = 6
	byteArray[144] = 6
	byteArray[96] = 6
	byteArray[80] = 6
	byteArray[48] = 6
	byteArray[224] = 4
	byteArray[176] = 4
	byteArray[208] = 4
	byteArray[112] = 4
	byteArray[240] = 3
	byteArray[244] = 3
	byteArray[242] = 3

	return byteArray
}

// MakeDecodingArray returns an array whose indices are Emmanual Paradis encodings
// of IUPAC codes and whose contents are IUPAC codes as strings
func MakeDecodingArray() [256]string {
	var byteArray [256]string

	byteArray[136] = "A"
	byteArray[72] = "G"
	byteArray[40] = "C"
	byteArray[24] = "T"
	byteArray[192] = "R"
	byteArray[160] = "M"
	byteArray[144] = "W"
	byteArray[96] = "S"
	byteArray[80] = "K"
	byteArray[48] = "Y"
	byteArray[224] = "V"
	byteArray[176] = "H"
	byteArray[208] = "D"
	byteArray[112] = "B"
	byteArray[240] = "N"
	byteArray[4] = "-"
	byteArray[244] = "-"
	byteArray[242] = "?"

	return byteArray
}
