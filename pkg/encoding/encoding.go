package encoding

// MakeByteDict is a map from the rune/int32 representation of IUPAC nucleotide
// codes to uint8 values according to the (modified) bit-level coding scheme
// for nucleotides developed by Emmanuel Paradis:
// http://ape-package.ird.fr/misc/BitLevelCodingScheme.html
func MakeByteDict() map[rune]uint8 {

	byteMap := make(map[rune]uint8)

	byteMap[65] = 136
	byteMap[71] = 72
	byteMap[67] = 40
	byteMap[84] = 24
	byteMap[82] = 192
	byteMap[77] = 160
	byteMap[87] = 144
	byteMap[83] = 96
	byteMap[75] = 80
	byteMap[89] = 48
	byteMap[86] = 224
	byteMap[72] = 176
	byteMap[68] = 208
	byteMap[66] = 112
	byteMap[78] = 240
	byteMap[45] = 244
	byteMap[63] = 242

	return byteMap
}

// MakeByteDict2 maps from bytes not runes
func MakeByteDict2() map[byte]uint8 {

	byteMap := make(map[byte]uint8)

	byteMap['A'] = 136
	byteMap['G'] = 72
	byteMap['C'] = 40
	byteMap['T'] = 24
	byteMap['R'] = 192
	byteMap['M'] = 160
	byteMap['W'] = 144
	byteMap['S'] = 96
	byteMap['K'] = 80
	byteMap['Y'] = 48
	byteMap['V'] = 224
	byteMap['H'] = 176
	byteMap['D'] = 208
	byteMap['B'] = 112
	byteMap['N'] = 240
	byteMap['-'] = 244
	byteMap['?'] = 242

	return byteMap
}

// MakeEncodingArray returns an array whose indices are the byte representations
// of IUPAC codes and whose contents are Emmanual Paradis encodings
// Lower case nucleotides are mapped to their upper case nucleotides's encoding
func MakeEncodingArray() []byte {
	byteArray := make([]byte, 256)

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

// MakeScoreDict is a map from uint8 values for IUPAC nucleotide codes to an
// integer for how unambiguous they are. The score is caculated as:
// 12 * 1/possible real nucleotides. E.g. an 'A' / 136::uint8 scores 12, but an
// 'N' (A, C, G or T) / 240::uint8 scores 3
func MakeScoreDict() map[uint8]int {

	scoreMap := make(map[uint8]int)

	scoreMap[136] = 12
	scoreMap[72] = 12
	scoreMap[40] = 12
	scoreMap[24] = 12
	scoreMap[192] = 6
	scoreMap[160] = 6
	scoreMap[144] = 6
	scoreMap[96] = 6
	scoreMap[80] = 6
	scoreMap[48] = 6
	scoreMap[224] = 4
	scoreMap[176] = 4
	scoreMap[208] = 4
	scoreMap[112] = 4
	scoreMap[240] = 3
	scoreMap[244] = 3
	scoreMap[242] = 3

	return scoreMap
}

// MakeNucDict maps from uint8 representation of nucleotides to their IUPAC code
func MakeNucDict() map[uint8]string {

	nucMap := make(map[uint8]string)

	nucMap[136] = "A"
	nucMap[72] = "G"
	nucMap[40] = "C"
	nucMap[24] = "T"
	nucMap[192] = "R"
	nucMap[160] = "M"
	nucMap[144] = "W"
	nucMap[96] = "S"
	nucMap[80] = "K"
	nucMap[48] = "Y"
	nucMap[224] = "V"
	nucMap[176] = "H"
	nucMap[208] = "D"
	nucMap[112] = "B"
	nucMap[240] = "N"
	nucMap[244] = "-"
	nucMap[242] = "?"

	return nucMap
}

// MakeDecodingArray returns an array whose indices are Emmanual Paradis encodings
// of IUPAC codes and whose contents are IUPAC codes as strings
func MakeDecodingArray() []string {
	byteArray := make([]string, 256)

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
	byteArray[244] = "-"
	byteArray[242] = "?"

	return byteArray
}
