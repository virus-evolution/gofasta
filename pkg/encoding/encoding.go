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
