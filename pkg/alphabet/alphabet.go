// Package alphabet provides mappings between codons and
// amino acids
package alphabet

import (
	"errors"

	"github.com/virus-evolution/gofasta/pkg/encoding"
)

var ErrorCDSNotModThree error = errors.New("Translation error: Nucleotide sequence not divisible by 3")

// Translate a nucleotide sequence to a protein sequence
// Codons with ambiguous nucleotides are resolved if it can only possibly
// represent one amino acid
func Translate(nuc string, strict bool) (string, error) {
	if len(nuc)%3 != 0 {
		return "", ErrorCDSNotModThree
	}
	translation := ""
	codon := ""
	counter := 0
	CD := MakeCodonDict()
	for i := 0; i < len(nuc); i++ {
		codon = codon + string(nuc[i])
		counter++
		if counter == 3 {
			if t, ok := CD[codon]; ok {
				translation = translation + t
			} else {
				if strict {
					return "", errors.New("Translation error: Untranslatable codon: " + codon)
				} else {
					translation = translation + "X"
				}
			}
			counter = 0
			codon = ""
		}
	}

	return translation, nil
}

func Complement(nuc string) string {
	CA := MakeCompArray()
	ba := make([]byte, len(nuc))
	for i := 0; i < len(nuc); i++ {
		ba[i] = CA[nuc[i]]
	}
	return string(ba)
}

func ReverseComplement(nuc string) string {
	temp := []byte(Complement(nuc))
	for i, j := 0, len(temp)-1; i < j; i, j = i+1, j-1 {
		temp[i], temp[j] = temp[j], temp[i]
	}
	return string(temp)
}

func MakeCompArray() [256]byte {
	var compArray [256]byte

	compArray['A'] = 'T'
	compArray['a'] = 't'
	compArray['G'] = 'C'
	compArray['g'] = 'c'
	compArray['C'] = 'G'
	compArray['c'] = 'g'
	compArray['T'] = 'A'
	compArray['t'] = 'a'
	compArray['R'] = 'Y'
	compArray['r'] = 'y'
	compArray['Y'] = 'R'
	compArray['y'] = 'r'
	compArray['S'] = 'S'
	compArray['s'] = 's'
	compArray['W'] = 'W'
	compArray['w'] = 'w'
	compArray['K'] = 'M'
	compArray['k'] = 'm'
	compArray['M'] = 'K'
	compArray['m'] = 'k'
	compArray['B'] = 'V'
	compArray['b'] = 'v'
	compArray['V'] = 'B'
	compArray['v'] = 'b'
	compArray['D'] = 'H'
	compArray['d'] = 'h'
	compArray['H'] = 'D'
	compArray['h'] = 'd'
	compArray['N'] = 'N'
	compArray['n'] = 'n'
	compArray['?'] = '?'
	compArray['-'] = '-'

	return compArray
}

func MakeEncodedCompArray() [256]byte {
	var compArray [256]byte

	EA := encoding.MakeEncodingArray()

	compArray[EA['A']] = EA['T']
	compArray[EA['a']] = EA['t']
	compArray[EA['G']] = EA['C']
	compArray[EA['g']] = EA['c']
	compArray[EA['C']] = EA['G']
	compArray[EA['c']] = EA['g']
	compArray[EA['T']] = EA['A']
	compArray[EA['t']] = EA['a']
	compArray[EA['R']] = EA['Y']
	compArray[EA['r']] = EA['y']
	compArray[EA['Y']] = EA['R']
	compArray[EA['y']] = EA['r']
	compArray[EA['S']] = EA['S']
	compArray[EA['s']] = EA['s']
	compArray[EA['W']] = EA['W']
	compArray[EA['w']] = EA['w']
	compArray[EA['K']] = EA['M']
	compArray[EA['k']] = EA['m']
	compArray[EA['M']] = EA['K']
	compArray[EA['m']] = EA['k']
	compArray[EA['B']] = EA['V']
	compArray[EA['b']] = EA['v']
	compArray[EA['V']] = EA['B']
	compArray[EA['v']] = EA['b']
	compArray[EA['D']] = EA['H']
	compArray[EA['d']] = EA['h']
	compArray[EA['H']] = EA['D']
	compArray[EA['h']] = EA['d']
	compArray[EA['N']] = EA['N']
	compArray[EA['n']] = EA['n']
	compArray[EA['?']] = EA['?']
	compArray[EA['-']] = EA['-']

	return compArray
}

// MakeCodonDict returns a map from codon (string) to amino acid code (string)
func MakeCodonDict() map[string]string {

	codonAA := make(map[string]string)

	codonAA["TTT"] = "F"
	codonAA["TTC"] = "F"
	codonAA["TTA"] = "L"
	codonAA["TTG"] = "L"
	codonAA["TCT"] = "S"
	codonAA["TCC"] = "S"
	codonAA["TCA"] = "S"
	codonAA["TCG"] = "S"
	codonAA["TAT"] = "Y"
	codonAA["TAC"] = "Y"
	codonAA["TAA"] = "*"
	codonAA["TAG"] = "*"
	codonAA["TGT"] = "C"
	codonAA["TGC"] = "C"
	codonAA["TGA"] = "*"
	codonAA["TGG"] = "W"

	codonAA["CTT"] = "L"
	codonAA["CTC"] = "L"
	codonAA["CTA"] = "L"
	codonAA["CTG"] = "L"
	codonAA["CCT"] = "P"
	codonAA["CCC"] = "P"
	codonAA["CCA"] = "P"
	codonAA["CCG"] = "P"
	codonAA["CAT"] = "H"
	codonAA["CAC"] = "H"
	codonAA["CAA"] = "Q"
	codonAA["CAG"] = "Q"
	codonAA["CGT"] = "R"
	codonAA["CGC"] = "R"
	codonAA["CGA"] = "R"
	codonAA["CGG"] = "R"

	codonAA["ATT"] = "I"
	codonAA["ATC"] = "I"
	codonAA["ATA"] = "I"
	codonAA["ATG"] = "M"
	codonAA["ACT"] = "T"
	codonAA["ACC"] = "T"
	codonAA["ACA"] = "T"
	codonAA["ACG"] = "T"
	codonAA["AAT"] = "N"
	codonAA["AAC"] = "N"
	codonAA["AAA"] = "K"
	codonAA["AAG"] = "K"
	codonAA["AGT"] = "S"
	codonAA["AGC"] = "S"
	codonAA["AGA"] = "R"
	codonAA["AGG"] = "R"

	codonAA["GTT"] = "V"
	codonAA["GTC"] = "V"
	codonAA["GTA"] = "V"
	codonAA["GTG"] = "V"
	codonAA["GCT"] = "A"
	codonAA["GCC"] = "A"
	codonAA["GCA"] = "A"
	codonAA["GCG"] = "A"
	codonAA["GAT"] = "D"
	codonAA["GAC"] = "D"
	codonAA["GAA"] = "E"
	codonAA["GAG"] = "E"
	codonAA["GGT"] = "G"
	codonAA["GGC"] = "G"
	codonAA["GGA"] = "G"
	codonAA["GGG"] = "G"

	// ambiguities:

	codonAA["AAR"] = "K"
	codonAA["AAY"] = "N"
	codonAA["ACR"] = "T"
	codonAA["ACY"] = "T"
	codonAA["ACS"] = "T"
	codonAA["ACW"] = "T"
	codonAA["ACK"] = "T"
	codonAA["ACM"] = "T"
	codonAA["ACB"] = "T"
	codonAA["ACD"] = "T"
	codonAA["ACH"] = "T"
	codonAA["ACV"] = "T"
	codonAA["ACN"] = "T"
	codonAA["AGR"] = "R"
	codonAA["AGY"] = "S"
	codonAA["ATY"] = "I"
	codonAA["ATW"] = "I"
	codonAA["ATM"] = "I"
	codonAA["ATH"] = "I"
	codonAA["CAR"] = "Q"
	codonAA["CAY"] = "H"
	codonAA["CCR"] = "P"
	codonAA["CCY"] = "P"
	codonAA["CCS"] = "P"
	codonAA["CCW"] = "P"
	codonAA["CCK"] = "P"
	codonAA["CCM"] = "P"
	codonAA["CCB"] = "P"
	codonAA["CCD"] = "P"
	codonAA["CCH"] = "P"
	codonAA["CCV"] = "P"
	codonAA["CCN"] = "P"
	codonAA["CGR"] = "R"
	codonAA["CGY"] = "R"
	codonAA["CGS"] = "R"
	codonAA["CGW"] = "R"
	codonAA["CGK"] = "R"
	codonAA["CGM"] = "R"
	codonAA["CGB"] = "R"
	codonAA["CGD"] = "R"
	codonAA["CGH"] = "R"
	codonAA["CGV"] = "R"
	codonAA["CGN"] = "R"
	codonAA["CTR"] = "L"
	codonAA["CTY"] = "L"
	codonAA["CTS"] = "L"
	codonAA["CTW"] = "L"
	codonAA["CTK"] = "L"
	codonAA["CTM"] = "L"
	codonAA["CTB"] = "L"
	codonAA["CTD"] = "L"
	codonAA["CTH"] = "L"
	codonAA["CTV"] = "L"
	codonAA["CTN"] = "L"
	codonAA["GAR"] = "E"
	codonAA["GAY"] = "D"
	codonAA["GCR"] = "A"
	codonAA["GCY"] = "A"
	codonAA["GCS"] = "A"
	codonAA["GCW"] = "A"
	codonAA["GCK"] = "A"
	codonAA["GCM"] = "A"
	codonAA["GCB"] = "A"
	codonAA["GCD"] = "A"
	codonAA["GCH"] = "A"
	codonAA["GCV"] = "A"
	codonAA["GCN"] = "A"
	codonAA["GGR"] = "G"
	codonAA["GGY"] = "G"
	codonAA["GGS"] = "G"
	codonAA["GGW"] = "G"
	codonAA["GGK"] = "G"
	codonAA["GGM"] = "G"
	codonAA["GGB"] = "G"
	codonAA["GGD"] = "G"
	codonAA["GGH"] = "G"
	codonAA["GGV"] = "G"
	codonAA["GGN"] = "G"
	codonAA["GTR"] = "V"
	codonAA["GTY"] = "V"
	codonAA["GTS"] = "V"
	codonAA["GTW"] = "V"
	codonAA["GTK"] = "V"
	codonAA["GTM"] = "V"
	codonAA["GTB"] = "V"
	codonAA["GTD"] = "V"
	codonAA["GTH"] = "V"
	codonAA["GTV"] = "V"
	codonAA["GTN"] = "V"
	codonAA["TAR"] = "*"
	codonAA["TAY"] = "Y"
	codonAA["TCR"] = "S"
	codonAA["TCY"] = "S"
	codonAA["TCS"] = "S"
	codonAA["TCW"] = "S"
	codonAA["TCK"] = "S"
	codonAA["TCM"] = "S"
	codonAA["TCB"] = "S"
	codonAA["TCD"] = "S"
	codonAA["TCH"] = "S"
	codonAA["TCV"] = "S"
	codonAA["TCN"] = "S"
	codonAA["TGY"] = "C"
	codonAA["TTR"] = "L"
	codonAA["TTY"] = "F"
	codonAA["TRA"] = "*"
	codonAA["YTA"] = "L"
	codonAA["YTG"] = "L"
	codonAA["YTR"] = "L"
	codonAA["MGA"] = "R"
	codonAA["MGG"] = "R"
	codonAA["MGR"] = "R"

	return codonAA
}
