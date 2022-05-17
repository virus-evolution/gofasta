// Package alphabet provides mappings between codons and
// amino acids
package alphabet

import "errors"

// Translate a nucleotide sequence to a protein sequence
// Codons with ambiguous nucleotides are resolved if it can only possibly
// represent one amino acid
func Translate(nuc string) (string, error) {
	if len(nuc)%3 != 0 {
		return "", errors.New("nucleotide string not divisible by 3")
	}
	translation := ""
	aa := ""
	counter := 0
	CD := MakeCodonDict()
	for i := 0; i < len(nuc); i++ {
		aa = aa + string(nuc[i])
		counter++
		if counter == 3 {
			if t, ok := CD[aa]; ok {
				translation = translation + t
			} else {
				translation = translation + "X"
			}
			counter = 0
			aa = ""
		}
	}

	return translation, nil
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

// MakeDegenDict returns a map from codon (string) to degeneracy code ([]byte, 3)
// TODO: needs testing
// func MakeDegenDict() map[string][]byte {
//
// 	degenAA := make(map[string][]byte)
//
// 	degenAA["TTT"] = []byte{0,0,2}
// 	degenAA["TTC"] = []byte{0,0,2}
// 	degenAA["TTA"] = []byte{2,0,2}
// 	degenAA["TTG"] = []byte{2,0,2}
// 	degenAA["TCT"] = []byte{0,0,4}
// 	degenAA["TCC"] = []byte{0,0,4}
// 	degenAA["TCA"] = []byte{0,0,4}
// 	degenAA["TCG"] = []byte{0,0,4}
// 	degenAA["TAT"] = []byte{0,0,2}
// 	degenAA["TAC"] = []byte{0,0,2}
// 	degenAA["TAA"] = []byte{0,2,3}
// 	degenAA["TAG"] = []byte{0,2,3}
// 	degenAA["TGT"] = []byte{0,0,2}
// 	degenAA["TGC"] = []byte{0,0,2}
// 	degenAA["TGA"] = []byte{0,2,3}
// 	degenAA["TGG"] = []byte{0,0,0}
//
// 	degenAA["CTT"] = []byte{0,0,4}
// 	degenAA["CTC"] = []byte{0,0,4}
// 	degenAA["CTA"] = []byte{2,0,4}
// 	degenAA["CTG"] = []byte{2,0,4}
// 	degenAA["CCT"] = []byte{0,0,4}
// 	degenAA["CCC"] = []byte{0,0,4}
// 	degenAA["CCA"] = []byte{0,0,4}
// 	degenAA["CCG"] = []byte{0,0,4}
// 	degenAA["CAT"] = []byte{0,0,2}
// 	degenAA["CAC"] = []byte{0,0,2}
// 	degenAA["CAA"] = []byte{0,0,2}
// 	degenAA["CAG"] = []byte{0,0,2}
// 	degenAA["CGT"] = []byte{0,0,4}
// 	degenAA["CGC"] = []byte{0,0,4}
// 	degenAA["CGA"] = []byte{2,0,4}
// 	degenAA["CGG"] = []byte{2,0,4}
//
// 	degenAA["ATT"] = []byte{0,0,3}
// 	degenAA["ATC"] = []byte{0,0,3}
// 	degenAA["ATA"] = []byte{0,0,3}
// 	degenAA["ATG"] = []byte{0,0,0}
// 	degenAA["ACT"] = []byte{0,0,4}
// 	degenAA["ACC"] = []byte{0,0,4}
// 	degenAA["ACA"] = []byte{0,0,4}
// 	degenAA["ACG"] = []byte{0,0,4}
// 	degenAA["AAT"] = []byte{0,0,2}
// 	degenAA["AAC"] = []byte{0,0,2}
// 	degenAA["AAA"] = []byte{0,0,2}
// 	degenAA["AAG"] = []byte{0,0,2}
// 	degenAA["AGT"] = []byte{0,0,2}
// 	degenAA["AGC"] = []byte{0,0,2}
// 	degenAA["AGA"] = []byte{2,0,2}
// 	degenAA["AGG"] = []byte{2,0,2}
//
// 	degenAA["GTT"] = []byte{0,0,4}
// 	degenAA["GTC"] = []byte{0,0,4}
// 	degenAA["GTA"] = []byte{0,0,4}
// 	degenAA["GTG"] = []byte{0,0,4}
// 	degenAA["GCT"] = []byte{0,0,4}
// 	degenAA["GCC"] = []byte{0,0,4}
// 	degenAA["GCA"] = []byte{0,0,4}
// 	degenAA["GCG"] = []byte{0,0,4}
// 	degenAA["GAT"] = []byte{0,0,2}
// 	degenAA["GAC"] = []byte{0,0,2}
// 	degenAA["GAA"] = []byte{0,0,2}
// 	degenAA["GAG"] = []byte{0,0,2}
// 	degenAA["GGT"] = []byte{0,0,4}
// 	degenAA["GGC"] = []byte{0,0,4}
// 	degenAA["GGA"] = []byte{0,0,4}
// 	degenAA["GGG"] = []byte{0,0,4}
//
// 	return degenAA
// }
