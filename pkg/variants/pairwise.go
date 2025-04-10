package variants

import (
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/alphabet"
	"github.com/virus-evolution/gofasta/pkg/encoding"
)

func getIndelsPair(ref, query []byte, offsetRefCoord []int, offsetMSACoord []int) []Variant {

	var (
		insOpen   bool
		insStart  int
		insLength int
		delOpen   bool
		delStart  int
		delLength int
	)

	variants := make([]Variant, 0)

	for pos := range ref {
		if ref[pos] == 244 { // insertion relative to reference (somewhere in the alignment)
			if query[pos] == 244 { // insertion is not in this seq
				continue
			} else { // insertion is in this seq
				if insOpen { // not the first position of an insertion
					insLength++ // we increment the length counter
				} else { // the first position of an insertion
					insStart = pos // we record the first position of the insertion 0-based in alignment coordinates
					insLength = 1
					insOpen = true
				}
			}
		} else { // not an insertion relative to the reference at this position
			if insOpen { // first base after an insertion, so we need to log the insertion
				variants = append(variants, Variant{Changetype: "ins", Position: (insStart - offsetMSACoord[insStart]), Length: insLength})
				insOpen = false
			}
			if query[pos] == 244 { // deletion in this seq
				if delOpen { // not the first position of a deletion
					delLength++ // we increment the length (there is not a deletion in the reference)
				} else { // the first position of a deletion
					delStart = pos
					delLength = 1
					delOpen = true
				}
			} else { // no deletion in this seq
				if delOpen { // first base after a deletion, so we need to log the deletion
					if delStart-offsetMSACoord[delStart] != 0 {
						variants = append(variants, Variant{Changetype: "del", Position: (delStart - offsetMSACoord[delStart]) + 1, Length: delLength})
					}
					delOpen = false // and reset things
				}
			}
		}
	}

	// don't want deletions at the end of the alignment either
	// if delOpen {
	// 	variants = append(variants, Variant{Changetype: "del", Position: delStart - offsetMSACoord[delStart], Length: delLength})
	// }
	// catch insertions that abut the end of the alignment
	if insOpen {
		variants = append(variants, Variant{Changetype: "ins", Position: (insStart - offsetMSACoord[insStart]) + 1, Length: insLength})
	}

	return variants
}

func getNucsPair(ref, query []byte, pos []int, offsetRefCoord []int, offsetMSACoord []int) []Variant {
	DA := encoding.MakeDecodingArray()
	variants := make([]Variant, 0)
	for _, p := range pos {
		alignPos := (p - 1) + offsetRefCoord[p-1]
		if (ref[alignPos] & query[alignPos]) < 16 { // check for SNPs
			variants = append(variants, Variant{Changetype: "nuc", RefAl: DA[ref[alignPos]], QueAl: DA[query[alignPos]], Position: p})
		}
	}
	return variants
}

// a version of the function that uses the Positions slice in the Region instead of the CodonStarts
func getAAsPair(ref, query []byte, region Region, offsetRefCoord []int, offsetMSACoord []int) []Variant {

	DA := encoding.MakeDecodingArray()
	CD := alphabet.MakeCodonDict()

	variants := make([]Variant, 0)
	codonSNPs := make([]Variant, 0, 3)

	var decodedCodon, refDecodedCodon, aa, refaa string // Add refDecodedCodon
	aaCounter := 0
	codonCounter := 0

	//
	for _, refPos := range region.Positions {
		// here is the actual position in the msa:
		alignmentPos := (refPos - 1) + offsetRefCoord[refPos-1]

		// skip insertions relative to the reference
		// TO DO = if they are in this record log/take them into account in terms of the protein sequence?
		if ref[alignmentPos] == 244 {
			continue
		}
		if (query[alignmentPos] & ref[alignmentPos]) < 16 {
			codonSNPs = append(codonSNPs, Variant{Changetype: "nuc", RefAl: DA[ref[alignmentPos]], QueAl: DA[query[alignmentPos]], Position: refPos})
			// IF ON THE REVERSE STRAND- WHICH NUCLEOTIDE SHOULD BE REPRESENTED IN THE NUC VARIANT?
		}
		decodedCodon = decodedCodon + DA[query[alignmentPos]]
		refDecodedCodon = refDecodedCodon + DA[ref[alignmentPos]] // Build reference codon

		codonCounter++

		if codonCounter == 3 {
			// if reverse strand, take complement (if we're using the slice of nucleotide positions,
			// we're already going in the reverse direction, so don't need the reverse complement)
			if region.Strand == -1 {
				decodedCodon = alphabet.Complement(decodedCodon)
				refDecodedCodon = alphabet.Complement(refDecodedCodon) // Complement reference codon if reverse strand
			}

			if _, ok := CD[decodedCodon]; ok {
				aa = CD[decodedCodon]
			} else {
				aa = "X"
			}

			refaa = string(region.Translation[aaCounter])

			if aa != refaa && aa != "X" {
				temp := []string{}
				for _, v := range codonSNPs {
					temp = append(temp, "nuc:"+v.RefAl+strconv.Itoa(v.Position)+v.QueAl)
				}
				variants = append(variants, Variant{Changetype: "aa", Feature: region.Name, RefAl: refaa, QueAl: aa, Position: refPos - (2 * region.Strand), Residue: aaCounter + 1, SNPs: strings.Join(temp, ";"), RefCodon: refDecodedCodon, QueCodon: decodedCodon}) // Add codons to Variant

			} else {
				for _, v := range codonSNPs {
					variants = append(variants, v)
				}
			}

			codonSNPs = make([]Variant, 0, 3)
			decodedCodon = ""
			refDecodedCodon = "" // Reset reference codon
			codonCounter = 0
			aaCounter++
		}
	}

	return variants
}
