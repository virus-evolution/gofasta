package variants

import (
	"errors"
	"sort"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/alphabet"
	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/genbank"
	"github.com/virus-evolution/gofasta/pkg/gff"
	"golang.org/x/exp/constraints"
)

type Region2 struct {
	Whichtype   string // int(ergenic) or protein-coding
	Name        string // name of feature, if it has one
	Start       int    // 1-based 5'-most position of region on the forward strand, inclusive
	Stop        int    // 1-based 3'-most position of region on the forward strand, inclusive
	Translation string // amino acid sequence of this region if it is CDS
	Strand      int    // values in the set {-1, +1} only (and "0" for a mixture?!)
	Positions   []int  // all the (1-based, unadjusted) positions in order, on the reverse strand if needs be
}

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
					insStart = pos // we record the position of the insertion 0-based
					insLength = 1
					insOpen = true
				}
			}
		} else { // not an insertion relative to the reference at this position
			if insOpen { // first base after an insertion, so we need to log the insertion
				variants = append(variants, Variant{Changetype: "ins", Position: insStart - offsetMSACoord[insStart], Length: insLength})
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
						variants = append(variants, Variant{Changetype: "del", Position: delStart - offsetMSACoord[delStart], Length: delLength})
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
		variants = append(variants, Variant{Changetype: "ins", Position: insStart - offsetMSACoord[insStart], Length: insLength})
	}

	return variants
}

func getNucsPair(ref, query []byte, pos []int, refToMSA []int, MSAToRef []int) []Variant {
	DA := encoding.MakeDecodingArray()
	variants := make([]Variant, 0)
	for _, p := range pos {
		alignPos := (p - 1) + refToMSA[p-1]
		if (ref[alignPos] & query[alignPos]) < 16 { // check for SNPs
			variants = append(variants, Variant{Changetype: "nuc", RefAl: DA[ref[alignPos]], QueAl: DA[query[alignPos]], Position: p - 1})
		}
	}
	return variants
}

// a version of the function that uses the Positions slice in the Region instead of the CodonStarts
func getAAsPair(ref, query []byte, region Region2, offsetRefCoord []int, offsetMSACoord []int) []Variant {

	DA := encoding.MakeDecodingArray()
	CD := alphabet.MakeCodonDict()

	variants := make([]Variant, 0)
	codonSNPs := make([]Variant, 0, 3)

	var decodedCodon, aa, refaa string
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
			codonSNPs = append(codonSNPs, Variant{Changetype: "nuc", RefAl: DA[ref[alignmentPos]], QueAl: DA[query[alignmentPos]], Position: refPos - 1})
			// IF ON THE REVERSE STRAND- WHICH NUCLEOTIDE SHOULD BE REPRESENTED IN THE NUC VARIANT?
		}
		decodedCodon = decodedCodon + DA[query[alignmentPos]]

		codonCounter++

		if codonCounter == 3 {
			// if reverse strand, take complement (if we're using the slice of nucleotide positions,
			// we're already going in the reverse direction, so don't need the reverse complement)
			if region.Strand == -1 {
				decodedCodon = alphabet.Complement(decodedCodon)
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
					temp = append(temp, "nuc:"+v.RefAl+strconv.Itoa(v.Position+1)+v.QueAl)
				}
				variants = append(variants, Variant{Changetype: "aa", Feature: region.Name, RefAl: refaa, QueAl: aa, Position: refPos, Residue: aaCounter, SNPs: strings.Join(temp, ";")})

			} else {
				for _, v := range codonSNPs {
					variants = append(variants, v)
				}
			}

			codonSNPs = make([]Variant, 0, 3)
			decodedCodon = ""
			codonCounter = 0
			aaCounter++
		}
	}

	return variants
}

func RegionsFromGFF(anno gff.GFF, refSeqDegapped string) ([]Region2, []int, error) {

	IDed := make(map[string][]gff.Feature)
	other := make([]gff.Feature, 0)
	for _, f := range anno.Features {
		if !(f.Type == "CDS" || f.Type == "mature_protein_region_of_CDS") {
			continue
		}
		if f.HasAttribute("ID") {
			id, ok := f.Attributes["ID"]
			if ok {
				IDed[id[0]] = append(IDed[id[0]], f)
			} else {
				IDed[id[0]] = []gff.Feature{f}
			}
		} else {
			other = append(other, f)
		}
	}

	tempcds := make([]Region2, 0)
	for _, f := range IDed {
		r, err := CDSRegion2fromGFF(f, refSeqDegapped)
		if err != nil {
			return []Region2{}, []int{}, err
		}
		tempcds = append(tempcds, r)
	}
	for _, f := range other {
		r, err := CDSRegion2fromGFF([]gff.Feature{f}, refSeqDegapped)
		if err != nil {
			return []Region2{}, []int{}, err
		}
		tempcds = append(tempcds, r)
	}

	// get a slide of positions that are not coding based on everything above
	inter := codes(tempcds, len(refSeqDegapped))

	// then make the final coding regions based on what has a name
	cds := make([]Region2, 0)
	for _, r := range tempcds {
		if r.Name == "" {
			continue
		}
		cds = append(cds, r)
	}

	// sort by start position
	sort.SliceStable(cds, func(j, k int) bool {
		return cds[j].Start < cds[k].Start
	})

	return cds, inter, nil
}

func CDSRegion2fromGFF(fs []gff.Feature, refSeqDegapped string) (Region2, error) {
	r := Region2{
		Whichtype: "protein-coding",
	}
	if fs[0].HasAttribute("Name") {
		r.Name = fs[0].Attributes["Name"][0]
	} else {
		r.Name = ""
	}
	pos := make([]int, 0)
	switch fs[0].Strand {
	case "+":
		for _, f := range fs {
			if f.Strand != "+" {
				return r, errors.New("Error parsing gff: mixed strands within a single ID")
			}
			for i := f.Start + f.Phase; i <= f.End; i++ {
				pos = append(pos, i)
			}
		}
		r.Strand = 1
	case "-":
		for j := len(fs); j >= 0; j-- {
			f := fs[j]
			if f.Strand != "-" {
				return r, errors.New("Error parsing gff: mixed strands within a single ID")
			}
			for i := f.End - f.Phase; i >= f.Start; i-- {
				pos = append(pos, i)
			}
		}
		r.Strand = -1
	case ".":
		return r, errors.New("Error parsing gff: protein coding feature needs a strand")
	}
	r.Positions = pos
	r.Start = gmin(r.Positions)
	r.Stop = gmax(r.Positions)
	refSeqFeat := ""
	for _, p := range r.Positions {
		refSeqFeat = refSeqFeat + string(refSeqDegapped[p-1])
	}
	t, err := alphabet.Translate(refSeqFeat, true)
	if err != nil {
		return r, err
	}
	r.Translation = t + "*"

	return r, nil
}

// Parses a genbank flat format file of genome annotations to extract information about the
// the positions of CDS and intergenic regions, in order to annotate mutations within each
func RegionsFromGenbank(gb genbank.Genbank, refLength int) ([]Region2, []int, error) {

	cds := make([]Region2, 0)
	for _, f := range gb.FEATURES {
		if f.Feature == "CDS" {
			REGION, err := CDSRegion2fromGenbank(f)
			if err != nil {
				return []Region2{}, []int{}, err
			}
			cds = append(cds, REGION)
		}
	}

	// Get a slice of the intergenic regions
	inter := codes(cds, refLength)

	return cds, inter, nil
}

func CDSRegion2fromGenbank(f genbank.GenbankFeature) (Region2, error) {

	if !f.HasAttribute("gene") {
		return Region2{}, errors.New("No \"gene\" attibute in Genbank CDS feature")
	}
	if !f.HasAttribute("codon_start") {
		return Region2{}, errors.New("No \"codon_start\" attibute in Genbank CDS feature")
	}

	r := Region2{
		Whichtype:   "protein-coding",
		Name:        f.Info["gene"],
		Translation: f.Info["translation"] + "*",
	}

	var err error

	temp, err := f.Location.GetPositions()
	if err != nil {
		return Region2{}, err
	}
	codon_start, err := strconv.Atoi(f.Info["codon_start"])
	r.Positions = temp[codon_start-1:]
	if len(r.Positions)%3 != 0 {
		return Region2{}, alphabet.ErrorCDSNotModThree
	}

	r.Start = gmin(r.Positions)
	r.Stop = gmax(r.Positions)

	reverse, err := f.Location.IsReverse()
	if err != nil {
		return Region2{}, err
	}
	if reverse {
		r.Strand = -1
	} else {
		r.Strand = 1
	}

	return r, nil
}

// get a single slice of intergenic positions after parsing genbank or gff for protein-coding regions
// TO DO - return it in MSA coordinates?
// TO DO - just give it the length of the reference sequence?
func codes(proteincoding []Region2, refLength int) []int {

	// true/false this site codes for a protein:
	codes := make([]bool, refLength) // initialises as falses
	for _, feature := range proteincoding {
		for _, pos := range feature.Positions {
			codes[pos-1] = true
		}
	}
	// then iterate over the genome and retain 1-based positions of sites that don't code for a protein
	intergenicregions := make([]int, 0)
	for i, b := range codes {
		if !b {
			intergenicregions = append(intergenicregions, i+1)
		}
	}

	return intergenicregions
}

// some generic functions because why not
func gmin[T constraints.Ordered](s []T) T {
	if len(s) == 0 {
		panic("Can't obtain the minimum item from a 0-length slice")
	}
	min := s[0]
	for _, x := range s {
		if x < min {
			min = x
		}
	}
	return min
}

func gmax[T constraints.Ordered](s []T) T {
	if len(s) == 0 {
		panic("Can't obtain the maximum item from a 0-length slice")
	}
	max := s[0]
	for _, x := range s {
		if x > max {
			max = x
		}
	}
	return max
}
