/*
Package variants implements functionality to annotate mutations relative
to a reference sequence for all records in a multiple sequence alignment
in fasta format.
*/
package variants

import (
	"bufio"
	"bytes"
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/alphabet"
	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/genbank"
	"github.com/virus-evolution/gofasta/pkg/gff"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
	"golang.org/x/exp/constraints"
)

type Region struct {
	Name        string // name of feature, if it has one
	Start       int    // 1-based 5'-most position of region on the forward strand, inclusive
	Stop        int    // 1-based 3'-most position of region on the forward strand, inclusive
	Translation string // amino acid sequence of this region if it is CDS
	Strand      int    // values in the set {-1, +1} only (and "0" for a mixture?!)
	Positions   []int  // all the (1-based, unadjusted) positions in order, on the reverse strand if needs be
}

// A Variant is a struct that contains information about one mutation (nuc, amino acid, indel) between
// reference and query
type Variant struct {
	Queryname      string
	RefAl          string
	QueAl          string
	Position       int    // (1-based) genomic location (for an amino acid change, this is the first position of the codon)
	Residue        int    // (1-based) amino acid location
	Changetype     string // one of {nuc,aa,ins,del}
	Feature        string // this should be, for example, the name of the CDS that the thing is in
	Length         int    // for indels
	SNPs           string // if this is an amino acid change, what are the snps
	Representation string
}

// AnnoStructs is for passing groups of Variant structs around with an index which is used to retain input
// order in the output
type AnnoStructs struct {
	Queryname string
	Vs        []Variant
	Idx       int
}

func Variants(msaIn io.Reader, stdin bool, refID string, annoIn io.Reader, annoSuffix string, out io.Writer, start int, end int, aggregate bool, threshold float64, appendSNP bool, threads int) error {

	var (
		ref fasta.EncodedRecord
		err error
	)

	// Find the reference
	// (Have to move the reader back to the beginning of the alignment, because we are scanning through it twice)
	if refID != "" {
		switch x := msaIn.(type) {
		case *os.File:
			if !stdin {
				ref, err = findReference(msaIn, refID)
				if err != nil {
					return err
				}
				_, err = x.Seek(0, io.SeekStart)
				if err != nil {
					return err
				}
			}
		case *bytes.Reader:
			ref, err = findReference(msaIn, refID)
			if err != nil {
				return err
			}
			_, err = x.Seek(0, io.SeekStart)
			if err != nil {
				return err
			}
		}
	}

	next := fasta.NewEncodedAlignmentSource(msaIn, false, false, false)

	// from stdin, the reference must be the first record
	if stdin && refID != "" {
		ref, err = next()
		if err != nil {
			return err
		}
		if ref.ID != refID {
			return errors.New("--reference is not the first record in --msa")
		}
	}

	var (
		cdsregions []Region
		intregions []int
		refToMSA   []int
	)

	switch annoSuffix {
	case "gb":
		gb, err := genbank.ReadGenBank(annoIn)
		if err != nil {
			return err
		}

		// get the reference from the genbank source if required
		if len(ref.Seq) == 0 {
			EA := encoding.MakeEncodingArray()
			encodedrefseq := make([]byte, len(gb.ORIGIN))
			for i := range gb.ORIGIN {
				encodedrefseq[i] = EA[gb.ORIGIN[i]]
			}
			ref = fasta.EncodedRecord{ID: "annotation_fasta", Seq: encodedrefseq}
			os.Stderr.WriteString("using --annotation fasta as reference\n")
		}

		refLenDegapped := len(ref.Decode().Degap().Seq)

		// get a list of CDS + intergenic regions from the genbank file
		cdsregions, intregions, err = RegionsFromGenbank(gb, refLenDegapped)
		if err != nil {
			return err
		}

		// get the offsets accounting for insertions relative to the reference
		refToMSA, _ = GetMSAOffsets(ref.Seq)

		// check that the reference sequence is in the same coordinates as the annotation
		if len(refToMSA) != len(gb.ORIGIN) {
			return errors.New("the degapped reference sequence (" + ref.ID + ") is not the same length as the genbank annotation")
		}

	case "gff":
		gff, err := gff.ReadGFF(annoIn)
		if err != nil {
			return err
		}

		// get the reference from the gff FASTA if required
		if len(ref.Seq) == 0 {
			switch len(gff.FASTA) {
			case 0:
				return errors.New("couldn't find a reference sequence in the --msa or the gff")
			case 1:
				var encodedrefseq []byte
				for _, v := range gff.FASTA {
					encodedrefseq = make([]byte, len(v.Seq))
					EA := encoding.MakeEncodingArray()
					for i := range v.Seq {
						encodedrefseq[i] = EA[v.Seq[i]]
					}
					ref = fasta.EncodedRecord{ID: "annotation_fasta", Seq: encodedrefseq}
				}
				os.Stderr.WriteString("using --annotation fasta as reference\n")
			default:
				return errors.New("more than one sequence in gff ##FASTA section")
			}

		}

		// get the offsets accounting for insertions relative to the reference
		refToMSA, _ = GetMSAOffsets(ref.Seq)

		refSeqDegapped := ref.Decode().Degap().Seq

		// check that the reference sequence is in the same coordinates as the annotation, if the gff
		// file has a ##sequence-region line
		if len(gff.SequenceRegions) > 1 {
			return errors.New("more than one sequence-region in gff header")
		} else if len(gff.SequenceRegions) == 1 {
			for key := range gff.SequenceRegions {
				region := key
				if len(refSeqDegapped) != gff.SequenceRegions[region].End {
					return errors.New("the degapped reference sequence (" + ref.ID + ") is not the same length as the gff annotation")
				}
			}
		}

		// get a list of CDS + intergenic regions from the gff file
		cdsregions, intregions, err = RegionsFromGFF(gff, refSeqDegapped)
		if err != nil {
			return err
		}

	default:
		return errors.New("couldn't tell if --annotation was a .gb or a .gff file")
	}

	work := func(record fasta.EncodedRecord) (AnnoStructs, error) {
		if len(record.Seq) != len(ref.Seq) {
			return AnnoStructs{}, errors.New("Gapped reference sequence and alignment are not the same width")
		}
		return GetVariantsPair(ref.Seq, record.Seq, ref.ID, record.ID, record.Idx, cdsregions, intregions, refToMSA)
	}

	if aggregate {
		agg, err := NewAggregator(out, start, end, appendSNP, ref.ID)
		if err != nil {
			return err
		}
		err = pipeline.Run(context.Background(), threads, next, work, pipeline.FirstError(agg.Add))
		if err != nil {
			return err
		}
		return agg.Write(threshold)
	}

	lw, err := NewLineWriter(out, start, end, appendSNP, ref.ID)
	if err != nil {
		return err
	}
	return pipeline.Run(context.Background(), threads, next, work, pipeline.FirstError(lw.Write))
}

// findReference gets the reference sequence from the msa if it is in there.
// If it isn't, we will try get it from the annotation (in which case there can
// be no insertions relative to the reference in the msa)
func findReference(msaIn io.Reader, referenceID string) (fasta.EncodedRecord, error) {

	var err error

	coding := encoding.MakeEncodingArray()

	s := bufio.NewScanner(msaIn)
	s.Buffer(make([]byte, 0), 1024*1024)

	first := true

	var id string
	var description string
	var seqBuffer []byte
	var line []byte
	var nuc byte
	var width int

	var refRec fasta.EncodedRecord
	refFound := false

	counter := 0

	for s.Scan() {
		line = s.Bytes()

		if first {

			if line[0] != '>' {
				return fasta.EncodedRecord{}, errors.New("badly formatted fasta file")
			}

			description = string(line[1:])
			id = strings.Fields(description)[0]

			if id == referenceID {
				refFound = true
			}

			first = false

		} else if line[0] == '>' {

			if refFound {
				refRec = fasta.EncodedRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
				return refRec, nil
			}

			if counter == 0 {
				width = len(seqBuffer)
			} else if len(seqBuffer) != width {
				return fasta.EncodedRecord{}, errors.New("different length sequences in input file: is this an alignment?")
			}

			counter++
			description = string(line[1:])
			id = strings.Fields(description)[0]
			seqBuffer = make([]byte, 0)

			if id == referenceID {
				refFound = true
			}

		} else {
			encodedLine := make([]byte, len(line))
			for i := range line {
				nuc = coding[line[i]]
				if nuc == 0 {
					return fasta.EncodedRecord{}, fmt.Errorf("invalid nucleotide in fasta file (%s)", string(line[i]))
				}
				encodedLine[i] = nuc
			}
			seqBuffer = append(seqBuffer, encodedLine...)
		}
	}

	err = s.Err()
	if err != nil {
		return fasta.EncodedRecord{}, err
	}

	if refFound {
		refRec = fasta.EncodedRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
	} else {
		return refRec, errors.New("Couldn't find reference (" + referenceID + ") in msa")
	}

	return refRec, nil
}

func RegionsFromGFF(anno gff.GFF, refSeqDegapped string) ([]Region, []int, error) {

	IDed := make(map[string][]gff.Feature)
	other := make([]gff.Feature, 0)
	for _, f := range anno.Features {
		if !(f.Type == "CDS" || f.Type == "mature_protein_region_of_CDS") {
			continue
		}
		if id := f.Attributes["ID"]; len(id) > 0 {
			IDed[id[0]] = append(IDed[id[0]], f)
		} else {
			other = append(other, f)
		}
	}

	tempcds := make([]Region, 0)
	for _, f := range IDed {
		r, err := CDSRegionfromGFF(f, refSeqDegapped)
		if err != nil {
			return []Region{}, []int{}, err
		}
		tempcds = append(tempcds, r)
	}
	for _, f := range other {
		r, err := CDSRegionfromGFF([]gff.Feature{f}, refSeqDegapped)
		if err != nil {
			return []Region{}, []int{}, err
		}
		tempcds = append(tempcds, r)
	}

	// get a slide of positions that are not coding based on everything above
	inter := codes(tempcds, len(refSeqDegapped))

	// then make the final coding regions based on what has a name
	cds := make([]Region, 0)
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

func CDSRegionfromGFF(fs []gff.Feature, refSeqDegapped string) (Region, error) {
	r := Region{}
	r.Name = gffFeatureName(fs[0])
	// every feature must lie within the reference
	for _, f := range fs {
		if f.Start < 1 || f.End < f.Start || f.End > len(refSeqDegapped) {
			return r, errors.New("Error parsing gff: CDS feature " + strconv.Itoa(f.Start) + ".." + strconv.Itoa(f.End) + " is not within the reference sequence (length " + strconv.Itoa(len(refSeqDegapped)) + ")")
		}
	}
	pos := make([]int, 0)
	switch fs[0].Strand {
	case "+":
		for _, f := range fs {
			if f.Strand != "+" {
				return r, errors.New("Error parsing gff: mixed strands within a single ID")
			}
			for i := f.Start; i <= f.End; i++ {
				pos = append(pos, i)
			}
		}
		r.Strand = 1
		r.Positions = pos
		r.Start = gmin(r.Positions)
		r.Stop = gmax(r.Positions)
		refSeqFeat := ""
		for _, p := range r.Positions {
			refSeqFeat = refSeqFeat + string(refSeqDegapped[p-1])
		}
		t, err := alphabet.Translate(strings.ToUpper(refSeqFeat), true)
		if err != nil {
			return r, err
		}
		r.Translation = t
	case "-":
		for j := len(fs) - 1; j >= 0; j-- {
			f := fs[j]
			if f.Strand != "-" {
				return r, errors.New("Error parsing gff: mixed strands within a single ID")
			}
			for i := f.End; i >= f.Start; i-- {
				pos = append(pos, i)
			}
		}
		r.Strand = -1
		r.Positions = pos
		r.Start = gmin(r.Positions)
		r.Stop = gmax(r.Positions)
		refSeqFeat := ""
		for _, p := range r.Positions {
			refSeqFeat = refSeqFeat + string(refSeqDegapped[p-1])
		}
		t, err := alphabet.Translate(alphabet.Complement(strings.ToUpper(refSeqFeat)), true)
		if err != nil {
			return r, err
		}
		r.Translation = t
	default:
		return r, errors.New("Error parsing gff: protein coding feature needs a strand (+ or -), got \"" + fs[0].Strand + "\"")
	}

	return r, nil
}

// gffFeatureName is the name used to report a CDS: its Name attribute, or "" if it has none
// (a CDS without a Name is not reported). Only the first segment of a multi-segment CDS is consulted.
func gffFeatureName(f gff.Feature) string {
	if v := f.Attributes["Name"]; len(v) > 0 {
		return v[0]
	}
	return ""
}

// Parses a genbank flat format file of genome annotations to extract information about the
// the positions of CDS and intergenic regions, in order to annotate mutations within each
func RegionsFromGenbank(gb genbank.Genbank, refLength int) ([]Region, []int, error) {

	if refLength < 0 {
		return []Region{}, []int{}, errors.New("negative reference length")
	}

	cds := make([]Region, 0)
	for _, f := range gb.FEATURES {
		if f.Feature == "CDS" {
			REGION, err := CDSRegionfromGenbank(f)
			if err != nil {
				return []Region{}, []int{}, err
			}
			// every position must lie within the reference
			if len(REGION.Positions) > 0 && (gmin(REGION.Positions) < 1 || gmax(REGION.Positions) > refLength) {
				return []Region{}, []int{}, errors.New("Genbank CDS feature \"" + REGION.Name + "\" is not within the reference sequence (length " + strconv.Itoa(refLength) + ")")
			}
			cds = append(cds, REGION)
		}
	}

	// Get a slice of the intergenic regions
	inter := codes(cds, refLength)

	return cds, inter, nil
}

func CDSRegionfromGenbank(f genbank.GenbankFeature) (Region, error) {

	if !f.HasAttribute("gene") {
		return Region{}, errors.New("No \"gene\" attibute in Genbank CDS feature")
	}
	if !f.HasAttribute("codon_start") {
		return Region{}, errors.New("No \"codon_start\" attibute in Genbank CDS feature")
	}

	if !f.HasAttribute("translation") {
		return Region{}, errors.New("No \"translation\" attibute in Genbank CDS feature")
	}

	r := Region{
		Name:        f.Info["gene"],
		Translation: f.Info["translation"] + "*",
	}

	var err error

	temp, err := f.Location.GetPositions()
	if err != nil {
		return Region{}, err
	}
	codon_start, err := strconv.Atoi(f.Info["codon_start"])
	if err != nil || codon_start < 1 || codon_start > 3 || codon_start > len(temp) {
		return Region{}, errors.New("Genbank CDS feature has an invalid \"codon_start\" attibute: " + f.Info["codon_start"])
	}
	r.Positions = temp[codon_start-1:]
	if len(r.Positions)%3 != 0 {
		return Region{}, alphabet.ErrorCDSNotModThree
	}
	if len(r.Translation) < len(r.Positions)/3 {
		return Region{}, errors.New("Genbank CDS feature \"" + r.Name + "\": the translation is shorter than the CDS")
	}

	r.Start = gmin(r.Positions)
	r.Stop = gmax(r.Positions)

	reverse, err := f.Location.IsReverse()
	if err != nil {
		return Region{}, err
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
func codes(proteincoding []Region, refLength int) []int {

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

// outsideRange is true if pos is before start or after end. A bound that is not > 0 is not applied, so either
// bound can be used on its own.
func outsideRange(pos, start, end int) bool {
	return (start > 0 && pos < start) || (end > 0 && pos > end)
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

// GetMSAOffsets returns two arrays which contain coordinate shifting information that can be used to
// convert reference to msa coordinates, and vice versa
func GetMSAOffsets(refseq []byte) ([]int, []int) {

	// get the length of the reference sequence without any gaps
	degappedLen := 0
	for _, nuc := range refseq {
		// if there is no alignment gap at this site, ++
		if nuc != 244 {
			degappedLen++
		}
	}

	// We get some offsets:
	// 1) refToMSA = the number of bases to add to convert each reference position to MSA coordinates
	// 2) MSAToRef = the number of bases to subtract to convert each MSA position to reference coordinates
	//			   - note that for MSAToRef positions that are gaps in the ref get a 0 in the MSAToRef slice
	gapsum := 0
	refToMSA := make([]int, degappedLen)
	MSAToRef := make([]int, len(refseq))
	for i, nuc := range refseq {
		if nuc == 244 { // 244 represents a gap ('-')
			gapsum++
			continue
		}
		refToMSA[i-gapsum] = gapsum
		MSAToRef[i] = gapsum
	}

	return refToMSA, MSAToRef
}

func GetVariantsPair(ref, query []byte, refID, queryID string, idx int, cdsregions []Region, intregions []int, offsetRefCoord []int) (AnnoStructs, error) {

	AS := AnnoStructs{}

	indels := getIndelsPair(ref, query)
	nucs := getNucsPair(ref, query, intregions, offsetRefCoord)
	AAs := make([]Variant, 0)
	for _, r := range cdsregions {
		AAs = append(AAs, getAAsPair(ref, query, r, offsetRefCoord)...)
	}

	variants := make([]Variant, 0)
	variants = append(variants, indels...)
	variants = append(variants, nucs...)
	variants = append(variants, AAs...)

	// sort the variants
	sort.SliceStable(variants, func(i, j int) bool {
		return variants[i].Position < variants[j].Position || (variants[i].Position == variants[j].Position && variants[i].Changetype < variants[j].Changetype)
	})

	// there might be dups if there was a snp in the region of a join()
	// (deletions that abut the start of the sequence are not reported, see getIndelsPair)
	finalVariants := make([]Variant, 0)
	for i, v := range variants {
		if i > 0 && v == variants[i-1] {
			continue
		}
		finalVariants = append(finalVariants, v)
	}

	// and we're done
	AS = AnnoStructs{Queryname: queryID, Vs: finalVariants, Idx: idx}

	return AS, nil
}

// FormatVariant returns a string representation of a single mutation, the format
// of which varies given its type (aa/nuc/indel)
func FormatVariant(v Variant, appendSNP bool) (string, error) {
	var s string

	switch v.Changetype {
	case "del":
		s = "del:" + strconv.Itoa(v.Position) + ":" + strconv.Itoa(v.Length)
	case "ins":
		s = "ins:" + strconv.Itoa(v.Position) + ":" + strconv.Itoa(v.Length)
	case "nuc":
		s = "nuc:" + v.RefAl + strconv.Itoa(v.Position) + v.QueAl
	case "aa":
		if appendSNP {
			s = "aa:" + v.Feature + ":" + v.RefAl + strconv.Itoa(v.Residue) + v.QueAl + "(" + v.SNPs + ")"
		} else {
			s = "aa:" + v.Feature + ":" + v.RefAl + strconv.Itoa(v.Residue) + v.QueAl
		}
	default:
		return "", errors.New("couldn't parse variant type")
	}

	return s, nil
}

// LineWriter writes each query's mutations as one line of CSV.
type LineWriter struct {
	w          io.Writer
	start, end int
	appendSNP  bool
	refID      string
}

// NewLineWriter writes the CSV header to w and returns a LineWriter that writes to w. Only mutations
// within start..end (either may be -1 for no limit) are written, and the reference itself (the
// record named refID) is not.
func NewLineWriter(w io.Writer, start, end int, appendSNP bool, refID string) (*LineWriter, error) {
	if _, err := w.Write([]byte("query,mutations\n")); err != nil {
		return nil, err
	}
	return &LineWriter{w: w, start: start, end: end, appendSNP: appendSNP, refID: refID}, nil
}

// Write writes one query's mutations.
func (lw *LineWriter) Write(AS AnnoStructs) error {
	if AS.Queryname == lw.refID {
		return nil
	}
	sa := make([]string, 0, len(AS.Vs))
	for _, v := range AS.Vs {
		if outsideRange(v.Position, lw.start, lw.end) {
			continue
		}
		newVar, err := FormatVariant(v, lw.appendSNP)
		if err != nil {
			return err
		}
		sa = append(sa, newVar)
	}
	_, err := lw.w.Write([]byte(AS.Queryname + "," + strings.Join(sa, "|") + "\n"))
	return err
}

// Aggregator counts how many queries have each mutation, and writes their frequencies.
type Aggregator struct {
	w          io.Writer
	start, end int
	appendSNP  bool
	refID      string
	counts     map[Variant]float64
	n          float64
}

// NewAggregator writes the CSV header to w and returns an Aggregator that will write to w. Only
// mutations within start..end (either may be -1 for no limit) are counted, and the reference itself
// (the record named refID) is not counted.
func NewAggregator(w io.Writer, start, end int, appendSNP bool, refID string) (*Aggregator, error) {
	if _, err := w.Write([]byte("mutation,frequency\n")); err != nil {
		return nil, err
	}
	return &Aggregator{w: w, start: start, end: end, appendSNP: appendSNP, refID: refID, counts: make(map[Variant]float64)}, nil
}

// Add counts one query's mutations.
func (a *Aggregator) Add(AS AnnoStructs) error {
	if AS.Queryname == a.refID {
		return nil
	}
	a.n++
	for _, v := range AS.Vs {
		if outsideRange(v.Position, a.start, a.end) {
			continue
		}
		rep, err := FormatVariant(v, a.appendSNP)
		if err != nil {
			return err
		}
		Vskinny := Variant{RefAl: v.RefAl, QueAl: v.QueAl, Position: v.Position, Residue: v.Residue, Changetype: v.Changetype, Feature: v.Feature, Length: v.Length, Representation: rep}
		a.counts[Vskinny]++
	}
	return nil
}

// Write writes the frequency of each mutation that is at least threshold, ordered by position and then
// by the other fields of the mutation.
func (a *Aggregator) Write(threshold float64) error {
	order := make([]Variant, 0, len(a.counts))
	for k := range a.counts {
		order = append(order, k)
	}

	// The entries came out of a map, so the order must be fully determined by the key: every field
	// that can differ between two entries has to be compared, otherwise ties come out in random order.
	sort.Slice(order, func(i, j int) bool {
		a, b := order[i], order[j]
		switch {
		case a.Position != b.Position:
			return a.Position < b.Position
		case a.Changetype != b.Changetype:
			return a.Changetype < b.Changetype
		case a.QueAl != b.QueAl:
			return a.QueAl < b.QueAl
		case a.Length != b.Length:
			return a.Length < b.Length
		case a.Feature != b.Feature:
			return a.Feature < b.Feature
		case a.Residue != b.Residue:
			return a.Residue < b.Residue
		case a.RefAl != b.RefAl:
			return a.RefAl < b.RefAl
		}
		return a.Representation < b.Representation
	})

	for _, V := range order {
		if a.counts[V]/a.n < threshold {
			continue
		}
		_, err := a.w.Write([]byte(V.Representation + "," + strconv.FormatFloat(a.counts[V]/a.n, 'f', 9, 64) + "\n"))
		if err != nil {
			return err
		}
	}
	return nil
}
