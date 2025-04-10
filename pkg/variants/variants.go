/*
Package variants implements functionality to annotate mutations relative
to a reference sequence for all records in a multiple sequence alignment
in fasta format.
*/
package variants

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/virus-evolution/gofasta/pkg/alphabet"
	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/genbank"
	"github.com/virus-evolution/gofasta/pkg/gff"
	"golang.org/x/exp/constraints"
)

type Region struct {
	Whichtype   string // only "protein-coding" for now
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
	RefCodon       string // Reference codon for aa change
	QueCodon       string // Query codon for aa change
}

// AnnoStructs is for passing groups of Variant structs around with an index which is used to retain input
// order in the output
type AnnoStructs struct {
	Queryname string
	Vs        []Variant
	Idx       int
}

func Variants(msaIn io.Reader, stdin bool, refID string, annoIn io.Reader, annoSuffix string, out io.Writer, start int, end int, aggregate bool, threshold float64, appendSNP bool, appendCodons bool, threads int) error {

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

	cMSA := make(chan fasta.EncodedRecord, 50+threads)
	cErr := make(chan error)
	cMSADone := make(chan bool)

	go fasta.StreamEncodeAlignment(msaIn, cMSA, cErr, cMSADone, false, false, false)

	firstmissing := false

	if stdin && refID != "" {
		select {
		case ref = <-cMSA:
			if ref.ID != refID {
				return errors.New("--reference is not the first record in --msa")
			}
			firstmissing = true
		case err := <-cErr:
			return err
		case <-cMSADone:
			return errors.New("is the pipe to --msa empty?") // TO DO - does this work/is this necessary?
		}
	}

	var (
		cdsregions         []Region
		intregions         []int
		refToMSA, MSAToRef []int
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
		refToMSA, MSAToRef = GetMSAOffsets(ref.Seq)

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
		refToMSA, MSAToRef = GetMSAOffsets(ref.Seq)

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

	cVariants := make(chan AnnoStructs, 50+threads)
	cVariantsDone := make(chan bool)
	cWriteDone := make(chan bool)

	switch aggregate {
	case true:
		go AggregateWriteVariants(out, start, end, appendSNP, appendCodons, threshold, ref.ID, cVariants, cWriteDone, cErr)
	case false:
		go WriteVariants(out, start, end, firstmissing, appendSNP, appendCodons, ref.ID, cVariants, cWriteDone, cErr)
	}

	var wgVariants sync.WaitGroup
	wgVariants.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			getVariants(ref, cdsregions, intregions, refToMSA, MSAToRef, cMSA, cVariants, cErr)
			wgVariants.Done()
		}()
	}

	go func() {
		wgVariants.Wait()
		cVariantsDone <- true
	}()

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cMSADone:
			close(cMSA)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cVariantsDone:
			close(cVariants)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cWriteDone:
			n--
		}
	}

	return nil
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
	r := Region{
		Whichtype: "protein-coding",
	}
	// TO DO - check that all CDS features in this group have the same "Name"
	// attribute (or none at all). At the moment only the first CDS line's Name
	// is used
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
		t, err := alphabet.Translate(refSeqFeat, true)
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
		t, err := alphabet.Translate(alphabet.Complement(refSeqFeat), true)
		if err != nil {
			return r, err
		}
		r.Translation = t
	case ".":
		return r, errors.New("Error parsing gff: protein coding feature needs a strand")
	}

	return r, nil
}

// Parses a genbank flat format file of genome annotations to extract information about the
// the positions of CDS and intergenic regions, in order to annotate mutations within each
func RegionsFromGenbank(gb genbank.Genbank, refLength int) ([]Region, []int, error) {

	cds := make([]Region, 0)
	for _, f := range gb.FEATURES {
		if f.Feature == "CDS" {
			REGION, err := CDSRegionfromGenbank(f)
			if err != nil {
				return []Region{}, []int{}, err
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

	r := Region{
		Whichtype:   "protein-coding",
		Name:        f.Info["gene"],
		Translation: f.Info["translation"] + "*",
	}

	var err error

	temp, err := f.Location.GetPositions()
	if err != nil {
		return Region{}, err
	}
	codon_start, err := strconv.Atoi(f.Info["codon_start"])
	r.Positions = temp[codon_start-1:]
	if len(r.Positions)%3 != 0 {
		return Region{}, alphabet.ErrorCDSNotModThree
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

// getVariants annotates mutations between query and reference sequences, one
// fasta record at a time. It reads each fasta record from a channel and passes
// all its mutations grouped together in one struct to another channel.
func getVariants(ref fasta.EncodedRecord, cdsregions []Region, intregions []int, offsetRefCoord []int, offsetMSACoord []int, cMSA chan fasta.EncodedRecord, cVariants chan AnnoStructs, cErr chan error) {

	for record := range cMSA {

		if len(record.Seq) != len(offsetMSACoord) {
			cErr <- errors.New("Gapped reference sequence and alignment are not the same width")
			break
		}

		AS, err := GetVariantsPair(ref.Seq, record.Seq, ref.ID, record.ID, record.Idx, cdsregions, intregions, offsetRefCoord, offsetMSACoord)
		if err != nil {
			cErr <- err
			break
		}

		cVariants <- AS
	}
}

func GetVariantsPair(ref, query []byte, refID, queryID string, idx int, cdsregions []Region, intregions []int, offsetRefCoord []int, offsetMSACoord []int) (AnnoStructs, error) {

	AS := AnnoStructs{}

	indels := getIndelsPair(ref, query, offsetRefCoord, offsetMSACoord)
	nucs := getNucsPair(ref, query, intregions, offsetRefCoord, offsetMSACoord)
	AAs := make([]Variant, 0)
	for _, r := range cdsregions {
		AAs = append(AAs, getAAsPair(ref, query, r, offsetRefCoord, offsetMSACoord)...)
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
	finalVariants := make([]Variant, 0)
	previousVariant := Variant{}
	for i, v := range variants {
		if i == 0 {
			// don't want deletions that abut the start of the sequence
			if v.Changetype == "del" && v.Position == 0 {
				continue
			}
			finalVariants = append(finalVariants, v)
			previousVariant = v
			continue
		}
		if v.Changetype == "del" && v.Position == 0 {
			continue
		}
		if v == previousVariant {
			continue
		}
		finalVariants = append(finalVariants, v)
		previousVariant = v
	}

	// and we're done
	AS = AnnoStructs{Queryname: queryID, Vs: finalVariants, Idx: idx}

	return AS, nil
}

// FormatVariant returns a string representation of a single mutation, the format
// of which varies given its type (aa/nuc/indel)
func FormatVariant(v Variant, appendSNP bool, appendCodons bool) (string, error) {
	var s string

	switch v.Changetype {
	case "del":
		s = "del:" + strconv.Itoa(v.Position) + ":" + strconv.Itoa(v.Length)
	case "ins":
		s = "ins:" + strconv.Itoa(v.Position) + ":" + strconv.Itoa(v.Length)
	case "nuc":
		s = "nuc:" + v.RefAl + strconv.Itoa(v.Position) + v.QueAl
	case "aa":
		if appendSNP && appendCodons {
			s = "aa:" + v.Feature + ":" + v.RefAl + strconv.Itoa(v.Residue) + v.QueAl + "(" + v.SNPs + ")" + "(" + v.RefCodon + "->" + v.QueCodon + ")"
		} else if appendSNP {
			s = "aa:" + v.Feature + ":" + v.RefAl + strconv.Itoa(v.Residue) + v.QueAl + "(" + v.SNPs + ")"
		} else if appendCodons {
			s = "aa:" + v.Feature + ":" + v.RefAl + strconv.Itoa(v.Residue) + v.QueAl + "(" + v.RefCodon + "->" + v.QueCodon + ")"
		} else {
			s = "aa:" + v.Feature + ":" + v.RefAl + strconv.Itoa(v.Residue) + v.QueAl
		}
	default:
		return "", errors.New("couldn't parse variant type")
	}

	return s, nil
}

// WriteVariants writes each query's mutations to file or stdout
func WriteVariants(w io.Writer, start, end int, firstmissing bool, appendSNP bool, appendCodons bool, refID string, cVariants chan AnnoStructs, cWriteDone chan bool, cErr chan error) {

	outputMap := make(map[int]AnnoStructs)

	var counter int
	switch firstmissing {
	case true:
		counter = 1
	case false:
		counter = 0
	}

	var err error
	var sa []string

	_, err = w.Write([]byte("query,mutations\n"))
	if err != nil {
		cErr <- err
		return
	}

	for variantLine := range cVariants {
		outputMap[variantLine.Idx] = variantLine

		for {
			if VL, ok := outputMap[counter]; ok {

				if VL.Queryname == refID {
					delete(outputMap, counter)
					counter++
					continue
				}

				_, err = w.Write([]byte(VL.Queryname + ","))
				if err != nil {
					cErr <- err
					return
				}
				sa = make([]string, 0)
				for _, v := range VL.Vs {
					if start > 0 && end > 0 {
						if v.Position < start || v.Position > end {
							continue
						}
					}
					newVar, err := FormatVariant(v, appendSNP, appendCodons)
					if err != nil {
						cErr <- err
						return
					}
					sa = append(sa, newVar)
				}
				_, err = w.Write([]byte(strings.Join(sa, "|") + "\n"))
				if err != nil {
					cErr <- err
					return
				}

				delete(outputMap, counter)
				counter++
			} else {
				break
			}
		}

	}

	cWriteDone <- true
}

// AggregateWriteOutput aggregates the mutations that are present greater than
// or equal to threshold, and writes their frequencies to file or stdout
func AggregateWriteVariants(w io.Writer, start, end int, appendSNP bool, appendCodons bool, threshold float64, refID string, cVariants chan AnnoStructs, cWriteDone chan bool, cErr chan error) {

	propMap := make(map[Variant]float64)

	var err error

	_, err = w.Write([]byte("mutation,frequency\n"))
	if err != nil {
		cErr <- err
		return
	}

	counter := 0.0

	for AS := range cVariants {
		if AS.Queryname == refID {
			continue
		}
		counter++
		for _, v := range AS.Vs {
			if start > 0 && end > 0 {
				if v.Position < start || v.Position > end {
					continue
				}
			}
			rep, err := FormatVariant(v, appendSNP, appendCodons)
			if err != nil {
				cErr <- err
				return
			}
			Vskinny := Variant{RefAl: v.RefAl, QueAl: v.QueAl, Position: v.Position, Residue: v.Residue, Changetype: v.Changetype, Feature: v.Feature, Length: v.Length, Representation: rep}
			propMap[Vskinny]++
		}
	}

	order := make([]Variant, 0)
	for k := range propMap {
		order = append(order, k)
	}

	sort.SliceStable(order, func(i, j int) bool {
		return order[i].Position < order[j].Position || (order[i].Position == order[j].Position && order[i].Changetype < order[j].Changetype) || (order[i].Position == order[j].Position && order[i].Changetype == order[j].Changetype && order[i].QueAl < order[j].QueAl)
	})

	for _, V := range order {
		if propMap[V]/counter < threshold {
			continue
		}
		_, err = w.Write([]byte(V.Representation + "," + strconv.FormatFloat(propMap[V]/counter, 'f', 9, 64) + "\n"))
		if err != nil {
			cErr <- err
			return
		}
	}

	cWriteDone <- true
}
