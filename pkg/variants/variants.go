/*
Package variants implements functionality to annotate mutations relative
to a reference sequence for all records in a multiple sequence alignment
(in fasta format).
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

	"github.com/cov-ert/gofasta/pkg/alphabet"
	"github.com/cov-ert/gofasta/pkg/encoding"
	"github.com/cov-ert/gofasta/pkg/fastaio"
	"github.com/cov-ert/gofasta/pkg/genbank"
)

// Region is a struct containing a part of the genome which might
// be a CDS or intergenic, for example
type Region struct {
	Whichtype   string // int(ergenic) or CDS
	Name        string // name of CDS, if it is one
	Start       int    // 1-based first position of region, inclusive
	Stop        int    // 1-based last position of region, inclusive
	Codonstarts []int  // a slice of the 1-based start positions of all its codons, if this region is a CDS
	Translation string // amino acid sequence of this region if it is a CDS
}

// A Variant is a struct that contains information about one mutation (nuc, amino acid, indel) between
// reference and query
type Variant struct {
	Queryname  string
	RefAl      string
	QueAl      string
	Position   int    // genomic location
	Residue    int    // feature location?
	Changetype string // one of {nuc,aa,ins,del}
	Feature    string // this should be, for example, the name of the CDS that the thing is in
	Length     int    // for indels
	SNPs       string // if this is an amino acid change, what are the snps
}

// AnnoStructs is for passing groups of Variant structs around with an index which is used to retain input
// order in the output
type AnnoStructs struct {
	Queryname string
	Vs        []Variant
	Idx       int
}

// Variants annotates amino acid, insertion, deletion, and nucleotide (anything outside of codons represented by an amino acid change)
// mutations relative to a reference sequence from a multiple sequence alignment in fasta format. Genome annotations are derived from a genbank flat file
func Variants(msaIn io.Reader, stdin bool, refID string, gbIn io.Reader, out io.Writer, aggregate bool, threshold float64, appendSNP bool, threads int) error {

	gb, err := genbank.ReadGenBank(gbIn)
	if err != nil {
		return err
	}

	var ref fastaio.EncodedFastaRecord

	// move the reader back to the beginning of the alignment, because we are scanning through it twice
	// NB - is there a more elegant way to do this?
	switch x := msaIn.(type) {
	case *os.File:
		if !stdin && refID != "" {
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

	if len(ref.Seq) == 0 {
		EA := encoding.MakeEncodingArray()
		encodedrefseq := make([]byte, len(gb.ORIGIN))
		for i := range gb.ORIGIN {
			encodedrefseq[i] = EA[gb.ORIGIN[i]]
		}
		ref = fastaio.EncodedFastaRecord{ID: "genbank_source", Seq: encodedrefseq}
	}

	// get a list of CDS + intergenic regions from the genbank file
	regions, err := GetRegions(gb)
	if err != nil {
		return err
	}

	cMSA := make(chan fastaio.EncodedFastaRecord)
	cErr := make(chan error)
	cMSADone := make(chan bool)

	cVariants := make(chan AnnoStructs)
	cVariantsDone := make(chan bool)
	cWriteDone := make(chan bool)

	go fastaio.ReadEncodeAlignment(msaIn, false, cMSA, cErr, cMSADone)

	if stdin && refID != "" {
		select {
		case ref = <-cMSA:
			if ref.ID != refID {
				return errors.New("--reference is not the first record in --msa (if --msa is reference-length you don't need to provide --reference)")
			}
		case err := <-cErr:
			return err
		case <-cMSADone:
			return errors.New("is the pipe to --msa empty?")
		}
	}

	switch aggregate {
	case true:
		go AggregateWriteVariants(out, appendSNP, threshold, refID, cVariants, cWriteDone, cErr)
	case false:
		go WriteVariants(out, appendSNP, refID, cVariants, cWriteDone, cErr)
	}

	// get the offsets accounting for insertions relative to the reference
	refToMSA, MSAToRef := GetMSAOffsets(ref.Seq)

	var wgVariants sync.WaitGroup
	wgVariants.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			getVariants(ref, regions, refToMSA, MSAToRef, cMSA, cVariants, cErr)
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

	// fmt.Println(ref.ID)

	return nil
}

// findReference gets the reference sequence from the msa if it is in there. If it isn't, we will try get it
// from the genbank record (in which case can be no insertions relative to the ref in the msa)
func findReference(msaIn io.Reader, referenceID string) (fastaio.EncodedFastaRecord, error) {

	var err error

	coding := encoding.MakeEncodingArray()

	s := bufio.NewScanner(msaIn)

	first := true

	var id string
	var description string
	var seqBuffer []byte
	var line []byte
	var nuc byte
	var width int

	var refRec fastaio.EncodedFastaRecord
	refFound := false

	counter := 0

	for s.Scan() {
		line = s.Bytes()

		if first {

			if line[0] != '>' {
				return fastaio.EncodedFastaRecord{}, errors.New("badly formatted fasta file")
			}

			description = string(line[1:])
			id = strings.Fields(description)[0]

			if id == referenceID {
				refFound = true
			}

			first = false

		} else if line[0] == '>' {

			if refFound {
				refRec = fastaio.EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
				return refRec, nil
			}

			if counter == 0 {
				width = len(seqBuffer)
			} else if len(seqBuffer) != width {
				return fastaio.EncodedFastaRecord{}, errors.New("different length sequences in input file: is this an alignment?")
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
					return fastaio.EncodedFastaRecord{}, fmt.Errorf("invalid nucleotide in fasta file (%s)", string(line[i]))
				}
				encodedLine[i] = nuc
			}
			seqBuffer = append(seqBuffer, encodedLine...)
		}
	}

	if counter == 0 {
		return fastaio.EncodedFastaRecord{}, errors.New("empty fasta file")
	}

	err = s.Err()
	if err != nil {
		return fastaio.EncodedFastaRecord{}, err
	}

	if refFound {
		refRec = fastaio.EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
	} else {
		os.Stderr.WriteString("no reference sequence found in --msa, attempting to use --genbank source as reference\n")
	}

	return refRec, nil
}

// GetRegions parses a genbank flat format file of genome annotations to extract information about the
// the positions of CDS and intergenic regions, in order to annotate mutations within each
func GetRegions(gb genbank.Genbank) ([]Region, error) {
	CDSFEATS := make([]genbank.GenbankFeature, 0)
	for _, F := range gb.FEATURES {
		if F.Feature == "CDS" {
			CDSFEATS = append(CDSFEATS, F)
		}
	}

	cdsregions := make([]Region, 0)

	// we get all the CDSes
	for _, feat := range CDSFEATS {
		REGION := Region{Whichtype: "CDS", Name: feat.Info["gene"], Codonstarts: make([]int, 0), Translation: feat.Info["translation"] + "*"}

		// these are genbank positions, so they are 1-based, inclusive
		positions, err := genbank.ParsePositions(feat.Pos)
		if err != nil {
			return make([]Region, 0), err
		}
		REGION.Start = positions[0]
		REGION.Stop = positions[len(positions)-1]

		// how many codons we have already defined if this CDS is two ranges Join()ed together:
		previouscodons := 0
		for i := 0; i < len(positions); i = i + 2 {
			start := positions[i]
			stop := positions[i+1]

			// start-1 to get the start position as 0-based
			if (stop-(start-1))%3 != 0 {
				return make([]Region, 0), errors.New("CDS position range is not a multiple of 3")
			}
			length := (stop - (start - 1))

			for j := 0; j < length; j = j + 3 {
				// we add another integer codon start to the slice
				REGION.Codonstarts = append(REGION.Codonstarts, start+j)
			}

			previouscodons = previouscodons + (length / 3)
		}
		cdsregions = append(cdsregions, REGION)
	}

	// then we add the intergenic regions (between the CDSes)
	regions := make([]Region, 0)
	newstart := 1
	for i, cdsregion := range cdsregions {
		start := newstart
		stop := cdsregion.Start - 1

		// hopefully this deals with any cases where there isn't an intergenic region:
		if !((stop - start) > 0) {
			continue
		}
		REGION := Region{Whichtype: "int", Start: start, Stop: stop}
		regions = append(regions, REGION)
		regions = append(regions, cdsregion)
		newstart = cdsregion.Stop + 1
		if i == len(cdsregions) {
			start := newstart
			stop := len(gb.ORIGIN)
			if !((stop - start) > 0) {
				continue
			}
			REGION := Region{Whichtype: "int", Start: start, Stop: stop}
			regions = append(regions, REGION)
		}
	}

	return regions, nil
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
	gapsum := 0
	refToMSA := make([]int, degappedLen)
	MSAToRef := make([]int, len(refseq))
	for i, nuc := range refseq {
		if nuc == 244 {
			gapsum++
			continue
		}
		refToMSA[i-gapsum] = gapsum
		MSAToRef[i] = gapsum
	}

	return refToMSA, MSAToRef
}

// getVariants annotates mutations between query and reference sequences, one fasta record at a time. It reads each fasta
// record from a channel and passes all its mutations to another channel grouped together in one struct.
func getVariants(ref fastaio.EncodedFastaRecord, regions []Region, offsetRefCoord []int, offsetMSACoord []int, cMSA chan fastaio.EncodedFastaRecord, cVariants chan AnnoStructs, cErr chan error) {

	for record := range cMSA {

		AS, err := GetVariantsPair(ref.Seq, record.Seq, ref.ID, record.ID, record.Idx, regions, offsetRefCoord, offsetMSACoord)
		if err != nil {
			cErr <- err
			break
		}

		cVariants <- AS
	}
}

// GetVariantsPair annotates mutations between one pair of query and reference sequences. ref and query are encoded according
// to EP's bitwise coding scheme
func GetVariantsPair(ref, query []byte, refID, queryID string, idx int, regions []Region, offsetRefCoord []int, offsetMSACoord []int) (AnnoStructs, error) {

	AS := AnnoStructs{}

	DA := encoding.MakeDecodingArray()
	CD := alphabet.MakeCodonDict()

	insOpen := false
	insStart := 0
	insLength := 0
	delOpen := false
	delStart := 0
	delLength := 0

	// check that the reference is the same length as this record
	// (might conceivably not be if the ref came from the genbank file and the msa has insertions in it)
	if len(query) != len(ref) {
		return AS, fmt.Errorf("sequence length for query %s (%d) is different to the sequence length of the reference %s (%d)", queryID, len(query), refID, len(ref))
	}

	// here is the slice of variants that we will populate, then sort, then put in an
	// annoStructs{} to write to file
	variants := make([]Variant, 0)

	// first, we get indels
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
					variants = append(variants, Variant{Changetype: "del", Position: delStart - offsetMSACoord[delStart], Length: delLength})
					delOpen = false // and reset things
				}
			}
		}
	}

	// catch things that abut the end of the alignment
	// don't want deletions at the end of the alignment (or at the beginning)
	// if delOpen {
	// 	variants = append(variants, Variant{Changetype: "del", Position: delStart - offsetMSACoord[delStart], Length: delLength})
	// }
	if insOpen {
		variants = append(variants, Variant{Changetype: "ins", Position: insStart - offsetMSACoord[insStart], Length: insLength})
	}

	// then we loop over the regions to get AAs and snps
	for _, region := range regions {
		// and switch on whether it is intergenic or CDS:
		switch region.Whichtype {
		case "int":
			adjStart := region.Start + offsetRefCoord[region.Start]
			adjStop := region.Stop + offsetRefCoord[region.Stop]
			for pos := adjStart - 1; pos < adjStop; pos++ {
				if (ref[pos] & query[pos]) < 16 { // check for SNPs
					variants = append(variants, Variant{Changetype: "nuc", RefAl: DA[ref[pos]], QueAl: DA[query[pos]], Position: pos - offsetMSACoord[pos]})
				}
			}
		case "CDS":
			codonSNPs := make([]Variant, 0, 3)
			decodedCodon := ""
			aa := ""
			refaa := ""
			// here are the 1-based start positions of each codon in reference coordinates
			for aaCounter, tempPos := range region.Codonstarts {
				// here is the actual position in the msa:
				pos := (tempPos - 1) + offsetRefCoord[tempPos]

				// for each position in this codon
				for codonCounter := 0; codonCounter < 3; codonCounter++ {
					// skip insertions relative to the reference
					// TO DO = if they are in this record log/take them into account
					if ref[pos+codonCounter] == 244 {
						continue
					}
					if (query[pos+codonCounter] & ref[pos+codonCounter]) < 16 {
						codonSNPs = append(codonSNPs, Variant{Changetype: "nuc", RefAl: DA[ref[pos+codonCounter]], QueAl: DA[query[pos+codonCounter]], Position: (pos + codonCounter) - offsetMSACoord[pos+codonCounter]})
					}
					decodedCodon = decodedCodon + DA[query[pos+codonCounter]]
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
					variants = append(variants, Variant{Changetype: "aa", Feature: region.Name, RefAl: refaa, QueAl: aa, Position: pos - offsetMSACoord[pos], Residue: aaCounter, SNPs: strings.Join(temp, ";")})

				} else {
					for _, v := range codonSNPs {
						variants = append(variants, v)
					}
				}

				codonSNPs = make([]Variant, 0, 3)
				decodedCodon = ""
			}
		}
	}

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

// FormatVariant returns a string representation of a single mutation the format of which varies
// given its type (aa/nuc/indel)
func FormatVariant(v Variant, appendSNP bool) (string, error) {
	var s string

	switch v.Changetype {
	case "del":
		s = "del:" + strconv.Itoa(v.Position+1) + ":" + strconv.Itoa(v.Length)
	case "ins":
		s = "ins:" + strconv.Itoa(v.Position) + ":" + strconv.Itoa(v.Length)
	case "nuc":
		s = "nuc:" + v.RefAl + strconv.Itoa(v.Position+1) + v.QueAl
	case "aa":
		if appendSNP {
			s = "aa:" + v.Feature + ":" + v.RefAl + strconv.Itoa(v.Residue+1) + v.QueAl + "(" + v.SNPs + ")"
		} else {
			s = "aa:" + v.Feature + ":" + v.RefAl + strconv.Itoa(v.Residue+1) + v.QueAl
		}
	default:
		return "", errors.New("couldn't parse variant type")
	}

	return s, nil
}

// WriteVariants writes each queries mutations to file or stdout
func WriteVariants(w io.Writer, appendSNP bool, refID string, cVariants chan AnnoStructs, cWriteDone chan bool, cErr chan error) {

	outputMap := make(map[int]AnnoStructs)

	counter := 0

	var err error
	var sa []string

	_, err = w.Write([]byte("query,mutations\n"))
	if err != nil {
		cErr <- err
		return
	}

	for variantLine := range cVariants {
		outputMap[variantLine.Idx] = variantLine

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
				newVar, err := FormatVariant(v, appendSNP)
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
			continue
		}
	}

	for n := 1; n > 0; {
		if len(outputMap) == 0 {
			n--
			break
		}

		VL := outputMap[counter]
		if VL.Queryname == "ref" {
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
			newVar, err := FormatVariant(v, appendSNP)
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
	}

	cWriteDone <- true
}

// AggregateWriteOutput aggregates the mutations that are present above a certain threshold, and
// writes their frequencies out to file or stdout
func AggregateWriteVariants(w io.Writer, appendSNP bool, threshold float64, refID string, cVariants chan AnnoStructs, cWriteDone chan bool, cErr chan error) {

	propMap := make(map[Variant]float64)

	var err error

	_, err = w.Write([]byte("mutation,frequency\n"))
	if err != nil {
		cErr <- err
		return
	}

	counter := 0.0

	for AS := range cVariants {
		counter++
		for _, V := range AS.Vs {
			Vskinny := Variant{RefAl: V.RefAl, QueAl: V.QueAl, Position: V.Position, Residue: V.Residue, Changetype: V.Changetype, Feature: V.Feature, Length: V.Length, SNPs: V.SNPs}
			if _, ok := propMap[Vskinny]; ok {
				propMap[Vskinny]++
			} else {
				propMap[Vskinny] = 1.0
			}
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
		Vformatted, err := FormatVariant(V, appendSNP)
		if err != nil {
			cErr <- err
			return
		}
		_, err = w.Write([]byte(Vformatted + "," + strconv.FormatFloat(propMap[V]/counter, 'f', 9, 64) + "\n"))
		if err != nil {
			cErr <- err
			return
		}
	}

	cWriteDone <- true
}
