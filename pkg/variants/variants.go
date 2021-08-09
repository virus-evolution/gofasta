package variants

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/spf13/cobra"

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

//
type variant struct {
	queryname  string
	refAl      string
	queAl      string
	position   int // genomic location
	residue    int // feature location?
	changetype string
	feature    string // this should be, for example, the name of the CDS that the thing is in
	length     int    // for indels
}

// for passing groups of annoStruct around with an index which is used to retain input
// order in the output
type annoStructs struct {
	queryname string
	vs        []variant
	idx       int
}

var msa string
var gbfile string
var reference string
var outfile string

var mainCmd = &cobra.Command{
	Use:   "variants",
	Short: "",
	Long: `variants ...

Example usage:

./variants ...
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		err = Variants(msa, reference, gbfile, outfile)

		return
	},
}

func init() {

	mainCmd.Flags().StringVarP(&msa, "msa", "", "", "mutliple sequence alignment in fasta format")
	mainCmd.Flags().StringVarP(&gbfile, "genbank", "", "", "genbank format annotation")
	mainCmd.Flags().StringVarP(&reference, "reference", "", "", "the name of the reference sequence in the msa")
	mainCmd.Flags().StringVarP(&outfile, "outfile", "o", "stdout", "the name of the reference sequence in the msa")

	mainCmd.Flags().SortFlags = false
}

func main() {
	mainCmd.Execute()
}

func Variants(msa string, reference string, gbfile string, outfile string) error {

	gb, err := genbank.ReadGenBank(gbfile)
	if err != nil {
		return err
	}

	ref, err := findReference(msa, reference)
	if err != nil {
		return err
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
	regions, err := getRegions(gbfile)
	if err != nil {
		return err
	}

	// get the offset accounting for insertions relative to the reference
	offsetRefCoord, offsetMSACoord := getOffsets(ref.Seq)

	cMSA := make(chan fastaio.EncodedFastaRecord)
	cErr := make(chan error)
	cMSADone := make(chan bool)

	cVariants := make(chan annoStructs)
	cVariantsDone := make(chan bool)
	cWriteDone := make(chan bool)

	go writeVariants(outfile, cVariants, cWriteDone, cErr)

	go fastaio.ReadEncodeAlignment(msa, cMSA, cErr, cMSADone)

	var wgVariants sync.WaitGroup
	wgVariants.Add(runtime.NumCPU())

	for n := 0; n < runtime.NumCPU(); n++ {
		go func() {
			getVariants(ref, regions, offsetRefCoord, offsetMSACoord, cMSA, cVariants, cErr)
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

// get the reference sequence from the msa if it is in there. If it isn't, we will try get it
// from the genbank record (in which case can be no insertions relative to the ref in the msa)
func findReference(msa string, reference string) (fastaio.EncodedFastaRecord, error) {

	var err error
	var f *os.File

	if msa != "stdin" {
		f, err = os.Open(msa)
		if err != nil {
			return fastaio.EncodedFastaRecord{}, err
		}
	} else {
		f = os.Stdin
	}

	defer f.Close()

	coding := encoding.MakeEncodingArray()

	s := bufio.NewScanner(f)

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

			if id == reference {
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

			if id == reference {
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

	err = s.Err()
	if err != nil {
		return fastaio.EncodedFastaRecord{}, err
	}

	if refFound {
		refRec = fastaio.EncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
	} else {
		fmt.Println("no reference sequence found in --msa, attempting to use --genbank source as reference")
	}

	return refRec, nil
}

func getRegions(genbankFileIn string) ([]Region, error) {
	gb, err := genbank.ReadGenBank(genbankFileIn)
	if err != nil {
		return make([]Region, 0), err
	}

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
		positions, err := parsePositions(feat.Pos)
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

func parsePositions(position string) ([]int, error) {
	var A []int
	if position[0:4] == "join" {
		A = make([]int, 0)
		position = strings.TrimLeft(position, "join(")
		position = strings.TrimRight(position, ")")
		ranges := strings.Split(position, ",")
		for _, x := range ranges {
			y := strings.Split(x, "..")
			for _, z := range y {
				temp, err := strconv.Atoi(z)
				if err != nil {
					return []int{}, err
				}
				A = append(A, temp)
			}
		}
	} else {
		A = make([]int, 0)
		y := strings.Split(position, "..")
		for _, z := range y {
			temp, err := strconv.Atoi(z)
			if err != nil {
				return []int{}, err
			}
			A = append(A, temp)
		}
	}

	return A, nil
}

func getOffsets(refseq []byte) ([]int, []int) {
	// we make an array of integers to offset the positions by.
	// this should be the same length as a degapped refseq?
	degappedLen := 0
	for _, nuc := range refseq {
		// if there is no alignment gap at this site, ++
		if nuc != 244 {
			degappedLen++
		}
	}
	// if there are no alignment gaps, we can return two slices of 0s
	if degappedLen == len(refseq) {
		return make([]int, len(refseq), len(refseq)), make([]int, len(refseq), len(refseq))
	}

	// otherwise, we get some offsets:
	// 1) offsetRefCoord = the number of bases to add to convert each position to msa coordinates
	// 2) offsetMSACoord = the number of bases to subtract to convert each position to reference coordinates
	gapsum := 0
	offsetRefCoord := make([]int, degappedLen)
	offsetMSACoord := make([]int, len(refseq))
	for i, nuc := range refseq {
		if nuc == 244 {
			gapsum++
			continue
		}
		offsetRefCoord[i-gapsum] = gapsum
		offsetMSACoord[i] = gapsum
	}

	return offsetRefCoord, offsetMSACoord
}

func getVariants(ref fastaio.EncodedFastaRecord, regions []Region, offsetRefCoord []int, offsetMSACoord []int, cMSA chan fastaio.EncodedFastaRecord, cVariants chan annoStructs, cErr chan error) {

	DA := encoding.MakeDecodingArray()
	CD := alphabet.MakeCodonDict()
	for record := range cMSA {

		insOpen := false
		insStart := 0
		insLength := 0
		delOpen := false
		delStart := 0
		delLength := 0

		// check that the reference is the same length as this record
		// (might conceivably not be if the ref came from the genbank file and the msa has insertions in it)
		if len(record.Seq) != len(ref.Seq) {
			cErr <- fmt.Errorf("sequence length for query %s (%d) is different to the sequence length of the reference %s (%d)", record.ID, len(record.Seq), ref.ID, len(ref.Seq))
			break
		}

		// here is the slice of variants that we will populate, then sort, then put in an
		// annoStructs{} to write to file
		variants := make([]variant, 0)

		// first, we get indels
		for pos := range record.Seq {
			if ref.Seq[pos] == 244 { // insertion relative to reference (somewhere in the alignment)
				if record.Seq[pos] == 244 { // insertion is not in this seq
					continue
				} else { // insertion is in this seq
					if insOpen { // not the first position of an insertion
						insLength++ // we increment the length counter
					} else { // the first position of an insertion
						insStart = pos - 1 // we record the position of the insertion as the left-neighbouring reference position
						insLength = 1
						insOpen = true
					}
				}
			} else { // not an insertion relative to the reference at this position
				if insOpen { // first base after an insertion, so we need to log the insertion
					variants = append(variants, variant{changetype: "ins", position: insStart - offsetMSACoord[insStart], length: insLength})
					insOpen = false
				}
				if record.Seq[pos] == 244 { // deletion in this seq
					if delOpen { // not the first position of a deletion
						delLength++ // we increment the length (there is not a deletion in the reference)
					} else { // the first position of a deletion
						delStart = pos
						delLength = 1
						delOpen = true
					}
				} else { // no deletion in this seq
					if delOpen { // first base after a deletion, so we need to log the deletion
						variants = append(variants, variant{changetype: "del", position: delStart - offsetMSACoord[delStart], length: delLength})
						delOpen = false // and reset things
					}
				}
			}
		}

		// catch things that abut the end of the alignment
		if delOpen {
			variants = append(variants, variant{changetype: "del", position: delStart - offsetMSACoord[delStart], length: delLength})
		}
		if insOpen {
			variants = append(variants, variant{changetype: "ins", position: insStart - offsetMSACoord[insStart], length: insLength})
		}

		// then we loop over the regions to get AAs and snps
		for _, region := range regions {
			// and switch on whether it is intergenic or CDS:
			switch region.Whichtype {
			case "int":
				adjStart := region.Start + offsetRefCoord[region.Start]
				adjStop := region.Stop + offsetRefCoord[region.Stop]
				for pos := adjStart - 1; pos < adjStop; pos++ {
					if (ref.Seq[pos] & record.Seq[pos]) < 16 { // check for SNPs
						variants = append(variants, variant{changetype: "nuc", refAl: DA[ref.Seq[pos]], queAl: DA[record.Seq[pos]], position: pos - offsetMSACoord[pos]})
					}
				}
			case "CDS":
				codonSNPs := make([]variant, 0, 3)
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
						if ref.Seq[pos+codonCounter] == 244 {
							continue
						}
						if (record.Seq[pos+codonCounter] & ref.Seq[pos+codonCounter]) < 16 {
							codonSNPs = append(codonSNPs, variant{changetype: "nuc", refAl: DA[ref.Seq[pos+codonCounter]], queAl: DA[record.Seq[pos+codonCounter]], position: (pos + codonCounter) - offsetMSACoord[pos+codonCounter]})
						}
						decodedCodon = decodedCodon + DA[record.Seq[pos+codonCounter]]
					}

					if _, ok := CD[decodedCodon]; ok {
						aa = CD[decodedCodon]
					} else {
						aa = "X"
					}

					refaa = string(region.Translation[aaCounter])

					if aa != refaa && aa != "X" {
						variants = append(variants, variant{changetype: "aa", feature: region.Name, refAl: refaa, queAl: aa, position: pos - offsetMSACoord[pos], residue: aaCounter})
					} else {
						for _, v := range codonSNPs {
							variants = append(variants, v)
						}
					}

					codonSNPs = make([]variant, 0, 3)
					decodedCodon = ""
				}
			}
		}

		// sort the variants
		sort.SliceStable(variants, func(i, j int) bool {
			return variants[i].position < variants[j].position || (variants[i].position == variants[j].position && variants[i].changetype < variants[j].changetype)
		})

		// there might be dups if there was a snp in the region of a join()
		finalVariants := make([]variant, 0)
		previousVariant := variant{}
		for i, v := range variants {
			if i == 0 {
				finalVariants = append(finalVariants, v)
				previousVariant = v
				continue
			}
			if v == previousVariant {
				continue
			}
			finalVariants = append(finalVariants, v)
			previousVariant = v
		}

		AS := annoStructs{queryname: record.ID, vs: finalVariants, idx: record.Idx}

		// if this is the reference, let the writer know
		if record.ID == ref.ID {
			AS.queryname = "ref"
		}

		// and we're done
		cVariants <- AS
	}
}

func formatVariant(v variant, cErr chan error) string {
	var s string

	switch v.changetype {
	case "del":
		s = "del:" + strconv.Itoa(v.position+1) + ":" + strconv.Itoa(v.length)
	case "ins":
		s = "ins:" + strconv.Itoa(v.position+1) + ":" + strconv.Itoa(v.length)
	case "nuc":
		s = "nuc:" + v.refAl + strconv.Itoa(v.position+1) + v.queAl
	case "aa":
		s = "aa:" + v.feature + ":" + v.refAl + strconv.Itoa(v.residue+1) + v.queAl
	default:
		cErr <- errors.New("couldn't parse variant type")
	}

	return s
}

func writeVariants(outfile string, cVariants chan annoStructs, cWriteDone chan bool, cErr chan error) {

	outputMap := make(map[int]annoStructs)

	counter := 0

	var f *os.File
	var err error
	var sa []string

	if outfile != "stdout" {
		f, err = os.Create(outfile)
		if err != nil {
			cErr <- err
		}
	} else {
		f = os.Stdout
	}

	defer f.Close()

	_, err = f.WriteString("query,variants\n")
	if err != nil {
		cErr <- err
	}

	for variantLine := range cVariants {
		outputMap[variantLine.idx] = variantLine

		if VL, ok := outputMap[counter]; ok {

			if VL.queryname == "ref" {
				delete(outputMap, counter)
				counter++
				continue
			}

			_, err = f.WriteString(VL.queryname + ",")
			if err != nil {
				cErr <- err
			}
			sa = make([]string, 0)
			for _, v := range VL.vs {
				sa = append(sa, formatVariant(v, cErr))
			}
			_, err = f.WriteString(strings.Join(sa, "|") + "\n")
			if err != nil {
				cErr <- err
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
		if VL.queryname == "ref" {
			delete(outputMap, counter)
			counter++
			continue
		}
		_, err = f.WriteString(VL.queryname + ",")
		if err != nil {
			cErr <- err
		}
		sa = make([]string, 0)
		for _, v := range VL.vs {
			sa = append(sa, formatVariant(v, cErr))
		}
		_, err = f.WriteString(strings.Join(sa, "|") + "\n")
		if err != nil {
			cErr <- err
		}

		delete(outputMap, counter)
		counter++
	}

	cWriteDone <- true
}
