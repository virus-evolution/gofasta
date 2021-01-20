package sam

import (
	"errors"
	"fmt"
	"sync"
	// "sort"
	"os"
	// "path"
	"strings"
	"strconv"

	"github.com/cov-ert/gofasta/pkg/genbank"
	"github.com/cov-ert/gofasta/pkg/fastaio"
	"github.com/cov-ert/gofasta/pkg/alphabet"
	"github.com/cov-ert/gofasta/pkg/encoding"

	biogosam "github.com/biogo/hts/sam"
)

//
type annoStruct struct {
	queryname string
	refAl string
	queAl string
	position int
	changetype string
	feature string // this should be, for example, the name of the CDS that the thing is in
}

// type alignPair struct {
// 	ref []byte
// 	query []byte
// 	refname string
// 	queryname string
// 	descriptor string
// }

// getSNPPos gets the position in reference coordinates from a subalignment
// of some particular feature
// TODO: use sam.getRefAdjustedPositions & sam.findOffsetPos to get the coordinates
// when insertions are included
func getSNPPos(i int, A []int) (int, error) {

	pos := -1
	leftadd := 0

	if (len(A) / 2) > 1 {

		// e.g. join(266..13468,13468..21555)
		// [266, 13468, 13468, 21555]

		for j := 1; j < len(A); j += 2 {
			if i <= (A[j] - (A[j-1] - leftadd)) {
				pos = A[j-1] - leftadd + i
				break
			}
			leftadd += (A[j] - A[j-1] + 1)
		}

	} else if len(A) == 2 {
		pos = A[0] + i
	} else {
		s := []string{}
		for _, j := range(A) {
			s = append(s, strconv.Itoa(j))
		}
		return 0, fmt.Errorf("bad(ly parsed?) annotation positions array: %s", strings.Join(s, " "))
	}

	if pos < 0 {
		return 0, errors.New("badly parsed SNP position")
	}

	return pos, nil
}

func getVariantsFromAlignPair(pair alignPair) ([]annoStruct, error) {

	codon_2_AA := alphabet.MakeCodonDict()
	rune_2_byte := encoding.MakeByteDict2()

	if len(pair.ref) != len(pair.query) {
		return []annoStruct{}, errors.New("ref and query sequences are different lengths")
	}

	annotation_array := make([]annoStruct, 0)

	ref_codon := make([]byte, 3)
	que_codon := make([]byte, 3)

	codon_snps := make([]annoStruct, 0)

	counter := 0

	// now compare ref and query
	for i, _ := range(pair.ref) {

		if pair.ref[i] != pair.query[i] {
			a := rune_2_byte[pair.ref[i]]
			b := rune_2_byte[pair.query[i]]

			if (a & b) < 16 {
				// need to calculate pos here because pair.ref is no longer in reference coordinates if the genbank feature is two stretches of sequence Join()ed together:
				pos, err := getSNPPos(i, pair.featPosArray)
				if err != nil {
					return []annoStruct{}, err
				}
				codon_snps = append(codon_snps, annoStruct{queryname: pair.queryname, refAl: string(pair.ref[i]), queAl: string(pair.query[i]), position: pos, changetype: "SNP", feature: pair.featName})
			}
		}

		ref_codon[counter] = pair.ref[i]
		que_codon[counter] = pair.query[i]

		counter += 1

		if counter == 3 {
			ref_AA, ok_ref := codon_2_AA[string(ref_codon)]
			que_AA, ok_que := codon_2_AA[string(que_codon)]

			if ok_ref && ok_que {

				if ref_AA != que_AA {
					annotation_array = append(annotation_array, annoStruct{queryname: pair.queryname, refAl: ref_AA, queAl: que_AA, position: (i + 1) / 3, changetype: "AA", feature: pair.featName})
				} else {
					if len(codon_snps) > 0 {
						for _, snp := range(codon_snps) {
							snp.changetype = "synSNP"
							annotation_array = append(annotation_array, snp)
						}
					}
				}
			}

			codon_snps = make([]annoStruct, 0)
			counter = 0
		}

	}

	return annotation_array, nil
}

// Apply some other function over the channel of align pairs
func getVariantsFromCDS(cPairParse chan []alignPair, cAnnotate chan []annoStruct, cErr chan error) {
	// this is what comes with the descriptor field of each alignPair struct from cPairParse:
	// subPair.descriptor = pair.queryname + "." + feature.Feature + "." + strings.ReplaceAll(feature.Info[anno], " ", "_")

	for A := range(cPairParse) {

		annoArray := make([]annoStruct, 0)

		for _, pair := range(A) {
			anno, err := getVariantsFromAlignPair(pair)
			if err != nil {
				cErr<- err
			}

			annoArray = append(annoArray, anno...)
		}

		cAnnotate<- annoArray
	}
}

func getAnnoLine(aS annoStruct) (string, error) {

	if aS.changetype == "synSNP" {
		s := aS.changetype + ":" + aS.refAl + strconv.Itoa(aS.position) + aS.queAl
		return s, nil
	}
	if aS.changetype == "AA" {
		s := aS.feature + ":" + aS.refAl + strconv.Itoa(aS.position) + aS.queAl
		return s, nil
	}

	return "", errors.New("couldn't parse variant for writing out; unrecognised variant type: needs to be one of AA or synSNP")
}

// write the annotation
func writeAnnotation(outfile string, cAnnotate chan []annoStruct, cWriteDone chan bool, cErr chan error) {

	_ = strings.Join

	var err error
	f := os.Stdout

	if outfile != "stdout" {
		f, err = os.Create(outfile)
		if err != nil {
			cErr<- err
		}
	}

	defer f.Close()

	_, err = f.WriteString("query,variants\n")
	if err != nil {
		cErr<- err
	}

	for A := range(cAnnotate) {

		if len(A) == 0 {
			continue
		}

		queryname := A[0].queryname

		f.WriteString(queryname + ",")

		temp := make([]string, 0)

		for _, aS := range(A) {
			AL, err := getAnnoLine(aS)
			if err != nil {
				cErr<- err
			}
			temp = append(temp, AL)
		}

		s := strings.Join(temp, "|")

		f.WriteString(s)
		f.WriteString("\n")
	}

	cWriteDone<- true
}

// Variants annotates variants wrt. a reference sequence
func Variants(samFile string, referenceFile string, genbankFile string,
	      outfile string, threads int) error {

	gb, err := genbank.ReadGenBank(genbankFile)
	if err != nil {
		return err
	}

	cErr := make(chan error)

	cRef := make(chan fastaio.FastaRecord)
	cRefDone := make(chan bool)

	go fastaio.ReadAlignment(referenceFile, cRef, cErr, cRefDone)

	var refSeq string

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case FR := <-cRef:
			refSeq = FR.Seq
			// refName = FR.ID
		case <-cRefDone:
			close(cRef)
			n--
		}
	}

	cSamRecords := make(chan samRecords, threads)
	cSH := make(chan biogosam.Header)
	cPairAlign := make(chan alignPair)
	cPairParse := make(chan []alignPair)
	cVariants := make(chan []annoStruct)

	cReadDone := make(chan bool)
	cAlignWaitGroupDone := make(chan bool)
	cParseWaitGroupDone := make(chan bool)
	cVariantsDone := make(chan bool)
	cWriteDone := make(chan bool)

	go groupSamRecords(samFile, cSH, cSamRecords, cReadDone, cErr)

	go writeAnnotation(outfile, cVariants, cWriteDone, cErr)

	var wgAlign sync.WaitGroup
	wgAlign.Add(threads)

	var wgParse sync.WaitGroup
	wgParse.Add(threads)

	var wgVar sync.WaitGroup
	wgVar.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			blockToPairwiseAlignment(cSamRecords, cPairAlign, cErr, []byte(refSeq), true)
			wgAlign.Done()
		}()
	}

	features := getFeaturesFromAnnotation(gb, "CDS")

	for n := 0; n < threads; n++ {
		go func() {
			parseAlignmentByAnnotation(features, cPairAlign, cPairParse, cErr)
			wgParse.Done()
		}()
	}

	for n := 0; n < threads; n++ {
		go func() {
			getVariantsFromCDS(cPairParse, cVariants, cErr)
			wgVar.Done()
		}()
	}

	go func() {
		wgAlign.Wait()
		cAlignWaitGroupDone<- true
	}()

	go func() {
		wgParse.Wait()
		cParseWaitGroupDone<- true
	}()

	go func() {
		wgVar.Wait()
		cVariantsDone<- true
	}()

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cReadDone:
			close(cSamRecords)
			close(cSH)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cAlignWaitGroupDone:
			close(cPairAlign)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cParseWaitGroupDone:
			close(cPairParse)
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
