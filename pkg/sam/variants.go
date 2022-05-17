package sam

import (
	"errors"
	"io"
	"os"
	"sync"

	biogosam "github.com/biogo/hts/sam"
	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fastaio"
	"github.com/virus-evolution/gofasta/pkg/genbank"
	"github.com/virus-evolution/gofasta/pkg/gff"
	"github.com/virus-evolution/gofasta/pkg/variants"
)

// Variants annotates amino acid, insertion, deletion, and nucleotide (anything outside of codons with an amino acid change)
// mutations relative to a reference sequence from pairwise alignments in sam format. Genome annotations are derived from a annotation file
// in genbank or gff version 3 format
func Variants(samIn, refIn io.Reader, refFromFile bool, annoIn io.Reader, annoSuffix string, out io.Writer, aggregate bool, threshold float64, appendSNP bool, threads int) error {

	var ref fastaio.EncodedFastaRecord
	if refFromFile {
		refs, err := fastaio.ReadEncodeAlignmentToList(refIn, false)
		if err != nil {
			return err
		}
		if len(refs) > 1 {
			return errors.New("more than one record in --reference")
		}
		ref = refs[0]
	}

	var regions []variants.Region
	switch annoSuffix {
	case "gb":
		gb, err := genbank.ReadGenBank(annoIn)
		if err != nil {
			return err
		}
		regions, err = variants.GetRegionsFromGenbank(gb)
		if err != nil {
			return err
		}
		if !refFromFile {
			temp := fastaio.FastaRecord{Seq: string(gb.ORIGIN), ID: "annotation_fasta"}
			ref = temp.Encode()
			os.Stderr.WriteString("using --annotation fasta as reference\n")
		}
	case "gff":
		gff, err := gff.ReadGFF(annoIn)
		if err != nil {
			return err
		}
		if !refFromFile {
			switch len(gff.FASTA) {
			case 0:
				return errors.New("couldn't find a reference sequence in the gff and none was provided to --reference")
			case 1:
				var encodedrefseq []byte
				for _, v := range gff.FASTA {
					encodedrefseq = make([]byte, len(v.Seq))
					EA := encoding.MakeEncodingArray()
					for i := range v.Seq {
						encodedrefseq[i] = EA[v.Seq[i]]
					}
					ref = fastaio.EncodedFastaRecord{ID: "annotation_fasta", Seq: encodedrefseq}
				}
				os.Stderr.WriteString("using --annotation fasta as reference\n")
			default:
				return errors.New("more that one sequence in gff ##FASTA section")
			}
		}
		regions, err = variants.GetRegionsFromGFF(gff, ref.Decode().Seq)
		if err != nil {
			return err
		}
	}

	cErr := make(chan error)

	// do some things that are basically just sam topairalign:
	cSR := make(chan samRecords, threads)
	cSH := make(chan biogosam.Header)
	cPairAlign := make(chan alignPair)

	cVariants := make(chan variants.AnnoStructs)

	cReadDone := make(chan bool)
	cAlignWaitGroupDone := make(chan bool)
	cVariantsDone := make(chan bool)
	cWriteDone := make(chan bool)

	switch aggregate {
	case true:
		go variants.AggregateWriteVariants(out, appendSNP, threshold, ref.ID, cVariants, cWriteDone, cErr)
	case false:
		go variants.WriteVariants(out, false, appendSNP, ref.ID, cVariants, cWriteDone, cErr)
	}

	go groupSamRecords(samIn, cSH, cSR, cReadDone, cErr)

	_ = <-cSH

	var wgAlign sync.WaitGroup
	wgAlign.Add(threads)

	var wgVariants sync.WaitGroup
	wgVariants.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			blockToPairwiseAlignment(cSR, cPairAlign, cErr, []byte(ref.Decode().Seq), false)
			wgAlign.Done()
		}()
	}

	for n := 0; n < threads; n++ {
		go func() {
			getVariantsSam(regions, cPairAlign, cVariants, cErr)
			wgVariants.Done()
		}()
	}

	go func() {
		wgAlign.Wait()
		cAlignWaitGroupDone <- true
	}()

	go func() {
		wgVariants.Wait()
		cVariantsDone <- true
	}()

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cReadDone:
			close(cSR)
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

// getVariantsSam gets the mutations for each pairwise alignment from a channel at a time, and passes them to a channel
// of annotated variants, given an array of annotated genome regions
func getVariantsSam(regions []variants.Region, cAlignPair chan alignPair, cVariants chan variants.AnnoStructs, cErr chan error) {

	EA := encoding.MakeEncodingArray()

	for pair := range cAlignPair {

		for i, nuc := range pair.query {
			pair.query[i] = EA[nuc]
		}

		for i, nuc := range pair.ref {
			pair.ref[i] = EA[nuc]
		}

		offsetRefCoord, offsetMSACoord := variants.GetMSAOffsets(pair.ref)

		AS, err := variants.GetVariantsPair(pair.ref, pair.query, pair.refname, pair.queryname, pair.idx, regions, offsetRefCoord, offsetMSACoord)
		if err != nil {
			cErr <- err
			break
		}

		// and we're done
		cVariants <- AS
	}
}
