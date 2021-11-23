package sam

import (
	"io"
	"sync"

	biogosam "github.com/biogo/hts/sam"
	"github.com/cov-ert/gofasta/pkg/encoding"
	"github.com/cov-ert/gofasta/pkg/fastaio"
	"github.com/cov-ert/gofasta/pkg/genbank"
	"github.com/cov-ert/gofasta/pkg/variants"
)

// Variants annotates amino acid, insertion, deletion, and nucleotide (anything outside of codons represented by an amino acid change)
// mutations relative to a reference sequence from pairwise alignments in sam format. Genome annotations are derived from a genbank flat file
func Variants(samIn, ref, genbankIn io.Reader, out io.Writer, threads int) error {

	cErr := make(chan error)
	cRef := make(chan fastaio.FastaRecord)
	cRefDone := make(chan bool)

	go fastaio.ReadAlignment(ref, cRef, cErr, cRefDone)

	var refSeq string
	var refID string

	// to do - do this more sensibly
	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case FR := <-cRef:
			refSeq = FR.Seq
			refID = FR.ID
			// refName = FR.ID
		case <-cRefDone:
			close(cRef)
			n--
		}
	}

	gb, err := genbank.ReadGenBank(genbankIn)
	if err != nil {
		return err
	}

	// get a list of CDS + intergenic regions from the genbank file
	regions, err := variants.GetRegions(gb)
	if err != nil {
		return err
	}

	// do some things that are basically just sam topairalign:
	cSR := make(chan samRecords, threads)
	cSH := make(chan biogosam.Header)
	cPairAlign := make(chan alignPair)

	cVariants := make(chan variants.AnnoStructs)

	cReadDone := make(chan bool)
	cAlignWaitGroupDone := make(chan bool)
	cVariantsDone := make(chan bool)
	cWriteDone := make(chan bool)

	go variants.WriteVariants(out, refID, cVariants, cWriteDone, cErr)

	go groupSamRecords(samIn, cSH, cSR, cReadDone, cErr)

	_ = <-cSH

	var wgAlign sync.WaitGroup
	wgAlign.Add(threads)

	var wgVariants sync.WaitGroup
	wgVariants.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			blockToPairwiseAlignment(cSR, cPairAlign, cErr, []byte(refSeq), false)
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
