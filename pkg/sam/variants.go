package sam

import (
	"fmt"
	"io"
	"sort"
	"sync"

	biogosam "github.com/biogo/hts/sam"
	"github.com/virus-evolution/gofasta/pkg/alphabet"
	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fastaio"
	"github.com/virus-evolution/gofasta/pkg/variants"
)

func Variants(samIn, ref, gb io.Reader, out io.Writer, threads int) error {

	cErr := make(chan error)
	cRef := make(chan fastaio.FastaRecord)
	cRefDone := make(chan bool)

	go fastaio.ReadAlignment(ref, cRef, cErr, cRefDone)

	var refSeq string

	// to do - do this more sensibly
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

	go variants.WriteVariants(out, cVariants, cWriteDone, cErr)

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

func getVariantsSam(regions []variants.Region, cAlignPair chan alignPair, cVariants chan variants.AnnoStructs, cErr chan error) {

	EA := encoding.MakeEncodingArray()
	DA := encoding.MakeDecodingArray()
	CD := alphabet.MakeCodonDict()

	for pair := range cAlignPair {

		for i, nuc := range pair.query {
			pair.query[i] = EA[nuc]
		}

		for i, nuc := range pair.ref {
			pair.ref[i] = EA[nuc]
		}

		offsetRefCoord, offsetMSACoord := variants.GetOffsets(pair.ref)

		insOpen := false
		insStart := 0
		insLength := 0
		delOpen := false
		delStart := 0
		delLength := 0

		// check that the reference is the same length as this record
		// (might conceivably not be if the ref came from the genbank file and the msa has insertions in it)
		if len(pair.ref) != len(pair.query) {
			cErr <- fmt.Errorf("sequence length for query %s (%d) is different to the sequence length of the reference %s (%d)", pair.queryname, len(pair.query), pair.refname, len(pair.ref))
			break
		}

		// here is the slice of variants that we will populate, then sort, then put in an
		// annoStructs{} to write to file
		vs := make([]variants.Variant, 0)

		// first, we get indels
		for pos := range pair.query {
			if pair.ref[pos] == 244 { // insertion relative to reference (somewhere in the alignment)
				if pair.query[pos] == 244 { // insertion is not in this seq
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
					vs = append(vs, variants.Variant{Changetype: "ins", Position: insStart - offsetMSACoord[insStart], Length: insLength})
					insOpen = false
				}
				if pair.query[pos] == 244 { // deletion in this seq
					if delOpen { // not the first position of a deletion
						delLength++ // we increment the length (there is not a deletion in the reference)
					} else { // the first position of a deletion
						delStart = pos
						delLength = 1
						delOpen = true
					}
				} else { // no deletion in this seq
					if delOpen { // first base after a deletion, so we need to log the deletion
						vs = append(vs, variants.Variant{Changetype: "del", Position: delStart - offsetMSACoord[delStart], Length: delLength})
						delOpen = false // and reset things
					}
				}
			}
		}

		// catch things that abut the end of the alignment
		if delOpen {
			vs = append(vs, variants.Variant{Changetype: "del", Position: delStart - offsetMSACoord[delStart], Length: delLength})
		}
		if insOpen {
			vs = append(vs, variants.Variant{Changetype: "ins", Position: insStart - offsetMSACoord[insStart], Length: insLength})
		}

		// then we loop over the regions to get AAs and snps
		for _, region := range regions {
			// and switch on whether it is intergenic or CDS:
			switch region.Whichtype {
			case "int":
				adjStart := region.Start + offsetRefCoord[region.Start]
				adjStop := region.Stop + offsetRefCoord[region.Stop]
				for pos := adjStart - 1; pos < adjStop; pos++ {
					if (pair.ref[pos] & pair.query[pos]) < 16 { // check for SNPs
						vs = append(vs, variants.Variant{Changetype: "nuc", RefAl: DA[pair.ref[pos]], QueAl: DA[pair.query[pos]], Position: pos - offsetMSACoord[pos]})
					}
				}
			case "CDS":
				codonSNPs := make([]variants.Variant, 0, 3)
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
						if pair.ref[pos+codonCounter] == 244 {
							continue
						}
						if (pair.query[pos+codonCounter] & pair.ref[pos+codonCounter]) < 16 {
							codonSNPs = append(codonSNPs, variants.Variant{Changetype: "nuc", RefAl: DA[pair.ref[pos+codonCounter]], QueAl: DA[pair.query[pos+codonCounter]], Position: (pos + codonCounter) - offsetMSACoord[pos+codonCounter]})
						}
						decodedCodon = decodedCodon + DA[pair.query[pos+codonCounter]]
					}

					if _, ok := CD[decodedCodon]; ok {
						aa = CD[decodedCodon]
					} else {
						aa = "X"
					}

					refaa = string(region.Translation[aaCounter])

					if aa != refaa && aa != "X" {
						vs = append(vs, variants.Variant{Changetype: "aa", Feature: region.Name, RefAl: refaa, QueAl: aa, Position: pos - offsetMSACoord[pos], Residue: aaCounter})
					} else {
						for _, v := range codonSNPs {
							vs = append(vs, v)
						}
					}

					codonSNPs = make([]variants.Variant, 0, 3)
					decodedCodon = ""
				}
			}
		}

		// sort the variants
		sort.SliceStable(vs, func(i, j int) bool {
			return vs[i].Position < vs[j].Position || (vs[i].Position == vs[j].Position && vs[i].Changetype < vs[j].Changetype)
		})

		// there might be dups if there was a snp in the region of a join()
		finalVariants := make([]variants.Variant, 0)
		previousVariant := variants.Variant{}
		for i, v := range vs {
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

		AS := variants.AnnoStructs{Queryname: pair.queryname, Vs: finalVariants, Idx: pair.idx}

		// and we're done
		cVariants <- AS
	}
}
