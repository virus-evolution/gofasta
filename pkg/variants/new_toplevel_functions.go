package variants

import (
	"bytes"
	"errors"
	"io"
	"os"
	"sort"
	"sync"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fastaio"
	"github.com/virus-evolution/gofasta/pkg/genbank"
	"github.com/virus-evolution/gofasta/pkg/gff"
)

func Variants2(msaIn io.Reader, stdin bool, refID string, annoIn io.Reader, annoSuffix string, out io.Writer, aggregate bool, threshold float64, appendSNP bool, threads int) error {

	var err error

	// Find the reference
	var ref fastaio.EncodedFastaRecord

	// Have to move the reader back to the beginning of the alignment, because we are scanning through it twice
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

	cMSA := make(chan fastaio.EncodedFastaRecord, 50+threads)
	cErr := make(chan error)
	cMSADone := make(chan bool)

	go fastaio.ReadEncodeAlignment(msaIn, false, cMSA, cErr, cMSADone)

	firstmissing := false

	if stdin && refID != "" {
		select {
		case ref = <-cMSA:
			if ref.ID != refID {
				return errors.New("--reference is not the first record in --msa (if --msa is reference-length you don't need to provide --reference)")
			}
			firstmissing = true
		case err := <-cErr:
			return err
		case <-cMSADone:
			return errors.New("is the pipe to --msa empty?") // TO DO - does this work/is this necessary?
		}
	}

	var cdsregions []Region2
	var intregions []int
	var refToMSA, MSAToRef []int

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
			ref = fastaio.EncodedFastaRecord{ID: "annotation_fasta", Seq: encodedrefseq}
		}

		refLenDegapped := len(ref.Decode().Degap().Seq)

		// get a list of CDS + intergenic regions from the genbank file
		cdsregions, intregions, err = RegionsFromGenbank(gb, refLenDegapped)
		if err != nil {
			return err
		}

		os.Stderr.WriteString("using --annotation fasta as reference\n")

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
					ref = fastaio.EncodedFastaRecord{ID: "annotation_fasta", Seq: encodedrefseq}
				}
				os.Stderr.WriteString("using --annotation fasta as reference\n")
			default:
				return errors.New("more that one sequence in gff ##FASTA section")
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
		go AggregateWriteVariants(out, appendSNP, threshold, ref.ID, cVariants, cWriteDone, cErr)
	case false:
		go WriteVariants(out, firstmissing, appendSNP, ref.ID, cVariants, cWriteDone, cErr)
	}

	var wgVariants sync.WaitGroup
	wgVariants.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			getVariants2(ref, cdsregions, intregions, refToMSA, MSAToRef, cMSA, cVariants, cErr)
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

func getVariants2(ref fastaio.EncodedFastaRecord, cdsregions []Region2, intregions []int, offsetRefCoord []int, offsetMSACoord []int, cMSA chan fastaio.EncodedFastaRecord, cVariants chan AnnoStructs, cErr chan error) {

	for record := range cMSA {

		AS, err := GetVariantsPair2(ref.Seq, record.Seq, ref.ID, record.ID, record.Idx, cdsregions, intregions, offsetRefCoord, offsetMSACoord)
		if err != nil {
			cErr <- err
			break
		}

		cVariants <- AS
	}
}

func GetVariantsPair2(ref, query []byte, refID, queryID string, idx int, cdsregions []Region2, intregions []int, offsetRefCoord []int, offsetMSACoord []int) (AnnoStructs, error) {

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
