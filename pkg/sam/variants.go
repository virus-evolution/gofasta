package sam

import (
	"context"
	"errors"
	"io"
	"os"
	"strconv"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/genbank"
	"github.com/virus-evolution/gofasta/pkg/gff"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
	"github.com/virus-evolution/gofasta/pkg/variants"
)

// Variants annotates amino acid, insertion, deletion, and nucleotide (anything
// outside of codons with an amino acid change) mutations relative to a reference
// sequence from pairwise alignments in sam format. Genome annotations are
// derived from a annotation file in genbank or gff version 3 format
func Variants(samIn, refIn io.Reader, refFromFile bool, annoIn io.Reader, annoSuffix string, out io.Writer, start, end int, aggregate bool, threshold float64, appendSNP bool, threads int) error {

	var ref fasta.EncodedRecord
	if refFromFile {
		refs, err := fasta.LoadEncodeAlignment(refIn, false, false, false)
		if err != nil {
			return err
		}
		if len(refs) > 1 {
			return errors.New("more than one record in --reference")
		}
		ref = refs[0]
	}

	var cdsregions []variants.Region
	var intregions []int
	switch annoSuffix {
	case "gb":
		gb, err := genbank.ReadGenBank(annoIn)
		if err != nil {
			return err
		}
		if !refFromFile {
			temp := fasta.Record{Seq: string(gb.ORIGIN), ID: "annotation_fasta"}
			ref, err = temp.Encode()
			if err != nil {
				return err
			}
			os.Stderr.WriteString("using --annotation fasta as reference\n")
		}
		refLenDegapped := len(ref.Decode().Degap().Seq)
		cdsregions, intregions, err = variants.RegionsFromGenbank(gb, refLenDegapped)
		if err != nil {
			return err
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
					ref = fasta.EncodedRecord{ID: "annotation_fasta", Seq: encodedrefseq}
				}
				os.Stderr.WriteString("using --annotation fasta as reference\n")
			default:
				return errors.New("more that one sequence in gff ##FASTA section")
			}
		}
		refSeqDegapped := ref.Decode().Degap().Seq
		cdsregions, intregions, err = variants.RegionsFromGFF(gff, refSeqDegapped)
		if err != nil {
			return err
		}
	}

	// the output header is written before the sam file is read
	var lw *variants.LineWriter
	var agg *variants.Aggregator
	var err error
	if aggregate {
		agg, err = variants.NewAggregator(out, start, end, appendSNP, ref.ID)
	} else {
		lw, err = variants.NewLineWriter(out, start, end, appendSNP, ref.ID)
	}
	if err != nil {
		return err
	}

	header, next, err := newSamSource(samIn)
	if err != nil {
		return err
	}
	refSeq := ref.Decode().Seq
	if header.Refs()[0].Len() != len(refSeq) {
		return errors.New("the length of the reference sequence (" + strconv.Itoa(len(refSeq)) + ") does not match the length in the sam header (" + strconv.Itoa(header.Refs()[0].Len()) + ")")
	}

	EA := encoding.MakeEncodingArray()
	work := func(group samRecords) (variants.AnnoStructs, error) {
		pair, err := blockToPairwiseAlignment(group, []byte(refSeq), false)
		if err != nil {
			return variants.AnnoStructs{}, err
		}
		for i, nuc := range pair.query {
			pair.query[i] = EA[nuc]
		}
		for i, nuc := range pair.ref {
			pair.ref[i] = EA[nuc]
		}
		offsetRefCoord, _ := variants.GetMSAOffsets(pair.ref)
		return variants.GetVariantsPair(pair.ref, pair.query, pair.refname, pair.queryname, pair.idx, cdsregions, intregions, offsetRefCoord)
	}

	if aggregate {
		err = pipeline.Run(context.Background(), threads, next, work, pipeline.FirstError(agg.Add))
		if err != nil {
			return err
		}
		return agg.Write(threshold)
	}
	return pipeline.Run(context.Background(), threads, next, work, pipeline.FirstError(lw.Write))
}
