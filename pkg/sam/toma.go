package sam

import (
	"context"
	"errors"
	"io"
	"os"

	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
)

// ToMultiAlign converts a SAM file containing pairwise alignments between assembled genomes to a fasta-format alignment.
// Insertions relative to the reference are discarded, so all the sequences are the same (=reference) length
// A query whose records can't be used is reported on stderr and skipped, unless strict, when it stops the run.
func ToMultiAlign(samIn io.Reader, out io.Writer, wrap int, trimstart int, trimend int, pad bool, threads int, strict bool) error {

	src, err := newSamSource(samIn)
	if err != nil {
		return err
	}
	header := src.header
	refLen := header.Refs()[0].Len()

	trimstart, trimend, trim, err := checkArgs(refLen, trimstart, trimend)
	if err != nil {
		return err
	}

	work := func(group samRecords) (fasta.Record, error) {
		group, err := prepareQuery(group, refLen)
		if err != nil {
			return fasta.Record{}, err
		}
		rawseq, err := getSeqFromBlock(group.records, refLen, false)
		if err != nil {
			return fasta.Record{}, pipeline.Named(group.records[0].Name, err)
		}
		return getRecord(rawseq, group.records[0].Name, group.idx, trim, pad, trimstart, trimend), nil
	}

	h := pipeline.NewHandler(strict, os.Stderr, queryNoun, func(record fasta.Record) error {
		return fasta.WriteRecord(out, record, wrap)
	})
	if err := pipeline.Run(context.Background(), threads, src.next, work, withNotes(src, h.Sink)); err != nil {
		return err
	}
	src.writeTail()
	return h.Finish()
}

// checkArgs sanity checks the trimming and padding arguments, given the length of the reference sequence
func checkArgs(refLen int, trimstart int, trimend int) (int, int, bool, error) {

	trim := false

	if trimstart == -1 {
		trimstart = 1
	} else {
		trim = true
	}

	if trimend == -1 {
		trimend = refLen
	} else {
		trim = true
	}

	if trimstart > refLen || trimstart < 1 {
		return 0, 0, false, errors.New("error parsing --start coordinate. --start must be > 0 && <= length(reference)")
	}
	if trimend > refLen || trimend < 1 {
		return 0, 0, false, errors.New("error parsing --end coordinate. --end must be > 0 &! > length(reference)")
	}
	if trimstart > trimend {
		return 0, 0, false, errors.New("error parsing trimming coordinates: --start must be <= --end")
	}

	return trimstart, trimend, trim, nil
}

// getRecord returns a Record struct with a sequence ID and a sequence
// that has been optionally trimmed and padded
func getRecord(rawseq []byte, id string, idx int, trim bool, pad bool, trimstart int,
	trimend int) fasta.Record {

	var seq []byte

	if pad {
		seq = swapInNs(rawseq)
	} else {
		seq = swapInGapsNs(rawseq)
	}

	if trim {
		if pad {
			for i, _ := range seq {
				if i < trimstart-1 || i >= trimend {
					seq[i] = 'N'
				}
			}
		} else {
			seq = seq[trimstart-1 : trimend]
		}
	}

	FR := fasta.Record{ID: id, Description: id, Seq: string(seq), Idx: idx}

	return FR
}
