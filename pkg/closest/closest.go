/*
Package closest provides routines to find the closest sequences
to a set of query sequences, by genetic distance
*/
package closest

import (
	"context"
	"errors"
	"fmt"
	"io"
	"math"
	"os"
	"runtime"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fasta"
	"github.com/virus-evolution/gofasta/pkg/pipeline"
)

// resultsStruct is a struct that contains information about a query sequence and its (current)
// single closest target by raw genetic distance, including the snps that distinguish them and
// the total snp-distance between them.
type resultsStruct struct {
	qname        string
	qidx         int
	tname        string
	completeness int64
	distance     float64
	snps         []string
}

// lessDistance reports whether distance a is strictly closer than distance b. An undefined
// distance (NaN, e.g. no comparable sites) is treated as further than every defined distance,
// and two undefined distances are never "less" than one another.
func lessDistance(a, b float64) bool {
	if math.IsNaN(a) {
		return false
	}
	if math.IsNaN(b) {
		return true
	}
	return a < b
}

func rawDistance(query, target fasta.EncodedRecord) float64 {
	n := 0
	d := 0
	for i, tNuc := range target.Seq {
		if (query.Seq[i] & tNuc) < 16 {
			n += 1
			d += 1
		}
		if (query.Seq[i]&8 == 8) && query.Seq[i] == tNuc {
			d += 1
		}
	}
	distance := float64(n) / float64(d)
	return distance
}

// TO DO - have this operate on the lists of snps not the entire sequences
func snpDistance(query, target fasta.EncodedRecord) float64 {
	n := 0
	for i, tNuc := range target.Seq {
		if (query.Seq[i] & tNuc) < 16 {
			n += 1
		}
	}
	return float64(n)
}

// See equation (7) in Tamura K, Nei M. Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees.
// Mol Biol Evol. 1993 May;10(3):512-26. doi: 10.1093/oxfordjournals.molbev.a040023. PMID: 8336541.
// See also ape: https://github.com/cran/ape/blob/c2fd899f66d6493a80484033772a3418e5d706a4/src/dist_dna.c
func tn93Distance(query, target fasta.EncodedRecord) float64 {

	// Total ATGC length of the two sequences
	L := float64(target.Count_A + target.Count_C + target.Count_G + target.Count_T + query.Count_A + query.Count_C + query.Count_G + query.Count_T)

	// estimates of the equilibrium base contents from the pair's sequence data
	g_A := float64(target.Count_A+query.Count_A) / L
	g_C := float64(target.Count_C+query.Count_C) / L
	g_G := float64(target.Count_G+query.Count_G) / L
	g_T := float64(target.Count_T+query.Count_T) / L

	g_R := float64(target.Count_A+query.Count_A+target.Count_G+query.Count_G) / L
	g_Y := float64(target.Count_C+query.Count_C+target.Count_T+query.Count_T) / L

	// tidies up the equations a bit, after ape
	k1 := 2.0 * g_A * g_G / g_R
	k2 := 2.0 * g_T * g_C / g_Y
	k3 := 2.0 * (g_R*g_Y - g_A*g_G*g_Y/g_R - g_T*g_C*g_R/g_Y)

	count_P1 := 0 // count of transitional differences between purines (A ⇄ G)
	count_P2 := 0 // count of transitional differences between pyramidines (C ⇄ T)

	count_d := 0 // total number of differences
	count_L := 0 // total length of resolved comparison

	// calculate the three types of change from the pairwise comparison
	for i, tNuc := range target.Seq {
		if (query.Seq[i]&tNuc) < 16 && query.Seq[i]&8 == 8 && tNuc&8 == 8 { // are the bases different (and known for sure)
			count_d++
			count_L++
			if (query.Seq[i] | tNuc) == 200 { // 1 if one of the bases is adenine and the other one is guanine, 0 otherwise
				count_P1++
			} else if (query.Seq[i] | tNuc) == 56 { // 1 if one of the bases is cytosine and the other one is thymine, 0 otherwise
				count_P2++
			}
		} else if query.Seq[i]&8 == 8 && query.Seq[i] == tNuc { // are the bases certainly the same
			count_L++
		}
	}

	// estimated rates from this pairwise comparison
	P1 := float64(count_P1) / float64(count_L)                   // rate of changes which are transitional differences between purines (A ⇄ G)
	P2 := float64(count_P2) / float64(count_L)                   // rate of changes which are transitional differences between pyramidines (C ⇄ T)
	Q := float64(count_d-(count_P1+count_P2)) / float64(count_L) // rate of changes which are transversional differences  (A ⇄ C || A ⇄ T || G ⇄ T || C ⇄ G) (i.e. everything else)

	// tidies up the equations a bit, after ape
	w1 := 1.0 - P1/k1 - Q/(2*g_R)
	w2 := 1.0 - P2/k2 - Q/(2*g_Y)
	w3 := 1.0 - Q/(2*g_R*g_Y)

	// tn93 distance:
	d := -k1*math.Log(w1) - k2*math.Log(w2) - k3*math.Log(w3)

	if d == 0.0 {
		d = 0.0
	}

	return d
}

// scoredTarget is one target with its distance to each query
type scoredTarget struct {
	target    fasta.EncodedRecord
	distances []float64
}

// loadQueries reads and encodes the queries. A query that can't be used (an invalid character, or a
// different width from the first usable query) is reported on stderr and left out, unless strict, when
// it is an error. It writes the number of queries read to stderr.
func loadQueries(query io.Reader, strict bool) ([]fasta.EncodedRecord, error) {
	width := -1 // of the first usable query; only one worker, so the queries are seen in order
	work := func(record fasta.Record) (fasta.EncodedRecord, error) {
		q, err := fasta.EncodeRecord(record, false, true, false)
		if err != nil {
			return q, pipeline.Named(record.ID, err)
		}
		if width == -1 {
			width = len(q.Seq)
		} else if len(q.Seq) != width {
			return q, pipeline.Named(record.ID, errors.New("not the same width as the first query"))
		}
		return q, nil
	}
	queries := make([]fasta.EncodedRecord, 0)
	h := pipeline.NewHandler(strict, os.Stderr, pipeline.Noun{One: "query", Many: "queries"}, func(q fasta.EncodedRecord) error {
		queries = append(queries, q)
		return nil
	})
	if err := pipeline.Run(context.Background(), 1, fasta.NewRecordSource(query), work, h.Sink); err != nil {
		return nil, err
	}
	if err := h.Finish(); err != nil {
		return nil, err
	}
	fmt.Fprintf(os.Stderr, "number of sequences in query alignment: %d\n", h.Total())
	return queries, nil
}

// scanTargets reads the targets from an alignment, works out each one's distance to every query (on
// threads goroutines), and passes the targets with their distances to visit in input order. A target
// that can't be used (an invalid character, or a different width from the queries) is reported on
// stderr and skipped, unless strict, when it is an error. It writes the number of targets read to stderr.
func scanTargets(queries []fasta.EncodedRecord, target io.Reader, measure string, threads int, strict bool, visit func(scoredTarget)) error {
	var distance func(query, target fasta.EncodedRecord) float64
	switch measure {
	case "raw":
		distance = rawDistance
	case "snp":
		distance = snpDistance
	case "tn93":
		distance = tn93Distance
	default:
		return errors.New("unknown distance measure: " + measure)
	}

	work := func(record fasta.Record) (scoredTarget, error) {
		t, err := fasta.EncodeRecord(record, false, true, true)
		if err != nil {
			return scoredTarget{}, pipeline.Named(record.ID, err)
		}
		if len(queries) > 0 && len(t.Seq) != len(queries[0].Seq) {
			return scoredTarget{}, pipeline.Named(record.ID, errors.New("query and target alignments are not the same width"))
		}
		ds := make([]float64, len(queries))
		for i, q := range queries {
			ds[i] = distance(q, t)
		}
		return scoredTarget{target: t, distances: ds}, nil
	}
	h := pipeline.NewHandler(strict, os.Stderr, pipeline.Noun{One: "target", Many: "targets"}, func(st scoredTarget) error {
		visit(st)
		return nil
	})
	if err := pipeline.Run(context.Background(), threads, fasta.NewRecordSource(target), work, h.Sink); err != nil {
		return err
	}
	if err := h.Finish(); err != nil {
		return err
	}
	fmt.Fprintf(os.Stderr, "number of sequences in target alignment: %d\n", h.Total())
	return nil
}

// snpsBetween lists the sites at which query and target certainly differ, as position, query base, target base
func snpsBetween(query, target fasta.EncodedRecord, decoding [256]string) []string {
	snps := make([]string, 0)
	for i, tNuc := range target.Seq {
		if (query.Seq[i] & tNuc) < 16 {
			snps = append(snps, strconv.Itoa(i+1)+decoding[query.Seq[i]]+decoding[tNuc])
		}
	}
	return snps
}

// consider makes target the closest to query if it is closer than the current closest, or as close and
// more complete. The first target considered always becomes the closest.
func (closest *resultsStruct) consider(query, target fasta.EncodedRecord, distance float64, first bool, decoding [256]string) {
	if first || lessDistance(distance, closest.distance) || (distance == closest.distance && target.Score > closest.completeness) {
		*closest = resultsStruct{qname: query.ID, qidx: query.Idx, tname: target.ID, completeness: target.Score, distance: distance, snps: snpsBetween(query, target, decoding)}
	}
}

// writeClosest parses an array of resultsStructs in order to write them, usually to stdout or file
func writeClosest(results []resultsStruct, measure string, w io.Writer) error {

	var err error

	_, err = w.Write([]byte("query,closest,distance,SNPs\n"))
	if err != nil {
		return err
	}

	for _, result := range results {
		switch measure {
		case "raw":
			w.Write([]byte(result.qname + "," + result.tname + "," + strconv.FormatFloat(result.distance, 'f', 9, 64) + "," + strings.Join(result.snps, ";") + "\n"))
		case "snp":
			w.Write([]byte(result.qname + "," + result.tname + "," + strconv.Itoa(int(result.distance)) + "," + strings.Join(result.snps, ";") + "\n"))
		case "tn93":
			w.Write([]byte(result.qname + "," + result.tname + "," + strconv.FormatFloat(result.distance, 'f', 9, 64) + "," + strings.Join(result.snps, ";") + "\n"))
		}
	}

	return nil
}

// Closest finds the single closest sequence by genetic distance to a query/queries. It writes the results
// to stdout or to file. Ties for distance are broken by genome completeness
func Closest(query, target io.Reader, measure string, out io.Writer, threads int, strict bool) error {

	if threads == 0 {
		threads = runtime.NumCPU()
	} else if threads < runtime.NumCPU() {
		runtime.GOMAXPROCS(threads)
	}

	queries, err := loadQueries(query, strict)
	if err != nil {
		return err
	}

	nQ := len(queries)

	QResultsArray := make([]resultsStruct, nQ)
	for i, q := range queries {
		QResultsArray[i] = resultsStruct{qname: q.ID, qidx: q.Idx}
	}

	decoding := encoding.MakeDecodingArray()
	first := true
	err = scanTargets(queries, target, measure, threads, strict, func(st scoredTarget) {
		for i, q := range queries {
			QResultsArray[i].consider(q, st.target, st.distances[i], first, decoding)
		}
		first = false
	})
	if err != nil {
		return err
	}

	err = writeClosest(QResultsArray, measure, out)
	if err != nil {
		return err
	}

	return nil
}
