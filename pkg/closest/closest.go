/*
Package closest provides routines to find the closest sequences
to a set of query sequences, by genetic distance
*/
package closest

import (
	"errors"
	"fmt"
	"io"
	"math"
	"os"
	"runtime"
	"strconv"
	"strings"
	"sync"

	"github.com/virus-evolution/gofasta/pkg/encoding"
	"github.com/virus-evolution/gofasta/pkg/fastaio"
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

// scoreEncodedAlignment scores a single fasta record by how complete its genome is.
func scoreEncodedAlignment(cIn chan fastaio.EncodedFastaRecord, cOut chan fastaio.EncodedFastaRecord, measure string) {
	scoring := encoding.MakeEncodedScoreArray()
	var score int64

	for EFR := range cIn {
		score = 0
		for _, nuc := range EFR.Seq {
			score += scoring[nuc]
		}
		EFR.Score = score

		if measure == "tn93" {
			EFR.CalculateBaseContent()
		}

		cOut <- EFR
	}

	return
}

func rawDistance(query, target fastaio.EncodedFastaRecord) float64 {
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
func snpDistance(query, target fastaio.EncodedFastaRecord) float64 {
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
func tn93Distance(query, target fastaio.EncodedFastaRecord) float64 {

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
		if (query.Seq[i] & tNuc) < 16 { // are the bases different
			count_d++
			count_L++
			if (query.Seq[i] | tNuc) == 200 { // 1 if one of the bases is adenine and the other one is guanine, 0 otherwise
				count_P1++
			} else if (query.Seq[i] | tNuc) == 56 { // 1 if one of the bases is cytosine and the other one is thymine, 0 otherwise
				count_P2++
			}
		} else if query.Seq[i]&8 == 8 && query.Seq[i] == tNuc { // at the bases certainly the same
			count_L++
		}
	}

	// estimated rates from this pairwise comparison
	P1 := float64(count_P1) / float64(count_L)                   // rate of changes which are transitional differences between purines (A ⇄ G)
	P2 := float64(count_P2) / float64(count_L)                   // rate of changes which are transitional differences between pyramidines (C ⇄ T)
	Q := float64(count_d-(count_P1+count_P2)) / float64(count_L) // rate of changes which are transversional differences  (A ⇄ C || A ⇄ T || C ⇄ A || C ⇄ G) (i.e. everything else)

	// tidies up the equations a bit, after ape
	w1 := 1.0 - P1/k1 - Q/(2*g_R)
	w2 := 1.0 - P2/k2 - Q/(2*g_Y)
	w3 := 1.0 - Q/(2*g_R*g_Y)

	// tn93 distance:
	d := -k1*math.Log(w1) - k2*math.Log(w2) - k3*math.Log(w3)

	return d
}

// findClosest finds the single closest sequence by genetic distance among a set of target sequences to a query sequence
func findClosest(query fastaio.EncodedFastaRecord, measure string, cIn chan fastaio.EncodedFastaRecord, cOut chan resultsStruct) {
	var closest resultsStruct
	var distance float64
	var snps []string

	first := true

	decoding := encoding.MakeDecodingArray()

	for target := range cIn {

		switch measure {
		case "raw":
			distance = rawDistance(query, target)
		case "snp":
			distance = snpDistance(query, target)
		case "tn93":
			distance = tn93Distance(query, target)
		}

		if first {
			snps = make([]string, 0)
			for i, tNuc := range target.Seq {
				if (query.Seq[i] & tNuc) < 16 {
					snps = append(snps, strconv.Itoa(i+1)+decoding[query.Seq[i]]+decoding[tNuc])
				}
			}
			closest = resultsStruct{tname: target.ID, completeness: target.Score, distance: distance, snps: snps}
			first = false
			continue
		}

		if distance < closest.distance {
			snps = make([]string, 0)
			for i, tNuc := range target.Seq {
				if (query.Seq[i] & tNuc) < 16 {
					snps = append(snps, strconv.Itoa(i+1)+decoding[query.Seq[i]]+decoding[tNuc])
				}
			}
			closest = resultsStruct{tname: target.ID, completeness: target.Score, distance: distance, snps: snps}

		} else if distance == closest.distance {
			if target.Score > closest.completeness {
				snps = make([]string, 0)
				for i, tNuc := range target.Seq {
					if (query.Seq[i] & tNuc) < 16 {
						snps = append(snps, strconv.Itoa(i+1)+decoding[query.Seq[i]]+decoding[tNuc])
					}
				}
				closest = resultsStruct{tname: target.ID, completeness: target.Score, distance: distance, snps: snps}
			}
		}
	}

	closest.qname = query.ID
	closest.qidx = query.Idx

	cOut <- closest
}

// splitInput fans out target sequences over an array of query sequences, so that each target is passed over each query.
func splitInput(queries []fastaio.EncodedFastaRecord, measure string, cIn chan fastaio.EncodedFastaRecord, cOut chan resultsStruct, cErr chan error, cSplitDone chan bool) {

	nQ := len(queries)

	// make an array of channels, one for each query
	QChanArray := make([]chan fastaio.EncodedFastaRecord, nQ)
	for i := 0; i < nQ; i++ {
		QChanArray[i] = make(chan fastaio.EncodedFastaRecord)
	}

	for i, q := range queries {
		go findClosest(q, measure, QChanArray[i], cOut)
	}

	targetCounter := 0
	for EFR := range cIn {
		if targetCounter == 0 {
			if len(EFR.Seq) != len(queries[0].Seq) {
				cErr <- errors.New("query and target alignments are not the same width")
			}
		}
		targetCounter++

		for i, _ := range QChanArray {
			QChanArray[i] <- EFR
		}
	}

	fmt.Fprintf(os.Stderr, "number of sequences in target alignment: %d\n", targetCounter)

	for i := range QChanArray {
		close(QChanArray[i])
	}

	cSplitDone <- true
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

// Closest finds the single closest sequence by raw genetic distance to a query/queries. It writes the results
// to stdout or to file. Ties for distance are broken by genome completeness
func Closest(query, target io.Reader, measure string, out io.Writer, threads int) error {

	if threads == 0 {
		threads = runtime.NumCPU()
	} else if threads < runtime.NumCPU() {
		runtime.GOMAXPROCS(threads)
	}

	queries, err := fastaio.ReadEncodeAlignmentToList(query, false)
	if err != nil {
		return err
	}

	nQ := len(queries)

	fmt.Fprintf(os.Stderr, "number of sequences in query alignment: %d\n", nQ)

	QResultsArray := make([]resultsStruct, nQ)

	cErr := make(chan error)

	cTEFR := make(chan fastaio.EncodedFastaRecord, runtime.NumCPU())
	cTEFRscored := make(chan fastaio.EncodedFastaRecord, runtime.NumCPU())
	cTEFRdone := make(chan bool)
	cTEFRscoreddone := make(chan bool)
	cSplitDone := make(chan bool)

	cResults := make(chan resultsStruct)

	go fastaio.ReadEncodeAlignment(target, false, cTEFR, cErr, cTEFRdone)

	var wgScore sync.WaitGroup
	wgScore.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			scoreEncodedAlignment(cTEFR, cTEFRscored, measure)
			wgScore.Done()
		}()
	}

	go splitInput(queries, measure, cTEFRscored, cResults, cErr, cSplitDone)

	go func() {
		wgScore.Wait()
		cTEFRscoreddone <- true
	}()

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cTEFRdone:
			close(cTEFR)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cTEFRscoreddone:
			close(cTEFRscored)
			n--
		}
	}

	for n := 1; n > 0; {
		select {
		case err := <-cErr:
			return err
		case <-cSplitDone:
			n--
		}
	}

	for i := 0; i < nQ; i++ {
		result := <-cResults
		QResultsArray[result.qidx] = result
	}

	err = writeClosest(QResultsArray, measure, out)
	if err != nil {
		return err
	}

	return nil
}
