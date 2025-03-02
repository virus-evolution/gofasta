/*
 */
package gff

import (
	"bufio"
	"bytes"
	"errors"
	"io"
	"regexp"
	"strconv"
	"strings"

	"github.com/virus-evolution/gofasta/pkg/fasta"
)

type GFF struct {
	GFF_version     string
	HeaderLines     []string
	CommentLines    []string
	SequenceRegions map[string]SequenceRegion
	Features        []Feature
	IDmap           map[string][]int // a map from ID attribute tag to that ID's Feature line(s)
	FASTA           map[string]fasta.Record
}

type SequenceRegion struct {
	Seqid string
	Start int // 1-based start position of this region
	End   int // 1-based end position of this region
}

var (
	errGFFParsingVersion    = errors.New("Error parsing gff version")
	errGFFParsingSeqReg     = errors.New("Error parsing gff sequence-region")
	errGFFParsingSeqID      = errors.New("Error parsing gff SeqID")
	errGFFParsingStrand     = errors.New("Error parsing gff strand")
	errGFFParsingPhase      = errors.New("Error parsing gff phase")
	errGFFParsingAttributes = errors.New("Error parsing gff attributes")
)

func errorBuilder(err error, s string) error {
	return errors.New(err.Error() + ": " + s)
}

func (g *GFF) versionStringFromHeader() error {
	for _, line := range g.HeaderLines {
		if strings.HasPrefix(line, "gff-version") {
			fields := strings.Fields(line)
			if len(fields) != 2 {
				return errGFFParsingVersion
			}
			g.GFF_version = fields[1]
			return nil
		}
	}
	return errGFFParsingVersion
}

func (g *GFF) setSequenceRegionsFromHeader() error {
	for _, line := range g.HeaderLines {
		if strings.HasPrefix(line, "sequence-region") {
			fields := strings.Fields(line)
			if len(fields) < 4 {
				return errorBuilder(errGFFParsingSeqReg, line)
			}
			start, err := strconv.Atoi(fields[2])
			if err != nil {
				return err
			}
			end, err := strconv.Atoi(fields[3])
			if err != nil {
				return err
			}
			g.SequenceRegions[fields[1]] = SequenceRegion{Seqid: fields[1], Start: start, End: end}
		}
	}

	return nil
}

func (g *GFF) populateIDMap() {
	g.IDmap = make(map[string][]int)
	for i, f := range g.Features {
		if f.HasAttribute("ID") {
			if len(g.IDmap[f.Attributes["ID"][0]]) > 0 {
				g.IDmap[f.Attributes["ID"][0]] = append(g.IDmap[f.Attributes["ID"][0]], i)
			} else {
				g.IDmap[f.Attributes["ID"][0]] = []int{i}
			}
		}
	}
}

func (g *GFF) featuresOfType(t string) []Feature {

	fa := make([]Feature, 0)

	for _, f := range g.Features {
		if f.Type == t {
			fa = append(fa, f)
		}
	}

	return fa
}

func (g *GFF) featuresFromID(id string) []Feature {

	fa := make([]Feature, 0)

	if indices, ok := g.IDmap[id]; ok {
		for _, idx := range indices {
			fa = append(fa, g.Features[idx])
		}
	}

	return fa
}

type Feature struct {

	// The ID of the landmark used to establish the coordinate system for the current feature.
	// IDs may contain any characters, but must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|].
	// In particular, IDs may not contain unescaped whitespace and must not begin with an unescaped ">".
	Seqid string

	// The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this feature.
	// Typically this is the name of a piece of software, such as "Genescan" or a database name, such as "Genbank."
	// In effect, the source is used to extend the feature ontology by adding a qualifier to the type creating a new composite type
	// that is a subclass of the type in the type column.
	Source string

	// The type of the feature (previously called the "method"). This is constrained to be either a term from
	// the Sequence Ontology or an SO accession number. The latter alternative is distinguished using the syntax SO:000000.
	// In either case, it must be sequence_feature (SO:0000110) or an is_a child of it.
	Type string

	// The start and end coordinates of the feature are given in positive 1-based integer coordinates,
	// relative to the landmark given in column one. Start is always less than or equal to end.
	// For features that cross the origin of a circular feature (e.g. most bacterial genomes, plasmids,
	// and some viral genomes), the requirement for start to be less than or equal to end is satisfied
	// by making end = the position of the end + the length of the landmark feature.
	// For zero-length features, such as insertion sites, start equals end and the implied site is to
	// the right of the indicated base in the direction of the landmark.
	Start int
	End   int

	// The score of the feature, a floating point number. As in earlier versions of the format,
	// the semantics of the score are ill-defined. It is strongly recommended that E-values be used for
	// sequence similarity features, and that P-values be used for ab initio gene prediction features.
	// Score float64
	Score string // NB this should be float (above)

	// The strand of the feature. + for positive strand (relative to the landmark), - for minus strand,
	// and . for features that are not stranded. In addition, ? can be used for features whose strandedness is relevant,
	// but unknown.
	Strand string

	// For features of type "CDS", the phase indicates where the next codon begins relative to the 5' end
	// (where the 5' end of the CDS is relative to the strand of the CDS feature) of the current CDS feature.
	// For clarification the 5' end for CDS features on the plus strand is the feature's start and and the
	// 5' end for CDS features on the minus strand is the feature's end. The phase is one of the integers 0, 1, or 2,
	// indicating the number of bases forward from the start of the current CDS feature the next codon begins.
	// A phase of "0" indicates that a codon begins on the first nucleotide of the CDS feature (i.e. 0 bases forward),
	// a phase of "1" indicates that the codon begins at the second nucleotide of this CDS feature and a phase of "2"
	// indicates that the codon begins at the third nucleotide of this region. Note that ‘Phase’ in the context of a
	// GFF3 CDS feature should not be confused with the similar concept of frame that is also a common concept in
	// bioinformatics. Frame is generally calculated as a value for a given base relative to the start of the complete
	// open reading frame (ORF) or the codon (e.g. modulo 3) while CDS phase describes the start of the next codon
	// relative to a given CDS feature.
	// The phase is REQUIRED for all CDS features.
	Phase int

	// A list of feature attributes in the format tag=value. Multiple tag=value pairs are separated by semicolons.
	// URL escaping rules are used for tags or values containing the following characters: ",=;".
	// Spaces are allowed in this field, but tabs must be replaced with the %09 URL escape.
	// Attribute values do not need to be and should not be quoted. The quotes should be included as part of the
	// value by parsers and not stripped.
	// Multiple attributes of the same type are indicated by separating the values with the comma "," character. as in: Parent=AF2312,AB2812,abc-3
	// In addition to Parent, the Alias, Note, Dbxref and Ontology_term attributes can have multiple values.
	// All attributes that begin with an uppercase letter are reserved for later use. Attributes that begin with a lowercase letter can be used freely by applications.
	Attributes map[string][]string
}

/*
Attributes Tags with predefined meanings:
ID; Name; Alias; Parent; Target; Gap; Derives_from; Note; Dbxref; Ontology_term; Is_circular
see: https://github.com/The-Sequence-Ontology/Specifications/blob/fe73505276dd324bf6a55773f3413fe2bed47af4/gff3.md
*/

func (F *Feature) HasAttribute(tag string) bool {
	if _, ok := F.Attributes[tag]; ok {
		return true
	}
	return false
}

// ReadGFF reads a gff version 3 format annotation file and returns a struct that contains
// parsed versions of the fields it contains.
func ReadGFF(f io.Reader) (GFF, error) {

	gff := GFF{}
	gff.SequenceRegions = make(map[string]SequenceRegion)
	gff.IDmap = make(map[string][]int)

	header := make([]string, 0)
	comments := make([]string, 0)
	features := make([]Feature, 0)

	firstAfterHeader := true
	inFasta := false
	fastaBuffer := new(bytes.Buffer)

	var err error

	s := bufio.NewScanner(f)
	for s.Scan() {
		line := s.Text()
		if inFasta {
			fastaBuffer.Write([]byte(line + "\n"))
		} else if strings.HasPrefix(line, "##FASTA") {
			inFasta = true
		} else if strings.HasPrefix(line, "##") {
			header = append(header, strings.TrimPrefix(line, "##"))
		} else if strings.HasPrefix(line, "#") {
			comments = append(comments, strings.TrimSpace(strings.TrimPrefix(line, "#")))
		} else {
			if firstAfterHeader {
				gff.HeaderLines = header
				gff.CommentLines = comments
				err = gff.versionStringFromHeader()
				if err != nil {
					return gff, err
				}
				err = gff.setSequenceRegionsFromHeader()
				if err != nil {
					return gff, err
				}
				firstAfterHeader = false
			}
			feature, err := featureFromLine(line)
			if err != nil {
				return gff, err
			}
			features = append(features, feature)
		}
	}

	gff.Features = features
	gff.populateIDMap()

	if len(fastaBuffer.Bytes()) > 0 {
		fastamap := make(map[string]fasta.Record)
		fastaReader := bytes.NewReader(fastaBuffer.Bytes())
		EFRs, err := fasta.LoadEncodeAlignment(fastaReader, false, false, false)
		if err != nil {
			return gff, err
		}
		for _, EFR := range EFRs {
			fastamap[EFR.ID] = EFR.Decode()
		}

		gff.FASTA = fastamap
	}

	return gff, err
}

func featureFromLine(l string) (Feature, error) {

	F := Feature{}

	fields := strings.Split(l, "\t")

	var err error

	if len(fields) != 9 {
		return F, errors.New("GFF parsing error: wrong number of fields: " + l)
	}

	F.Seqid, err = seqidFromField(fields[0])
	if err != nil {
		return F, err
	}

	F.Source = fields[1]

	F.Type = fields[2]

	F.Start, err = strconv.Atoi(fields[3])
	if err != nil {
		return F, err
	}

	F.End, err = strconv.Atoi(fields[4])
	if err != nil {
		return F, err
	}

	// TO DO - this properly
	F.Score = fields[5]

	F.Strand, err = strandFromField(fields[6], l)
	if err != nil {
		return F, err
	}

	F.Phase, err = phaseFromField(F.Type, fields[7], l)
	if err != nil {
		return F, err
	}

	F.Attributes, err = attributesFromField(fields[8], l)
	if err != nil {
		return F, err
	}

	return F, nil
}

func isEscapedCorrectly(f string) error {

	if strings.HasPrefix(f, ">") {
		return errorBuilder(errGFFParsingSeqID, f)
	}

	matched, err := regexp.MatchString("[^\\][^a-zA-Z0-9.:^*$@!+_?-|]", f)
	if err != nil {
		return err
	}

	if matched {
		return errorBuilder(errGFFParsingSeqID, f)
	}

	return nil
}

func seqidFromField(f string) (string, error) {

	err := isEscapedCorrectly(f)
	if err != nil {
		return f, err
	}

	return f, nil
}

func strandFromField(f, l string) (string, error) {
	if f == "+" || f == "-" || f == "." || f == "?" {
		return f, nil
	}
	return f, errorBuilder(errGFFParsingStrand, l)
}

func phaseFromField(t, f, l string) (int, error) {
	if t == "CDS" {
		if result, err := strconv.Atoi(f); err == nil {
			if result >= 0 && result <= 2 {
				return result, err
			}
		}
		return 0, errorBuilder(errGFFParsingPhase, l)
	} else {
		if result, err := strconv.Atoi(f); err == nil {
			if result >= 0 && result <= 2 {
				return result, err
			}
		} else if f == "." {
			return 0, nil
		}
	}
	return 0, errorBuilder(errGFFParsingPhase, l)
}

func attributesFromField(f, l string) (map[string][]string, error) {
	m := make(map[string][]string)

	tagvaluepairs := strings.Split(f, ";")

	for _, tvp := range tagvaluepairs {
		tagvalues := strings.Split(tvp, "=")
		if len(tagvalues) != 2 {
			return m, errorBuilder(errGFFParsingAttributes, l)
		}
		m[tagvalues[0]] = strings.Split(tagvalues[1], ",")
	}

	return m, nil
}
