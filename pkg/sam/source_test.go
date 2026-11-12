package sam

import (
	"errors"
	"io"
	"reflect"
	"strings"
	"testing"

	"github.com/virus-evolution/gofasta/pkg/pipeline"
)

// samLineFor returns a sam alignment line for query name with the given flag and position
func samLineFor(name, flag, pos string) string {
	return name + "\t" + flag + "\tr\t" + pos + "\t60\t10M\t*\t0\t0\tATGAAACCCG\t*"
}

// readAllQueries calls src.next until io.EOF and returns the queries in order
func readAllQueries(t *testing.T, src *samSource) []samRecords {
	t.Helper()
	var groups []samRecords
	for {
		g, err := src.next()
		if err == io.EOF {
			return groups
		}
		if err != nil {
			t.Fatalf("next: unexpected error: %v", err)
		}
		groups = append(groups, g)
	}
}

func TestSamSourceGroupsConsecutiveRecordsByName(t *testing.T) {
	in := "@SQ\tSN:r\tLN:45\n" + strings.Join([]string{
		samLineFor("q1", "0", "1"),
		samLineFor("q1", "2048", "20"),
		samLineFor("q2", "0", "5"),
		samLineFor("q3", "0", "1"),
		samLineFor("q3", "2048", "30"),
		samLineFor("q3", "2048", "10"),
	}, "\n") + "\n"
	src, err := newSamSource(strings.NewReader(in))
	if err != nil {
		t.Fatal(err)
	}
	groups := readAllQueries(t, src)

	type summary struct {
		idx  int
		name string
		pos  []int
	}
	var got []summary
	for _, g := range groups {
		if g.err != nil {
			t.Errorf("%s: unexpected error: %v", g.name, g.err)
		}
		s := summary{idx: g.idx, name: g.name}
		for _, r := range g.records {
			s.pos = append(s.pos, r.Pos)
		}
		got = append(got, s)
	}
	want := []summary{
		{0, "q1", []int{0, 19}},
		{1, "q2", []int{4}},
		{2, "q3", []int{0, 29, 9}},
	}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, want %v", got, want)
	}
}

// A line that can't be parsed becomes the error of the query it belongs to (by QNAME); that query's other
// lines are still consumed, so the queries after it are read normally.
func TestSamSourceParseErrorStaysWithItsQuery(t *testing.T) {
	in := "@SQ\tSN:r\tLN:45\n" + strings.Join([]string{
		samLineFor("q1", "0", "1"),
		samLineFor("q2", "0", "1"),
		samLineFor("q2", "notaflag", "1"),
		samLineFor("q2", "2048", "20"),
		samLineFor("q3", "0", "1"),
	}, "\n") + "\n"
	src, err := newSamSource(strings.NewReader(in))
	if err != nil {
		t.Fatal(err)
	}
	groups := readAllQueries(t, src)
	if len(groups) != 3 {
		t.Fatalf("got %d queries, want 3", len(groups))
	}
	if groups[0].name != "q1" || groups[0].err != nil || len(groups[0].records) != 1 {
		t.Errorf("q1: got %+v", groups[0])
	}
	if groups[1].name != "q2" || groups[1].err == nil {
		t.Errorf("q2: want a parse error, got name %q err %v", groups[1].name, groups[1].err)
	}
	if len(groups[1].records) != 2 {
		t.Errorf("q2: got %d parsed records, want 2", len(groups[1].records))
	}
	if groups[2].name != "q3" || groups[2].err != nil || groups[2].idx != 2 {
		t.Errorf("q3: got %+v", groups[2])
	}

	// the error makes prepareQuery fail with the query's name attached
	_, err = prepareQuery(groups[1], 45)
	var recErr *pipeline.RecordError
	if !errors.As(err, &recErr) || recErr.Name != "q2" {
		t.Errorf("prepareQuery(q2): want a RecordError named q2, got %#v", err)
	}
}

// Blank lines are ignored, CRLF line endings are accepted, and a last line without a newline is read.
func TestSamSourceBlankLinesCRLFAndNoFinalNewline(t *testing.T) {
	in := "@SQ\tSN:r\tLN:45\r\n" +
		"\n" + samLineFor("q1", "0", "1") + "\r\n" +
		"\r\n\n" + samLineFor("q2", "0", "1") + "\n" +
		"\n" + samLineFor("q3", "0", "1")
	src, err := newSamSource(strings.NewReader(in))
	if err != nil {
		t.Fatal(err)
	}
	groups := readAllQueries(t, src)
	var names []string
	for _, g := range groups {
		if g.err != nil {
			t.Errorf("%s: unexpected error: %v", g.name, g.err)
		}
		names = append(names, g.name)
	}
	if want := []string{"q1", "q2", "q3"}; !reflect.DeepEqual(names, want) {
		t.Errorf("got queries %v, want %v", names, want)
	}
}

func TestSamSourceHeaderOnly(t *testing.T) {
	src, err := newSamSource(strings.NewReader("@SQ\tSN:r\tLN:45\n"))
	if err != nil {
		t.Fatal(err)
	}
	if groups := readAllQueries(t, src); len(groups) != 0 {
		t.Errorf("got %d queries, want 0", len(groups))
	}
}

// The messages about skipped lines are kept with the first query that comes after them in the file, and
// those after the last query go to the tail.
func TestSamSourceNotesOrdering(t *testing.T) {
	in := "@SQ\tSN:r\tLN:45\n@SQ\tSN:s\tLN:90\n" + strings.Join([]string{
		"u0\t4\t*\t0\t0\t*\t*\t0\t0\tATGAAACCCG\t*", // before q1: noted with q1
		samLineFor("q1", "0", "1"),
		samLineFor("q1", "256", "5"), // after q1 has started: noted with the next query
		samLineFor("q1", "2048", "20"),
		"o1\t0\ts\t1\t60\t10M\t*\t0\t0\tATGAAACCCG\t*", // other @SQ: noted with q2
		samLineFor("q2", "0", "1"),
		"u2\t4\t*\t0\t0\t*\t*\t0\t0\tATGAAACCCG\t*", // after the last query: tail
	}, "\n") + "\n"
	src, err := newSamSource(strings.NewReader(in))
	if err != nil {
		t.Fatal(err)
	}
	groups := readAllQueries(t, src)
	if len(groups) != 2 || len(groups[0].records) != 2 || len(groups[1].records) != 1 {
		t.Fatalf("unexpected queries: %+v", groups)
	}

	wantNotes := map[int][]string{
		0: {"skipping unmapped read: u0\n"},
		1: {
			"ignoring secondary mapping: q1\n",
			"skipping record mapped to a reference other than the first @SQ (s): o1\n",
		},
	}
	if !reflect.DeepEqual(src.notes, wantNotes) {
		t.Errorf("notes: got %q, want %q", src.notes, wantNotes)
	}
	if want := []string{"skipping unmapped read: u2\n"}; !reflect.DeepEqual(src.tail, want) {
		t.Errorf("tail: got %q, want %q", src.tail, want)
	}
}
