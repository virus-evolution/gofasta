package cmd

import (
	"bufio"
	"errors"
	"os"
	"path/filepath"

	"github.com/spf13/cobra"

	"github.com/cov-ert/gofasta/pkg/gfio"
	"github.com/cov-ert/gofasta/pkg/updown"
)

var TRquery string
var TRtarget string
var TRignore string
var TRoutfile string

var TRsizetotal int
var TRsizeup int
var TRsizedown int
var TRsizeside int
var TRsizesame int

var TRnofill bool

var TRdistall int
var TRdistup int
var TRdistdown int
var TRdistside int

var TRdistpush int

var TRthresholdpair float32
var TRthresholdtarget int

func init() {
	updownCmd.AddCommand(toprankingCmd)

	toprankingCmd.Flags().StringVarP(&TRquery, "query", "q", "", "File with sequences to find neighbours for. Either the CSV output of gofasta updown list, or an alignment in fasta format")
	toprankingCmd.Flags().StringVarP(&TRtarget, "target", "t", "", "File of sequences to look for neighbours in. Either the CSV output of gofasta updown list, or an alignment in fasta format")
	toprankingCmd.Flags().StringVarP(&TRoutfile, "outfile", "o", "stdout", "CSV-format file of closest neighbours to write")
	toprankingCmd.Flags().StringVarP(&udReference, "reference", "r", "", "Reference sequence, in fasta format - only required if --query and --target are fasta files")
	toprankingCmd.Flags().StringVarP(&TRignore, "ignore", "", "", "Optional plain text file of IDs to ignore in the target file when searching for neighbours")

	toprankingCmd.Flags().IntVarP(&TRsizetotal, "size-total", "", 0, "Max number of neighbours to find (attempts to split equally between same/up/down/side). A hard limit")
	toprankingCmd.Flags().IntVarP(&TRsizeup, "size-up", "", 0, "Max number of closest parent sequences to find, if size-total not specified. A soft limit unless --no-fill")
	toprankingCmd.Flags().IntVarP(&TRsizedown, "size-down", "", 0, "Max number of closest child sequences to find, if size-total not specified. A soft limit unless --no-fill")
	toprankingCmd.Flags().IntVarP(&TRsizeside, "size-side", "", 0, "Max number of closest sibling sequences to find, if size-total not specified. A soft limit unless --no-fill")
	toprankingCmd.Flags().IntVarP(&TRsizesame, "size-same", "", 0, "Max number of identical sequences to find, if size-total not specified. A soft limit unless --no-fill")

	toprankingCmd.Flags().IntVarP(&TRdistall, "dist-all", "", 0, "Maximum allowed SNP-distance between target and query sequence in any direction. Overrides the settings below")
	toprankingCmd.Flags().IntVarP(&TRdistup, "dist-up", "", 0, "Maximum allowed SNP-distance from query for sequences in the parent bin")
	toprankingCmd.Flags().IntVarP(&TRdistdown, "dist-down", "", 0, "Maximum allowed SNP-distance from query for sequences in the child bin")
	toprankingCmd.Flags().IntVarP(&TRdistside, "dist-side", "", 0, "Maximum allowed SNP-distance from query for sequences in the sibling bin")

	toprankingCmd.Flags().Float32VarP(&TRthresholdpair, "threshold-pair", "", 0.1, "Up to this proportion of consequential sites is allowed to be ambiguous in either sequence for each pairwise comparison")
	toprankingCmd.Flags().IntVarP(&TRthresholdtarget, "threshold-target", "", 10000, "Target can have at most this number of ambiguities to be considered")

	toprankingCmd.Flags().BoolVarP(&TRnofill, "no-fill", "", false, "Don't make up for a shortfall in any of --size-up, -down, -side or -same by increasing the count for other bins")
	toprankingCmd.Flags().IntVarP(&TRdistpush, "dist-push", "", 0, "Push the SNP-distance boundaries for any empty bins to the closest n distances for which there are neighbours, where possible")

	toprankingCmd.Flags().Lookup("no-fill").NoOptDefVal = "true"
	toprankingCmd.Flags().Lookup("dist-push").NoOptDefVal = "1"

	toprankingCmd.Flags().SortFlags = false
}

var toprankingCmd = &cobra.Command{
	Use:   "topranking",
	Short: "get pseudo-tree-aware catchments for query sequences from alignments",
	Long: `get pseudo-tree-aware catchments for query sequences from alignments

Example usage:
	gofasta updown topranking -q smallquery.fasta -r WH04.fasta -t mutationlist.csv --size-total 1000 -o catchment.csv

For each sequence in --query, this routine finds the closest sequences by SNP-distance in --target, binned according to
whether they are likely children, parents, or siblings of, or on a polytomy with, the query sequence. It does this by comparing
SNPs relative to a common reference sequence which is imagined to be the root of the tree.

--query and --target can either be alignments in fasta format, or the CSV output of gofasta updown list, or one of each. They
must have file extensions .csv .fasta or .fa . If either is an alignment, you must provide --reference, and this should be the
same sequence that was used by gofasta updown list.

If you provide a number to --size-total, the output will try to include this many closest sequences to the query, split evenly
between the four bins. If one or more bins has a shortfall, the sizes of the other bins will increase until --size-total is met,
if possible, unless you use --no-fill.

You can provide --size-up, -down, -side and -same, instead of --size-total, if you want different proportions for each bin.
The program will aim to provide the sum of these numbers in total in the output, and will make up for a shortfall in one bin
by increasing the count of the other bins where possible, unless --no-fill.

You can also filter on SNP-distance instead of returning the closest n sequences. Use the --dist flags to do this. If --dist-push 
is invoked, the program will push the SNP-distance boundaries for any empty bins until it finds a neighbour, and then it will
return all the neighbours at that distance in that bin, and print what it's done to stderr.

You can combine the two types of flag (size and dist), to return only the closest n sequences under a set distance.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		var qtype string
		var ttype string

		switch filepath.Ext(TRquery) {
		case ".csv":
			qtype = "csv"
		case ".fasta":
			qtype = "fasta"
		case ".fa":
			qtype = "fasta"
		default:
			return errors.New("couldn't tell if --query was a .csv or a .fasta file")
		}

		switch filepath.Ext(TRtarget) {
		case ".csv":
			ttype = "csv"
		case ".fasta":
			ttype = "fasta"
		case ".fa":
			ttype = "fasta"
		default:
			return errors.New("couldn't tell if --target was a .csv or a .fasta file")
		}

		if (qtype == "fasta" || ttype == "fasta") && len(udReference) == 0 {
			return errors.New("if your either of your input files are fastas, you must provide a --reference")
		}

		// targets to ignore (potentially):
		ignoreArray := make([]string, 0)
		if len(TRignore) != 0 {
			f, err := os.Open(TRignore)
			if err != nil {
				return err
			}
			defer f.Close()
			s := bufio.NewScanner(f)
			for s.Scan() {
				ignoreArray = append(ignoreArray, s.Text())
			}
			err = s.Err()
			if err != nil {
				return err
			}
		}

		query, err := gfio.OpenIn(TRquery)
		if err != nil {
			return err
		}
		defer query.Close()

		target, err := gfio.OpenIn(TRtarget)
		if err != nil {
			return err
		}
		defer target.Close()

		ref, err := gfio.OpenIn(udReference)
		if err != nil {
			return err
		}
		defer ref.Close()

		ignore, err := gfio.OpenIn(TRignore)
		if err != nil {
			return err
		}
		defer ignore.Close()

		out, err := gfio.OpenOut(TRoutfile)
		if err != nil {
			return err
		}
		defer out.Close()

		err = updown.TopRanking(query, target, ref, out,
			qtype, ttype, ignoreArray,
			TRsizetotal, TRsizeup, TRsizedown, TRsizeside, TRsizesame,
			TRdistall, TRdistup, TRdistdown, TRdistside,
			TRthresholdpair, TRthresholdtarget, TRnofill, TRdistpush)

		return
	},
}
