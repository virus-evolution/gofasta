package cmd

import (
	"errors"
	"os"
	"path/filepath"

	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/sam"
)

var samVariantsAnnotation string
var samVariantsOutfile string
var samVariantsAggregate bool
var samVariantsThreshold float64
var samVariantsAppendSNP bool
var samVariantsStart int
var samVariantsEnd int

// for backwards compatibility:
var samVariantsGenbank string

func init() {
	samCmd.AddCommand(samVariantsCmd)

	samVariantsCmd.Flags().StringVarP(&samVariantsAnnotation, "annotation", "a", "", "Genbank or GFF3 format annotation file. Must have suffix .gb or .gff")
	samVariantsCmd.Flags().StringVarP(&samVariantsOutfile, "outfile", "o", "stdout", "Where to write the variants")
	samVariantsCmd.Flags().IntVarP(&samVariantsStart, "start", "", -1, "Only report variants after (and including) this position")
	samVariantsCmd.Flags().IntVarP(&samVariantsEnd, "end", "", -1, "Only report variants before (and including) this position")
	samVariantsCmd.Flags().BoolVarP(&samVariantsAggregate, "aggregate", "", false, "Report the proportions of each change")
	samVariantsCmd.Flags().Float64VarP(&samVariantsThreshold, "threshold", "", 0.0, "If --aggregate, only report changes with a freq greater than or equal to this value")
	samVariantsCmd.Flags().BoolVarP(&samVariantsAppendSNP, "append-snps", "", false, "Report the codon's SNPs in parenthesis after each amino acid mutation")

	samVariantsCmd.Flags().Lookup("aggregate").NoOptDefVal = "true"
	samVariantsCmd.Flags().Lookup("append-snps").NoOptDefVal = "true"

	samVariantsCmd.Flags().SortFlags = false

	samVariantsCmd.Flags().StringVarP(&samVariantsGenbank, "genbank", "g", "", "Genbank format annotation of a sequence in the same coordinates as the alignment")
	samVariantsCmd.Flags().MarkHidden("genbank")
}

var samVariantsCmd = &cobra.Command{
	Use:   "variants",
	Short: "Annotate mutations relative to a reference from an alignment in sam format",
	Long: `Annotate mutations relative to a reference from an alignment in sam format

Example usage:
	gofasta sam variants -s aligned.sam -r reference.fasta -a annotation.gb -o variants.csv
	gofasta sam variants -s aligned.sam -r reference.fasta -a annotation.gff -o variants.csv

If input sam and output csv files are not specified, the behaviour is to read the sam from stdin and write
the variants to stdout.

--reference should be the same sequence that was used to generate the sam file, and should be in the same coordinates
as the --annotation. You don't have to provide a file to --reference if your annotation has the fasta record in it.

gff-format annotations must be valid version 3 files. See github.com/virus-evolution/gofasta for more details
of the format.

Mutations are annotated with ins (insertion), del (deletion), aa (amino acid change) or nuc (a nucleotide change that
isn't in a codon that is represented by an amino acid change). The formats are:

	ins:2028:3 - a 3-base insertion immediately after (1-based) position 2028 in reference coordinates
	del:11288:9 - a 9-base deletion whose first missing nucleotide is at (1-based) position 11288 in reference coordinates
	aa:s:D614G - the amino acid at (1-based) residue 614 in the S gene is a D in the reference and a G in this sequence
	nuc:C3037T - the nucleotide at (1-based) position 3037 is a C in the reference and a T in this sequence

Frame-shifting mutations in coding sequence are reported as indels but are ignored for subsequent amino-acids in the alignment.
`,

	RunE: func(cmd *cobra.Command, args []string) (err error) {

		samIn, err := gfio.OpenIn(*cmd.Flag("samfile"))
		if err != nil {
			return err
		}
		defer samIn.Close()

		refFromFile := false
		var ref *os.File
		if samReference != "" {
			ref, err = gfio.OpenIn(*cmd.Flag("reference"))
			if err != nil {
				return err
			}
			refFromFile = true
		}
		defer ref.Close()

		// some backwards compatibility wrangling of --genbank vs --annotation
		var anno *os.File
		var annoSuffix string
		if samVariantsGenbank != "" {
			if samVariantsAnnotation != "" {
				return errors.New("--annotation replaces --genbank")
			}
			anno, err = gfio.OpenIn(*cmd.Flag("genbank"))
			if err != nil {
				return err
			}
			annoSuffix = "gb"
		} else {
			anno, err = gfio.OpenIn(*cmd.Flag("annotation"))
			if err != nil {
				return err
			}
			switch filepath.Ext(samVariantsAnnotation) {
			case ".gb":
				annoSuffix = "gb"
			case ".gff":
				annoSuffix = "gff"
			default:
				return errors.New("couldn't tell if --annotation was a .gb or a .gff file")
			}
		}
		defer anno.Close()

		out, err := gfio.OpenOut(*cmd.Flag("outfile"))
		if err != nil {
			return err
		}
		defer out.Close()

		err = sam.Variants(samIn, ref, refFromFile, anno, annoSuffix, out, samVariantsStart, samVariantsEnd, samVariantsAggregate, samVariantsThreshold, samVariantsAppendSNP, false, samThreads)

		return err
	},
}
