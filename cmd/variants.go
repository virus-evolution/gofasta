package cmd

import (
	"errors"
	"os"
	"path/filepath"

	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/variants"
)

var variantsMSA string
var variantsReference string
var variantsAnnotation string
var variantsOutfile string
var variantsThreads int
var variantsAggregate bool
var variantsThreshold float64
var variantsAppendSNP bool
var variantsStart int
var variantsEnd int

// for backwards compatibility:
var variantsGenbank string

func init() {
	rootCmd.AddCommand(variantsCmd)

	variantsCmd.Flags().StringVarP(&variantsMSA, "msa", "", "stdin", "Multiple sequence alignment in fasta format")
	variantsCmd.Flags().StringVarP(&variantsReference, "reference", "r", "", "The ID of the reference record in the msa")
	variantsCmd.Flags().StringVarP(&variantsAnnotation, "annotation", "a", "", "Genbank or GFF3 format annotation file. Must have suffix .gb or .gff")
	variantsCmd.Flags().StringVarP(&variantsOutfile, "outfile", "o", "stdout", "Name of the file of variants to write")
	variantsCmd.Flags().IntVarP(&variantsStart, "start", "", -1, "Only report variants after (and including) this position")
	variantsCmd.Flags().IntVarP(&variantsEnd, "end", "", -1, "Only report variants before (and including) this position")
	variantsCmd.Flags().BoolVarP(&variantsAggregate, "aggregate", "", false, "Report the proportions of each change")
	variantsCmd.Flags().Float64VarP(&variantsThreshold, "threshold", "", 0.0, "If --aggregate, only report changes with a freq greater than or equal to this value")
	variantsCmd.Flags().BoolVarP(&variantsAppendSNP, "append-snps", "", false, "Report the codon's SNPs in parenthesis after each amino acid mutation")
	variantsCmd.Flags().IntVarP(&variantsThreads, "threads", "t", 1, "Number of threads to use")

	variantsCmd.Flags().Lookup("aggregate").NoOptDefVal = "true"
	variantsCmd.Flags().Lookup("append-snps").NoOptDefVal = "true"

	variantsCmd.Flags().StringVarP(&variantsGenbank, "genbank", "", "", "Genbank format annotation")
	variantsCmd.Flags().MarkHidden("genbank")

	variantsCmd.Flags().SortFlags = false
}

var variantsCmd = &cobra.Command{
	Use:   "variants",
	Short: "Annotate mutations relative to a reference from a multiple sequence alignment in fasta format",
	Long: `Annotate mutations relative to a reference from a multiple sequence alignment in fasta format

Example usage:

	./gofasta variants --msa alignment.fasta --annotation MN908947.gb --reference MN908947.3 > variants.csv
	./gofasta variants --msa another.fasta --annotation MN908947.gff --reference Wuhan-Hu-1 > variants.csv

--reference is the name of the reference record in --msa. If you are reading the --msa from stdin, it must
be the first sequence. If you don't provide a --reference the program will try to use the fasta record in the
annotation file, in which case the --msa must be in the same coordinates.

gff-format annotations must be valid version 3 files. See github.com/virus-evolution/gofasta for more details
of the format.

If input --msa and output csv files are not specified, the behaviour is to read the alignment from stdin and write
the variants to stdout.

You can use --aggregate to report the overall proportions of each mutation in the --msa, and --threshold to filter on 
frequency.

Mutations are annotated with ins (insertion), del (deletion), aa (amino acid change) or nuc (a nucleotide change that
isn't in a codon that is represented by an amino acid change). The formats are:

	ins:2028:3 - a 3-base insertion immediately after (1-based) position 2028 in reference coordinates
	del:11288:9 - a 9-base deletion whose first missing nucleotide is at (1-based) position 11288 in reference coordinates
	aa:s:D614G - the amino acid at (1-based) residue 614 in the S gene is a D in the reference and a G in this sequence
	nuc:C3037T - the nucleotide at (1-based) position 3037 in reference coordinates is a C in the reference and a T in this sequence

Frame-shifting mutations in coding sequence are reported as indels but are ignored for subsequent amino-acids in the alignment.	
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {

		msa, err := gfio.OpenIn(*cmd.Flag("msa"))
		if err != nil {
			return err
		}
		defer msa.Close()

		stdin := false
		if variantsMSA == "stdin" {
			stdin = true
		}

		// some backwards compatibility wrangling of --genbank vs --annotation
		var anno *os.File
		var annoSuffix string
		if variantsGenbank != "" {
			if variantsAnnotation != "" {
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
			switch filepath.Ext(variantsAnnotation) {
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

		err = variants.Variants(msa, stdin, variantsReference, anno, annoSuffix, out, variantsStart, variantsEnd, variantsAggregate, variantsThreshold, variantsAppendSNP, variantsThreads)

		return
	},
}
