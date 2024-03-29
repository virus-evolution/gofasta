package cmd

import (
	"errors"

	"github.com/spf13/cobra"

	"github.com/virus-evolution/gofasta/pkg/gfio"
	"github.com/virus-evolution/gofasta/pkg/sam"
)

var toMultiAlignOutfile string
var toMultiAlignStart int
var toMultiAlignEnd int
var toMultiAlignPad bool
var toMultiAlignWrap int

// junk:
var toMultiAlignTrim bool
var toMultiAlignTrimStart int
var toMultiAlignTrimEnd int

func init() {
	samCmd.AddCommand(toMultiAlignCmd)

	toMultiAlignCmd.Flags().IntVarP(&toMultiAlignStart, "start", "", -1, "1-based first nucleotide position to retain in the output. Bases before this position are omitted, or are replaced with N if --pad")
	toMultiAlignCmd.Flags().IntVarP(&toMultiAlignEnd, "end", "", -1, "1-based last nucleotide position to retain the in output. Bases after this position are omitted, or are replaced with N if --pad")
	toMultiAlignCmd.Flags().BoolVarP(&toMultiAlignPad, "pad", "", false, "If --start and/or --end, replace the trimmed-out regions with Ns, else replace external deletions with Ns")
	toMultiAlignCmd.Flags().StringVarP(&toMultiAlignOutfile, "fasta-out", "o", "stdout", "Where to write the alignment")
	toMultiAlignCmd.Flags().IntVarP(&toMultiAlignWrap, "wrap", "w", -1, "Wrap the output alignment to this number of nucleotides wide. Omit this option not to wrap the output.")

	toMultiAlignCmd.Flags().BoolVarP(&toMultiAlignTrim, "trim", "", false, "Trim the alignment")
	toMultiAlignCmd.Flags().IntVarP(&toMultiAlignTrimStart, "trimstart", "", -1, "Start coordinate for trimming (0-based, half open)")
	toMultiAlignCmd.Flags().IntVarP(&toMultiAlignTrimEnd, "trimend", "", -1, "End coordinate for trimming (0-based, half open)")
	toMultiAlignCmd.Flags().MarkHidden("trim")
	toMultiAlignCmd.Flags().MarkHidden("trimstart")
	toMultiAlignCmd.Flags().MarkHidden("trimend")

	toMultiAlignCmd.Flags().SortFlags = false
}

var toMultiAlignCmd = &cobra.Command{
	Use:     "toMultiAlign",
	Aliases: []string{"tomultialign", "toma"},
	Short:   "Convert a SAM file to a multiple alignment in fasta format",
	Long: `Convert a SAM file to a multiple alignment in fasta format

Insertions relative to the reference are omitted, so all sequences in the output are the same ( = reference) length.

Example usage:
	gofasta sam toMultiAlign -s aligned.sam -o aligned.fasta

If you want, you can trim (and optionally pad) the output alignment to coordinates of your choosing:
	gofasta sam toMultiAlign -s aligned.sam --start 266 --end 29674 --pad -o aligned.fasta

If input and output files are not specified, the behaviour is to read the sam file from stdin and write
the fasta file to stdout, e.g.:
	minimap2 -a -x asm20 --score-N=0 reference.fasta unaligned.fasta | gofasta sam toMultiAlign > aligned.fasta`,

	RunE: func(cmd *cobra.Command, args []string) (err error) {

		/*
			Start of trimming argument reconciliation to maintain backwards compatibility
		*/
		if toMultiAlignTrim || toMultiAlignTrimStart != -1 || toMultiAlignTrimEnd != -1 {
			if toMultiAlignStart != -1 || toMultiAlignEnd != -1 {
				return errors.New(`As of version 1.0.0 --start and --end replace --trim, --trimstart and --trimend 
The old flags and their behaviour are retained for backwards compatibility, but you shouldn't combine the two sets of flags.
Note the change of coordinate system if moving from old to new flags.

`)
			}
		}

		if toMultiAlignTrimStart != -1 {
			toMultiAlignStart = toMultiAlignTrimStart + 1
		}

		if toMultiAlignTrimEnd != -1 {
			toMultiAlignEnd = toMultiAlignTrimEnd
		}
		/*
			End of trimming argument reconciliation to maintain backwards compatibility
		*/

		samIn, err := gfio.OpenIn(*cmd.Flag("samfile"))
		if err != nil {
			return err
		}
		defer samIn.Close()

		out, err := gfio.OpenOut(*cmd.Flag("fasta-out"))
		if err != nil {
			return err
		}
		defer out.Close()

		err = sam.ToMultiAlign(samIn, out, toMultiAlignWrap, toMultiAlignStart, toMultiAlignEnd, toMultiAlignPad, samThreads)

		return
	},
}
