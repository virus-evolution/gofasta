# gofasta

A command-line utility for genomic epidemiology, developed to handle SARS-CoV-2 alignments as part of the [COG-UK project](https://www.cogconsortium.uk/).

### Third party licences / acknowledgements

Gofasta uses a slightly modified version of the bit-level coding scheme for nucleotides by Emmanuel Paradis (described [here](http://ape-package.ird.fr/misc/BitLevelCodingScheme.html), and implemented in the R package [ape](https://doi.org/10.1093/bioinformatics/btg412)).

Gofasta also incorporates [bíogo](https://github.com/biogo/biogo), which is distributed under licence. Its licence is reproduced under `THIRD_PARTY_LICENCES/biogo` or run `gofasta licences` to print it.

The functions with SAM files as input were written exclusively to handle pairwise alignments between assembled SARS-CoV-2 genomes from [minimap2](https://github.com/lh3/minimap2).

### Installation

Binaries are available for Mac OS and Linux under the [latest release](https://github.com/virus-evolution/gofasta/releases/latest).

Or you can get them from Conda:

`conda install bioconda::gofasta`

Or if you have Go installed, you can run `go install github.com/virus-evolution/gofasta@latest` to build a binary of the latest release locally (maybe in `~/go/bin/`).

You can also build the current contents of this repository:

```
git clone git@github.com:virus-evolution/gofasta.git
cd gofasta
go build -o gofasta
```


### Commands

For a full list of commands and options, run `gofasta` with the `-h` flag, for example: `gofasta -h`,  `gofasta sam -h`, `gofasta sam variants -h`, etc.


| subcommand         | description                                                                                                                                                     |
|--------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| licences           | Print gofasta's and   third-party licence information                                                                                                           |
| closest            | Find the closest   sequence(s) to a query by raw genetic distance. Ties are broken by genome   completeness (including for 0-length distances between genomes). |
| sam   toMultiAlign | Convert a SAM file   to a multiple alignment in fasta format. Insertions relative to the reference   are discarded.                                             |
| sam   toPairAlign  | Convert a SAM file   to pairwise alignments in fasta format, by default including insertions   relative to the reference sequence.                              |
| sam   variants     | As `variants` but with a SAM file as input                                                                                                                      |
| snps               | List all nucleotide   changes relative to a reference sequence.                                                                                                 |
| updown             | Tools for   pseudo-tree-aware SNP distances between sequences                                                                                                   |
| variants           | For each sequence in a multiple sequence alignment, list the   amino acid, nucleotide and indel changes relative to an annotated reference   sequence.          |

