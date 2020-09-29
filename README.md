# gofasta

Some of [Julialign](https://github.com/cov-ert/julialign), but in Go, for llama and civet.


Gofasta uses a slightly modified version of the bit-level coding scheme for nucleotides by Emmanuel Paradis (described [here](http://ape-package.ird.fr/misc/BitLevelCodingScheme.html), and implemented in the R package [ape](https://doi.org/10.1093/bioinformatics/btg412)).

### Installation

Binaries are available for Mac OS and Linux under the [latest release](https://github.com/cov-ert/gofasta/releases/latest)

Or you can get them from Conda: `conda install -c benjamincjackson gofasta`

Or if you have Go installed, you can run `go get github.com/cov-ert/gofasta` to build a binary in your `$GOPATH`


### Commands



| run              | description                                                                        |
|------------------|------------------------------------------------------------------------------------|
| `./gofasta closest`   | Find the closest sequence to a query by SNP-distance. Ties are broken by genome completeness (including for 0-length distances between genomes).                             |


