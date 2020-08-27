# gofasta

Some of [Julialign](https://github.com/cov-ert/julialign), but in Go, for llama and civet.


Gofasta uses a slightly modified version of the bit-level coding scheme for nucleotides by Emmanuel Paradis (described [here](http://ape-package.ird.fr/misc/BitLevelCodingScheme.html), and implemented in the R package [ape](https://doi.org/10.1093/bioinformatics/btg412)).

### Commands



| run              | description                                                                        |
|------------------|------------------------------------------------------------------------------------|
| `./gofasta closest`   | Find the closest sequence to a query by SNP-distance. Ties are broken by genome completeness (including for 0-length distances between genomes).                             |


