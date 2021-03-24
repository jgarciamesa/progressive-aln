# Nucleotide Progressive Aligner

Progressive alignment of nucleotide data using Needleman-Wunsch pairwise algorithm with sequences
encoded as profile matrices written in R.

## Dependencies (libraries)

* `Biostrings` : use Needleman-Wunsch to construct initial distance matrix
* `seqinr`: read, write, and general manipulation of FASTA sequences
* `ape`: create guide tree

## Usage

```Rscript --vanilla nuc_progressive_aln.R file.fasta```
