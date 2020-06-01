
# Python for Bioinformatics

Contains scripts coded in base python to solve basic biological data problems.


### Find Orthologs Script

Script for finding orthologs using reciprocal BLAST hits.

Command-line arguments:
- -i1: Input fasta file 1
- -i2: Input fasta file 2
- -o: Output file
- -t: Sequence type (n: nucleotide sequence; p: protein sequence)

**Input files:** Input fasta file 1, Input fasta file 2, Output file name

**Execution:** ./find_orthologs.py -i1 <Input file 1> -i2 <Input file 2> -o <Output file name> –t <Sequence type – n/p>

*blast+ needs to be installed to run the script*

**Output:** The output file contains a list of reciprocal BLAST hits, with each element as a tab separated pair of gene names.


### Kmer-Counter Script

Calculates the number of k-mers of a specific length in a given nucleotide sequence.

**Input:**
- FASTA file containing the nucleotide sequence
- k-mer length (k)

**Execution:** ./kmer_counter.py <k> <input FASTA file>

**Output:** Printed on standard output as tab separated columns, the first column represents the k-mer sequence and the second column has the number of times it occurs in the input file

## Sequence Alignment

### Global Alignment

nwAlign.py script performs global sequence alignment using the Needle-Wunsch Algorithm

**Input:**
- FASTA file containing sequence 1
- FASTA file containing sequence 2

**Execution:** ./nwAlign.py <seq1.fa> <seq2.fa>

**Output:** The aligned sequence and the score is printed to standard output

### Local alignment

swAlign.py script performs global sequence alignment using the Smith-Waterman Algorithm

**Input:**
- FASTA file containing sequence 1
- FASTA file containing sequence 2

**Execution:** ./swAlign.py <seq1.fa> <seq2.fa>

**Output:** The aligned sequence and the score is printed to standard output
