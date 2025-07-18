# codon_optimization
Codon optimization

## Installation
To install:

```
# Clone the repository
git clone https://github.com/szsctt/codon_optimization.git
cd codon_optimization

# Install the package
pip install .
```

To install dev version (for testing)

```
# Clone the repository
git clone https://github.com/szsctt/codon_optimization.git
cd codon_optimization

# Install in development mode with test dependencies
pip install -e ".[test]"
```

Verify installation

```
codon_optimise --help
```

## Usage

`codon_optimise` provides two main commands:

### Optimize Codons

```bash
codon_optimise opt [OPTIONS] INPUT_FILE
```

Optimizes codons in the input sequence file.

Arguments:

INPUT_FILE: Fasta file containing sequence to optimize
Options:

 - `--output-file TEXT`: Output file (default: codon_optimised.fasta)
 - `--type TEXT`: Type of sequence: 'dna' or 'protein' (default: dna)
 - `--before TEXT`: Sequence to add before the optimized sequence
 - `--after TEXT`: Sequence to add after the optimized sequence
 - `--avoid TEXT`: List of sequences to avoid in the optimized sequence. Use 'N' for any nucleotide, e.g. 'GGCC' or 'GG.NN'. Can be specified multiple times.

Example:

```
# Basic codon optimization of a DNA sequence
codon_optimise opt input.fasta --output-file optimized.fasta

# Optimize a protein sequence and avoid specific motifs
codon_optimise opt protein.fasta --type protein --avoid GGCC --avoid AATAAA
```
### Check optmised sequences

```
codon_optimise check [OPTIONS] REFERENCE_FILE OPTIMISED_FILES...
```

Checks optimized sequences against a reference sequence.

Arguments:

`REFERENCE_FILE`: Reference fasta file containing original sequences
`OPTIMISED_FILES`: One or more optimized fasta files to check

```
# Check a single optimized file against reference
codon_optimise check reference.fasta optimized.fasta

# Check multiple optimized files with detailed output
codon_optimise check reference.fasta opt1.fasta opt2.fasta opt3.fasta --detailed
```

### Input Format
Input files should be in FASTA format:

```
>Sequence_name
ATGCAGTCGTCGTAGCTAGCTAGCTAGCTAG
```

For protein sequences, use standard single-letter amino acid codes:

```
>Protein_name
MSRRSSSKLVPRGS
```