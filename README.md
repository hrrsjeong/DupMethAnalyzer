# DupMethAnalyzer

A tool for analyzing DNA methylation differences between segmental duplications in the genome.

Developed by Harris Jeong (2025)

## Overview

DupMethAnalyzer examines how DNA methylation patterns differ between highly identical duplicated regions in the genome. It provides detailed methylation comparisons with visualizations to help identify epigenetic differences that may be functionally significant.

## Installation

### Prerequisites

- Python 3.7 or higher
- minimap2

### Install from GitHub

```bash
# Clone the repository
git clone https://github.com/hrrsjeong/DupMethAnalyzer.git
cd DupMethAnalyzer

# Install dependencies
pip install -r requirements.txt
```
## Quick Start

```bash
python dupmethanalyzer.py --genome genome.fasta --methyl methylation.bed --dups segmental_duplications.bed 
```

## Input Files

DupMethAnalyzer requires three input files:

1. **Genome FASTA file** (`--genome`)
   - Contains the reference genome sequences
   - Standard FASTA format with chromosome sequences

2. **Methylation BED file** (`--methyl`)
   - Contains DNA methylation information for all CpG loci
   - Format: 5 columns (chromosome, start, end, methylated reads, total reads)
   - Example:
     ```
     chr1  10  11  8  10
     chr1  25  26  7  10
     ```

3. **Segmental Duplication BED file** (`--dups`)
   - Contains genomic coordinates of segmental duplicates
   - Format: 6 columns (chrom1, start1, end1, chrom2, start2, end2)
   - Example:
     ```
     chr1  1000  2000  chr2  3000  4000
     ```

## Detailed Usage Options

```
usage: dupmethanalyzer.py [-h] [--genome GENOME] [--methyl METHYL] [--dups DUPS] 
                          [--window WINDOW] [--output OUTPUT] [--minimap2 MINIMAP2]
                          [--threads THREADS] [--minimap2-threads MINIMAP2_THREADS]
                          [--keep-alignments] [--alignment-dir ALIGNMENT_DIR]
                          [--windows-bed WINDOWS_BED] [--version]
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--genome` | Path to genome FASTA file |
| `--methyl` | Path to methylation BED file |
| `--dups`   | Path to segmental duplications BED file |

### Optional Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--window` | Sliding window size in base pairs | 100 |
| `--output` | Output file path | methyl_diff_output.tsv |
| `--minimap2` | Path to minimap2 executable | minimap2 |
| `--threads` | Number of threads for parallel processing | 4 |
| `--minimap2-threads` | Number of threads per minimap2 process | 1 |
| `--keep-alignments` | Keep alignment files from minimap2 | False |
| `--alignment-dir` | Directory to store alignment files | alignments |
| `--windows-bed` | Output BED file with analyzed windows | None |
| `--version` | Show program version and exit | |

## Example Commands

### Basic Analysis
```bash
python dupmethanalyzer.py --genome genome.fasta --methyl methylation.bed --dups segmental_duplications.bed
```

## Output Files

### Main Output
The primary output is a TSV file containing:
- Window coordinates for each duplicated pair
- Sequence identity (overall and window-specific)
- Methylation levels for both regions
- Absolute methylation differences

### Additional Outputs
If requested, the tool also generates:
- Windows BED file with all analyzed regions
- Alignment files in PAF format
- Visualization plots:
  - Sequence identity vs. methylation difference
  - Methylation levels comparison
  - Sequence identity distribution

## Example Data

We provide example input files to test DupMethAnalyzer:

```bash
# Run with example data
cd examples
python ../dupmethanalyzer.py --genome toy_genome.fasta --methyl toy_methylation.bed --dups toy_segmental_dups.bed --window 100
```

## Contact

Harris Jeong  
hrrsjeong@gmail.com
