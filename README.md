# DupMethAnalyzer

Tool for analyzing DNA methylation differences between segmental duplication pairs in the genome.

## Description

DupMethAnalyzer examines DNA methylation patterns between segmental duplication pairs to identify epigenetic 
differences between highly identical duplicated regions. It uses minimap2 for sequence alignment and 
provides detailed methylation comparisons with visualizations.

## Installation

```bash
# Clone the repository
git clone https://github.com/hrrsjeong/DupMethAnalyzer.git
cd DupMethAnalyzer

# Install dependencies
pip install -r requirements.txt

# Ensure minimap2 is installed and in your PATH

# Usage
`python dupmethanalyzer.py --genome genome.fasta --methyl methylation.bed --dups segmental_duplications.bed --window 100 --threads 8
`
