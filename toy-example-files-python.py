#!/usr/bin/env python3
"""
Generate more realistic toy example files for segmental duplication methylation analysis.
This script creates:
1. A genome with actual segmental duplications (with mutations)
2. Methylation data with realistic patterns
3. A segmental duplication BED file
"""

import random
import os

def generate_random_sequence(length):
    """Generate a random DNA sequence of the given length."""
    return ''.join(random.choice('ACGT') for _ in range(length))

def introduce_mutations(sequence, mutation_rate=0.05):
    """Introduce random mutations in the sequence at the given rate."""
    mutated = list(sequence)
    for i in range(len(mutated)):
        if random.random() < mutation_rate:
            # Replace with a different nucleotide
            current = mutated[i]
            options = [b for b in 'ACGT' if b != current]
            mutated[i] = random.choice(options)
    return ''.join(mutated)

def introduce_indels(sequence, indel_rate=0.02, max_indel_size=5):
    """Introduce insertions and deletions in the sequence."""
    mutated = list(sequence)
    i = 0
    while i < len(mutated):
        if random.random() < indel_rate/2:  # Insertion
            size = random.randint(1, max_indel_size)
            insertion = generate_random_sequence(size)
            mutated[i:i] = insertion
            i += size
        elif random.random() < indel_rate/2:  # Deletion
            size = random.randint(1, min(max_indel_size, len(mutated) - i))
            mutated[i:i+size] = []
        else:
            i += 1
    return ''.join(mutated)

def generate_methylation_pattern(length, cpg_density=0.05, mean_methylation=0.7, std_dev=0.2):
    """Generate methylation data with a realistic correlation pattern."""
    # First, determine CpG positions based on density
    cpg_positions = []
    for i in range(length):
        if random.random() < cpg_density:
            cpg_positions.append(i)
    
    # Generate correlated methylation values
    methylation_values = []
    current_value = random.gauss(mean_methylation, std_dev)
    
    for _ in cpg_positions:
        # Ensure value stays between 0 and 1
        current_value = max(0, min(1, current_value + random.gauss(0, 0.1)))
        methylation_values.append(current_value)
    
    return list(zip(cpg_positions, methylation_values))

def create_toy_files(output_dir):
    """Create toy example files for testing."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Parameters
    chromosome_length = 10000
    num_reads = 20  # Total reads for methylation data
    
    # Create segmental duplications
    duplications = [
        # (chrom1, start1, end1, chrom2, start2, end2, identity)
        ("chr1", 1000, 3000, "chr2", 2000, 4000, 0.95),
        ("chr1", 5000, 7000, "chr3", 1000, 3000, 0.90),
        ("chr2", 6000, 8000, "chr3", 5000, 7000, 0.85)
    ]
    
    # Generate chromosome sequences
    sequences = {
        "chr1": generate_random_sequence(chromosome_length),
        "chr2": generate_random_sequence(chromosome_length),
        "chr3": generate_random_sequence(chromosome_length)
    }
    
    # Introduce duplications with specified identity
    for chrom1, start1, end1, chrom2, start2, end2, identity in duplications:
        seq_to_duplicate = sequences[chrom1][start1:end1]
        
        # Introduce mutations to achieve desired sequence identity
        mutation_rate = 1.0 - identity
        mutated_seq = introduce_mutations(seq_to_duplicate, mutation_rate/2)
        
        # Introduce indels
        mutated_seq = introduce_indels(mutated_seq, mutation_rate/2)
        
        # Insert into target location, handling length differences
        target_len = end2 - start2
        if len(mutated_seq) > target_len:
            mutated_seq = mutated_seq[:target_len]
        
        sequences[chrom2] = sequences[chrom2][:start2] + mutated_seq + sequences[chrom2][start2+len(mutated_seq):]
    
    # Write genome FASTA
    with open(os.path.join(output_dir, "toy_genome.fasta"), "w") as f:
        for chrom, seq in sequences.items():
            f.write(f">{chrom}\n")
            # Write in lines of 60 characters
            for i in range(0, len(seq), 60):
                f.write(f"{seq[i:i+60]}\n")
    
    # Generate methylation data
    methylation_data = {}
    for chrom, seq in sequences.items():
        methylation_data[chrom] = generate_methylation_pattern(len(seq))
    
    # Write methylation BED file
    with open(os.path.join(output_dir, "toy_methylation.bed"), "w") as f:
        for chrom, positions in methylation_data.items():
            for pos, meth_value in positions:
                methylated_reads = int(meth_value * num_reads)
                f.write(f"{chrom}\t{pos}\t{pos+1}\t{methylated_reads}\t{num_reads}\n")
    
    # Write segmental duplication BED file
    with open(os.path.join(output_dir, "toy_segmental_dups.bed"), "w") as f:
        for chrom1, start1, end1, chrom2, start2, end2, _ in duplications:
            f.write(f"{chrom1}\t{start1}\t{end1}\t{chrom2}\t{start2}\t{end2}\n")
    
    # Write README
    with open(os.path.join(output_dir, "README.txt"), "w") as f:
        f.write("""
Realistic Toy Example Input Files for Segmental Duplication Methylation Analysis
===============================================================================

This directory contains realistic example files for testing the segmental duplication 
methylation analysis script.

Files included:
--------------

1. toy_genome.fasta
   - A simulated genome with 3 chromosomes (10,000 bp each)
   - Contains realistic segmental duplications with varying sequence identities

2. toy_methylation.bed
   - Realistic methylation data for CpG sites on all 3 chromosomes
   - Format: chromosome, start, end, methylated reads, total reads
   - Methylation patterns have realistic correlation structures

3. toy_segmental_dups.bed
   - Contains 3 segmental duplication pairs with varying sequence identities
   - Format: chrom1, start1, end1, chrom2, start2, end2
   - Pair 1: ~95% identity
   - Pair 2: ~90% identity
   - Pair 3: ~85% identity

Usage:
------
Run the segmental duplication analysis script with these files:

python segmental_dup_analysis.py \\
  --genome toy_genome.fasta \\
  --methyl toy_methylation.bed \\
  --dups toy_segmental_dups.bed \\
  --window 100 \\
  --threads 2 \\
  --output toy_results.tsv
        """)

    print(f"Realistic toy example files have been created in the {output_dir} directory.")

if __name__ == "__main__":
    create_toy_files("realistic_toy_example")
