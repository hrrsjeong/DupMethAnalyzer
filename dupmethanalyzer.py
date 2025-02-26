#!/usr/bin/env python3
"""
DupMethAnalyzer

Analyzes DNA methylation differences between segmental duplications in the genome.
Requires minimap2 for alignment.

Author: Harris Jeong
Version: 0.1.0
Date: February 25, 2025
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import subprocess
import tempfile
import os
import concurrent.futures
from tqdm import tqdm
from Bio import SeqIO
import logging
import sys

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Define version
__version__ = "0.1.0"

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='DupMethAnalyzer: Analyze methylation differences between segmental duplications.')
    parser.add_argument('--genome', help='Path to genome FASTA file')
    parser.add_argument('--methyl', help='Path to methylation BED file')
    parser.add_argument('--dups', help='Path to segmental duplications BED file')
    parser.add_argument('--window', type=int, default=100, help='Sliding window size (default: 100bp)')
    parser.add_argument('--output', default='methyl_diff_output.tsv', help='Output file path')
    parser.add_argument('--minimap2', default='minimap2', help='Path to minimap2 executable')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads for parallel processing')
    parser.add_argument('--minimap2-threads', type=int, default=1, 
                        help='Number of threads per minimap2 process (default: 1)')
    parser.add_argument('--keep-alignments', action='store_true', 
                        help='Keep alignment files (PAF format) from minimap2')
    parser.add_argument('--alignment-dir', default='alignments', 
                        help='Directory to store alignment files if --keep-alignments is set')
    parser.add_argument('--windows-bed', default=None, 
                        help='Output BED file containing all analyzed windows coordinates')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}',
                        help='Show program version and exit')
    
    args = parser.parse_args()
    
    # Check required arguments if not showing version
    if '--version' not in sys.argv and '-v' not in sys.argv:
        if not args.genome:
            parser.error("--genome argument is required")
        if not args.methyl:
            parser.error("--methyl argument is required")
        if not args.dups:
            parser.error("--dups argument is required")
    
    return args

def load_genome(genome_file):
    """Load genome sequences from FASTA file."""
    logger.info("Loading genome...")
    genome_dict = {}
    for record in SeqIO.parse(genome_file, "fasta"):
        genome_dict[record.id] = str(record.seq)
    logger.info(f"Loaded {len(genome_dict)} sequences from genome")
    return genome_dict

def load_methylation_data(methyl_file):
    """Load methylation data from BED file."""
    logger.info("Loading methylation data...")
    methyl_data = defaultdict(list)
    count = 0
    with open(methyl_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 5:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                meth_reads, total_reads = int(parts[3]), int(parts[4])
                
                # Calculate methylation level
                meth_level = meth_reads / total_reads if total_reads > 0 else 0
                
                # Store as (position, methylation_level)
                methyl_data[chrom].append((start, meth_level))
                count += 1
    
    # Sort methylation data by position for each chromosome
    for chrom in methyl_data:
        methyl_data[chrom].sort(key=lambda x: x[0])
    
    logger.info(f"Loaded {count} methylation sites across {len(methyl_data)} chromosomes")
    return methyl_data

def load_segmental_duplications(dups_file):
    """Load segmental duplication coordinates from BED file."""
    logger.info("Loading segmental duplications...")
    dup_pairs = []
    with open(dups_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                pair1 = (parts[0], int(parts[1]), int(parts[2]))  # chrom, start, end for pair 1
                pair2 = (parts[3], int(parts[4]), int(parts[5]))  # chrom, start, end for pair 2
                dup_pairs.append((pair1, pair2))
    logger.info(f"Loaded {len(dup_pairs)} segmental duplication pairs")
    return dup_pairs

def get_sequence(genome_dict, chrom, start, end):
    """Extract sequence from genome given coordinates."""
    if chrom in genome_dict:
        return genome_dict[chrom][start:end]
    return None

def run_minimap2(seq1, seq2, minimap2_path, threads=1, keep_alignment=False, alignment_dir=None, 
              alignment_prefix=None, return_mapping=False):
    """
    Align two sequences using minimap2 and return alignment statistics.
    
    Parameters:
    -----------
    seq1, seq2 : str
        Sequences to align
    minimap2_path : str
        Path to minimap2 executable
    threads : int
        Number of threads for minimap2
    keep_alignment : bool
        Whether to keep the alignment file
    alignment_dir : str
        Directory to save alignment files
    alignment_prefix : str
        Prefix for alignment filename
    return_mapping : bool
        Whether to return coordinate mapping between sequences
        
    Returns:
    --------
    If not return_mapping:
        float or (float, str): Sequence identity or (identity, paf_path)
    If return_mapping:
        (float, dict) or (float, dict, str): (identity, coordinate_mapping) or 
                                            (identity, coordinate_mapping, paf_path)
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as seq1_file:
        seq1_file.write(">seq1\n" + seq1 + "\n")
        seq1_path = seq1_file.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as seq2_file:
        seq2_file.write(">seq2\n" + seq2 + "\n")
        seq2_path = seq2_file.name
    
    # Create alignment directory if keeping alignments
    if keep_alignment and alignment_dir:
        os.makedirs(alignment_dir, exist_ok=True)
        if alignment_prefix:
            paf_output = os.path.join(alignment_dir, f"{alignment_prefix}.paf")
        else:
            # Create a unique filename if no prefix provided
            paf_output = os.path.join(alignment_dir, f"alignment_{os.path.basename(tempfile.mktemp())}.paf")
    else:
        paf_output = tempfile.NamedTemporaryFile(suffix='.paf', delete=False).name
    
    try:
        # Run minimap2 with parameters suitable for highly similar sequences
        # -x asm5 is suitable for sequences with ~95% identity (typical for segmental duplications)
        # -c outputs CIGAR string which is needed for detailed alignment analysis
        # -t sets the number of threads
        cmd = [minimap2_path, "-x", "asm5", "-c", "--eqx", "-t", str(threads), seq2_path, seq1_path]
        with open(paf_output, 'w') as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
        
        # Parse PAF output
        identity = 0.0
        matches = 0
        total_aligned = 0
        
        # Create coordinate mapping if requested
        coordinate_mapping = {}
        best_alignment = None
        best_alignment_score = -1
        
        with open(paf_output, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 12:
                    # PAF format: 
                    # query_name, query_length, query_start, query_end, strand, target_name, 
                    # target_length, target_start, target_end, num_matches, alignment_length, mapping_quality, ...
                    query_name = fields[0]
                    query_length = int(fields[1])
                    query_start = int(fields[2])
                    query_end = int(fields[3])
                    strand = fields[4]
                    target_name = fields[5]
                    target_length = int(fields[6])
                    target_start = int(fields[7])
                    target_end = int(fields[8])
                    num_matches = int(fields[9])
                    alignment_length = int(fields[10])
                    mapping_quality = int(fields[11])
                    
                    # Keep track of the best alignment (highest number of matches)
                    current_score = num_matches
                    if current_score > best_alignment_score:
                        best_alignment_score = current_score
                        best_alignment = {
                            'query_start': query_start,
                            'query_end': query_end,
                            'target_start': target_start,
                            'target_end': target_end,
                            'strand': strand,
                            'num_matches': num_matches,
                            'alignment_length': alignment_length
                        }
                    
                    # Parse CIGAR string for detailed coordinate mapping
                    cigar_str = None
                    for field in fields[12:]:
                        if field.startswith('cg:Z:'):
                            cigar_str = field[5:]
                            break
                    
                    # Process CIGAR string to build coordinate mapping
                    if return_mapping and cigar_str:
                        import re
                        cigar_parts = re.findall(r'(\d+)([MIDNSHPX=])', cigar_str)
                        
                        # Current positions in query (seq2) and target (seq1)
                        q_pos = query_start
                        t_pos = target_start
                        
                        for length, op in cigar_parts:
                            length = int(length)
                            
                            if op == '=' or op == 'M':  # Match or mismatch
                                # Create mapping for each position in this segment
                                for i in range(length):
                                    coordinate_mapping[q_pos + i] = t_pos + i
                                
                                # Update positions
                                q_pos += length
                                t_pos += length
                                
                                # Count matches
                                matches += length
                                total_aligned += length
                                
                            elif op == 'I':  # Insertion in query
                                q_pos += length
                                total_aligned += length
                                
                            elif op == 'D':  # Deletion in query (insertion in target)
                                t_pos += length
                                total_aligned += length
                            
                            elif op == 'X':  # Mismatch
                                # Create mapping for mismatched positions too
                                for i in range(length):
                                    coordinate_mapping[q_pos + i] = t_pos + i
                                
                                q_pos += length
                                t_pos += length
                                total_aligned += length
                    
                    # If CIGAR parsing failed, use other fields
                    if total_aligned == 0:
                        matches = num_matches
                        total_aligned = alignment_length
        
        if total_aligned > 0:
            identity = matches / total_aligned
        
        # Return appropriate results based on parameters
        if return_mapping:
            if keep_alignment:
                return identity, coordinate_mapping, paf_output
            else:
                return identity, coordinate_mapping, best_alignment
        else:
            if keep_alignment:
                return identity, paf_output
            else:
                return identity
    
    finally:
        # Clean up temporary files unless we're keeping the alignment
        files_to_remove = [seq1_path, seq2_path]
        if not keep_alignment:
            files_to_remove.append(paf_output)
            
        for file_path in files_to_remove:
            try:
                os.unlink(file_path)
            except:
                pass

def get_methylation_level_in_region(methyl_data, chrom, start, end):
    """Calculate average methylation level in a genomic region."""
    if chrom not in methyl_data:
        return None
    
    # Find all methylation sites within the region
    region_meth = [m[1] for m in methyl_data[chrom] if start <= m[0] < end]
    
    if not region_meth:
        return None
    
    return sum(region_meth) / len(region_meth)

def process_pair(pair_data):
    """Process a single duplication pair."""
    i, ((chrom1, start1, end1), (chrom2, start2, end2)), genome_dict, methyl_data, window_size, minimap2_path, minimap2_threads, keep_alignments, alignment_dir = pair_data
    
    pair_results = []
    
    # Get sequences
    seq1 = get_sequence(genome_dict, chrom1, start1, end1)
    seq2 = get_sequence(genome_dict, chrom2, start2, end2)
    
    if not seq1 or not seq2:
        logger.warning(f"Could not retrieve sequences for pair {i+1}")
        return pair_results
    
    # Calculate overall sequence identity using minimap2 and get the coordinate mapping
    if keep_alignments:
        # Create alignment prefix for the whole pair
        pair_prefix = f"pair_{i+1}_overall"
        overall_identity, coord_mapping, overall_paf = run_minimap2(
            seq1, seq2, minimap2_path, threads=minimap2_threads, return_mapping=True,
            keep_alignment=True, alignment_dir=alignment_dir, alignment_prefix=pair_prefix
        )
    else:
        overall_identity, coord_mapping, best_alignment = run_minimap2(
            seq1, seq2, minimap2_path, threads=minimap2_threads, return_mapping=True
        )
    
    # Length of both sequences
    len1, len2 = len(seq1), len(seq2)
    
    # Extract alignment information
    if not keep_alignments and best_alignment:
        query_start = best_alignment['query_start']  # seq2 start
        query_end = best_alignment['query_end']      # seq2 end
        target_start = best_alignment['target_start']  # seq1 start
        target_end = best_alignment['target_end']      # seq1 end
        strand = best_alignment['strand']
    
    # Analyze methylation in sliding windows
    for w_start1 in range(0, len1 - window_size + 1, window_size // 2):  # 50% overlap
        w_end1 = w_start1 + window_size
        if w_end1 > len1:
            w_end1 = len1
        
        # Calculate proportion of window size
        window_prop = (w_end1 - w_start1) / window_size
        
        # Only process windows that are at least 80% of the desired window size
        if window_prop < 0.8:
            continue
        
        # Get window sequence for pair 1
        window_seq1 = seq1[w_start1:w_end1]
        
        # Find the corresponding window in pair 2 using the coordinate mapping
        mapped_positions = []
        
        # Collect all mapped positions for this window
        for pos1 in range(w_start1, w_end1):
            # Find positions in seq2 that map to this position in seq1
            # We need to reverse the mapping (query→target to target→query)
            for pos2, mapped_pos1 in coord_mapping.items():
                if mapped_pos1 == pos1:
                    mapped_positions.append(pos2)
        
        # If we couldn't find mappings using the coordinate map, use the alignment boundaries
        if not mapped_positions and best_alignment:
            # Calculate relative position within the aligned region
            if w_start1 >= target_start and w_start1 < target_end:
                rel_start = (w_start1 - target_start) / (target_end - target_start)
                rel_end = (min(w_end1, target_end) - target_start) / (target_end - target_start)
                
                # Map to query coordinates
                w_start2 = int(query_start + rel_start * (query_end - query_start))
                w_end2 = int(query_start + rel_end * (query_end - query_start))
                
                # Check if window is out of bounds
                if w_end2 > len2:
                    w_end2 = len2
            else:
                # Window is outside the aligned region, skip it
                continue
        elif mapped_positions:
            # Use the minimum and maximum mapped positions to define the window in seq2
            w_start2 = min(mapped_positions)
            w_end2 = max(mapped_positions) + 1  # +1 to make it inclusive
            
            # Check if window is too small
            if w_end2 - w_start2 < window_size * 0.8:
                # Try to expand the window
                expansion = (window_size - (w_end2 - w_start2)) // 2
                w_start2 = max(0, w_start2 - expansion)
                w_end2 = min(len2, w_end2 + expansion)
        else:
            # No mapping information available, use approximation based on sequence length
            logger.warning(f"No mapping found for window {w_start1}-{w_end1} in pair {i+1}, using approximation")
            position_ratio = len2 / len1
            w_start2 = int(w_start1 * position_ratio)
            w_end2 = int(w_end1 * position_ratio)
            
            if w_end2 > len2:
                w_end2 = len2
        
        # Get window sequence for pair 2
        window_seq2 = seq2[w_start2:w_end2]
        
        # Skip if window is too small in either sequence
        if len(window_seq1) < window_size * 0.8 or len(window_seq2) < window_size * 0.8:
            continue
        
        # Window-specific alignment file
        if keep_alignments:
            # Create alignment prefix for this window
            window_prefix = f"pair_{i+1}_window_{w_start1}-{w_end1}"
            window_identity, window_paf = run_minimap2(
                window_seq1, window_seq2, minimap2_path, threads=minimap2_threads,
                keep_alignment=True, alignment_dir=alignment_dir, alignment_prefix=window_prefix
            )
        else:
            window_identity = run_minimap2(window_seq1, window_seq2, minimap2_path, threads=minimap2_threads)
        
        # Calculate methylation levels
        meth_level1 = get_methylation_level_in_region(methyl_data, chrom1, start1 + w_start1, start1 + w_end1)
        meth_level2 = get_methylation_level_in_region(methyl_data, chrom2, start2 + w_start2, start2 + w_end2)
        
        # Calculate methylation difference
        meth_diff = abs(meth_level1 - meth_level2) if meth_level1 is not None and meth_level2 is not None else None
        
        # Create unique window ID for tracking
        window_id = f"pair{i+1}_w{w_start1}-{w_end1}"
        
        # Store results
        result_entry = {
            'window_id': window_id,
            'dup_pair': i + 1,
            'chrom1': chrom1,
            'start1': start1 + w_start1,
            'end1': start1 + w_end1,
            'chrom2': chrom2,
            'start2': start2 + w_start2,
            'end2': start2 + w_end2,
            'overall_identity': overall_identity,
            'window_identity': window_identity,
            'meth_level1': meth_level1,
            'meth_level2': meth_level2,
            'meth_diff': meth_diff,
            'mapped_positions_count': len(mapped_positions) if mapped_positions else 0
        }
        
        # Add alignment file paths if keeping alignments
        if keep_alignments:
            result_entry['window_alignment'] = window_paf
        
        pair_results.append(result_entry)
    
    return pair_results

def analyze_duplicated_pairs_parallel(genome_dict, methyl_data, dup_pairs, window_size, minimap2_path, 
                             threads, minimap2_threads, keep_alignments=False, alignment_dir=None):
    """Analyze methylation differences between duplicated pairs using parallel processing."""
    logger.info(f"Starting parallel analysis with {threads} worker threads")
    
    all_results = []
    
    # Prepare input data for parallel processing
    input_data = [
        (i, dup_pair, genome_dict, methyl_data, window_size, minimap2_path, minimap2_threads, keep_alignments, alignment_dir) 
        for i, dup_pair in enumerate(dup_pairs)
    ]
    
    # Process pairs in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        # Use tqdm to show progress
        futures = list(tqdm(executor.map(process_pair, input_data), total=len(input_data), desc="Processing duplicate pairs"))
        
        # Collect results
        for result in futures:
            all_results.extend(result)
    
    logger.info(f"Completed analysis of {len(dup_pairs)} pairs, generated {len(all_results)} window results")
    return all_results

def write_windows_bed(results_df, output_file):
    """Write window coordinates to BED file."""
    with open(output_file, 'w') as f:
        # Write header
        f.write("# BED file containing all analyzed windows\n")
        f.write("# Format: chrom start end name score strand\n")
        
        # Write data for pair 1 windows
        for _, row in results_df.iterrows():
            # Create name with window ID and identity/methylation info
            name = f"{row['window_id']}_id={row['window_identity']:.3f}_meth={row['meth_level1']:.3f}"
            # Use methylation difference as score (scaled 0-1000)
            score = int(row['meth_diff'] * 1000) if pd.notnull(row['meth_diff']) else 0
            # Write BED line for pair 1
            f.write(f"{row['chrom1']}\t{row['start1']}\t{row['end1']}\t{name}\t{score}\t+\n")
            
            # Write BED line for pair 2
            name2 = f"{row['window_id']}_pair2_id={row['window_identity']:.3f}_meth={row['meth_level2']:.3f}"
            f.write(f"{row['chrom2']}\t{row['start2']}\t{row['end2']}\t{name2}\t{score}\t-\n")
    
    logger.info(f"Window coordinates saved to BED file: {output_file}")

def main():
    args = parse_arguments()
    
    # Display banner
    logger.info(f"DupMethAnalyzer v{__version__}")
    logger.info(f"================================================")
    
    # Check if minimap2 is available
    try:
        minimap2_version = subprocess.run(
            [args.minimap2, "--version"], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            check=True, 
            text=True
        ).stdout.strip()
        logger.info(f"Using minimap2: {minimap2_version}")
    except (subprocess.SubprocessError, FileNotFoundError):
        logger.error(f"Could not run minimap2. Please ensure it is installed and provide the correct path.")
        return
    
    # Create alignment directory if keeping alignments
    if args.keep_alignments:
        logger.info(f"Alignment files will be kept in directory: {args.alignment_dir}")
        os.makedirs(args.alignment_dir, exist_ok=True)
    
    # Set maximum threads for minimap2 to avoid oversubscription
    minimap2_threads = min(args.minimap2_threads, 4)  # Limit minimap2 threads to avoid excessive resource usage
    logger.info(f"Using {args.threads} worker processes with {minimap2_threads} threads per minimap2 process")
    
    # Log analysis parameters
    logger.info(f"Window size: {args.window}bp")
    logger.info(f"Output file: {args.output}")
    if args.windows_bed:
        logger.info(f"Windows BED file: {args.windows_bed}")
    
    # Load data
    genome_dict = load_genome(args.genome)
    methyl_data = load_methylation_data(args.methyl)
    dup_pairs = load_segmental_duplications(args.dups)
    
    logger.info(f"Starting analysis of {len(dup_pairs)} segmental duplication pairs")
    
    # Analyze duplicated pairs using parallel processing
    results = analyze_duplicated_pairs_parallel(
        genome_dict, methyl_data, dup_pairs, args.window, args.minimap2, 
        args.threads, minimap2_threads, args.keep_alignments, args.alignment_dir
    )
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    if results_df.empty:
        logger.warning("No results generated. Check input files and parameters.")
        return
    
    # Save results to file
    results_df.to_csv(args.output, sep='\t', index=False)
    logger.info(f"Results saved to {args.output}")
    
    # Write window coordinates to BED file if requested
    if args.windows_bed:
        write_windows_bed(results_df, args.windows_bed)
    
    # Generate summary plots
    if not results_df.empty:
        # Plot 1: Sequence Identity vs Methylation Difference
        plt.figure(figsize=(10, 6))
        valid_data = results_df.dropna(subset=['window_identity', 'meth_diff'])
        
        if len(valid_data) > 0:
            plt.scatter(valid_data['window_identity'], valid_data['meth_diff'], alpha=0.6)
            plt.xlabel('Sequence Identity')
            plt.ylabel('Methylation Difference')
            plt.title('Relationship between Sequence Identity and Methylation Difference')
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.savefig(args.output.replace('.tsv', '_identity_vs_diff.png'))
            logger.info(f"Created sequence identity vs methylation difference plot")
        else:
            logger.warning("No valid data for sequence identity vs methylation difference plot")
        
        # Plot 2: Methylation Levels Comparison
        plt.figure(figsize=(10, 6))
        valid_data = results_df.dropna(subset=['meth_level1', 'meth_level2'])
        
        if len(valid_data) > 0:
            plt.scatter(valid_data['meth_level1'], valid_data['meth_level2'], alpha=0.6)
            max_val = max(valid_data['meth_level1'].max(), valid_data['meth_level2'].max())
            plt.plot([0, max_val], [0, max_val], 'r--')  # Diagonal line
            plt.xlabel('Methylation Level (Pair 1)')
            plt.ylabel('Methylation Level (Pair 2)')
            plt.title('Comparison of Methylation Levels between Duplicated Pairs')
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.savefig(args.output.replace('.tsv', '_meth_comparison.png'))
            logger.info(f"Created methylation levels comparison plot")
        else:
            logger.warning("No valid data for methylation levels comparison plot")
        
        # Plot 3: Histogram of sequence identity distribution
        plt.figure(figsize=(10, 6))
        if 'window_identity' in results_df.columns and not results_df['window_identity'].isna().all():
            plt.hist(results_df['window_identity'], bins=20, alpha=0.7)
            plt.xlabel('Sequence Identity')
            plt.ylabel('Frequency')
            plt.title('Distribution of Sequence Identity in Segmental Duplications')
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.savefig(args.output.replace('.tsv', '_identity_dist.png'))
            logger.info(f"Created sequence identity distribution plot")
        else:
            logger.warning("No valid data for sequence identity distribution plot")
    
    logger.info("Analysis completed successfully")
    
if __name__ == "__main__":
    main()
