#!/usr/bin/env python3
"""
Calculate FPKM and TPM from featureCounts count matrix.
Supports extracting gene lengths from a GTF file.
"""

import sys
import os
import pandas as pd
import argparse
import logging

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def extract_gene_lengths_from_gtf(gtf_file, output_file=None):
    """Extract gene lengths (merge all exons) from a GTF file"""
    logger.info(f"Extracting gene lengths from GTF file: {gtf_file}")
    
    gene_exons = {}
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            if feature_type != 'exon':
                continue
            
            start = int(parts[3])
            end = int(parts[4])
            attributes = parts[8]
            
            # Extract gene_id
            gene_id = None
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1]
                    break
            
            if gene_id:
                if gene_id not in gene_exons:
                    gene_exons[gene_id] = []
                gene_exons[gene_id].append((start, end))
    
    # Merge overlapping exons and calculate total gene length
    gene_lengths = {}
    for gene_id, exons in gene_exons.items():
        exons.sort()
        merged_exons = []
        
        for start, end in exons:
            if not merged_exons or merged_exons[-1][1] < start:
                merged_exons.append((start, end))
            else:
                merged_exons[-1] = (merged_exons[-1][0], max(merged_exons[-1][1], end))
        
        total_length = sum(end - start + 1 for start, end in merged_exons)
        gene_lengths[gene_id] = total_length
    
    logger.info(f"Extracted gene lengths for {len(gene_lengths)} genes")
    
    # Save to file (if specified)
    if output_file:
        with open(output_file, 'w') as f:
            f.write("Geneid\tLength\n")
            for gene_id, length in sorted(gene_lengths.items()):
                f.write(f"{gene_id}\t{length}\n")
        logger.info(f"Gene lengths saved to: {output_file}")
    
    return gene_lengths


def load_gene_lengths(length_file):
    """Load gene lengths from a file"""
    logger.info(f"Loading gene length file: {length_file}")
    
    lengths = {}
    with open(length_file, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                gene, length = parts[0], int(parts[1])
                lengths[gene] = length
    
    logger.info(f"Loaded gene lengths for {len(lengths)} genes")
    return lengths


def calculate_fpkm_tpm(count_matrix_file, gene_lengths, output_prefix="output"):
    """
    Calculate FPKM and TPM from a count matrix
    
    Parameters:
        count_matrix_file: Count matrix file (first column: Geneid, others: samples)
        gene_lengths: Dictionary of gene lengths {gene_id: length}
        output_prefix: Output file prefix
    """
    logger.info(f"Calculating FPKM and TPM: {count_matrix_file}")
    
    # Read count matrix
    df = pd.read_csv(count_matrix_file, sep='\t', index_col=0)
    logger.info(f"Loaded matrix: {len(df)} genes × {len(df.columns)} samples")
    
    # Get sample names
    samples = df.columns.tolist()
    
    # Initialize output DataFrames
    fpkm_df = pd.DataFrame(index=df.index)
    tpm_df = pd.DataFrame(index=df.index)
    
    # Track genes without length information
    missing_lengths = []
    
    # Process each sample separately
    for sample in samples:
        logger.info(f"  Processing sample: {sample}")
        counts = df[sample]
        
        # Calculate RPK (Reads Per Kilobase)
        rpk = {}
        for gene in counts.index:
            if gene not in gene_lengths:
                if gene not in missing_lengths:
                    missing_lengths.append(gene)
                continue
            
            gene_len = gene_lengths[gene]
            if gene_len <= 0:
                continue
            
            # RPK = count / (gene_length / 1000)
            rpk[gene] = counts[gene] / (gene_len / 1000.0)
        
        # Calculate total RPK and total counts
        total_rpk = sum(rpk.values())
        total_counts = counts.sum()
        
        if total_rpk == 0 or total_counts == 0:
            logger.warning(f"    Sample {sample}: total RPK or total counts = 0, skipping")
            continue
        
        # Calculate FPKM and TPM
        fpkm_values = {}
        tpm_values = {}
        
        for gene in rpk:
            # FPKM = RPK / (total_reads / 1,000,000)
            fpkm_values[gene] = rpk[gene] / (total_counts / 1e6)
            
            # TPM = (RPK / total_RPK) * 1,000,000
            tpm_values[gene] = (rpk[gene] / total_rpk) * 1e6
        
        fpkm_df[sample] = pd.Series(fpkm_values)
        tpm_df[sample] = pd.Series(tpm_values)
    
    # Warn about missing gene lengths
    if missing_lengths:
        logger.warning(f"{len(missing_lengths)} genes lack length information and were skipped")
        if len(missing_lengths) <= 10:
            logger.warning(f"Missing gene lengths: {', '.join(missing_lengths)}")
    
    # Save results
    fpkm_file = f"{output_prefix}_fpkm_matrix.tsv"
    tpm_file = f"{output_prefix}_tpm_matrix.tsv"
    
    fpkm_df.to_csv(fpkm_file, sep='\t')
    tpm_df.to_csv(tpm_file, sep='\t')
    
    logger.info(f"FPKM file saved: {fpkm_file}")
    logger.info(f"TPM file saved: {tpm_file}")
    logger.info(f"  - Number of genes: {len(fpkm_df)}")
    logger.info(f"  - Number of samples: {len(fpkm_df.columns)}")
    
    return fpkm_df, tpm_df


def main():
    parser = argparse.ArgumentParser(
        description='Calculate FPKM and TPM from a count matrix',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:

1. Extract gene lengths from GTF and calculate FPKM/TPM:
   python3 calculate_fpkm_tpm.py -c gene_count_matrix.tsv -g annotation.gtf -o results

2. Use an existing gene length file:
   python3 calculate_fpkm_tpm.py -c gene_count_matrix.tsv -l gene_lengths.txt -o results

3. Only extract gene lengths:
   python3 calculate_fpkm_tpm.py --extract-length annotation.gtf -o gene_lengths.txt
        '''
    )
    
    parser.add_argument('-c', '--count', help='Count matrix file (TSV format)')
    parser.add_argument('-g', '--gtf', help='GTF annotation file (used to extract gene lengths)')
    parser.add_argument('-l', '--length', help='Gene length file (two columns: Geneid\\tLength)')
    parser.add_argument('-o', '--output', default='output', help='Output file prefix (default: output)')
    parser.add_argument('--extract-length', help='Only extract gene lengths from GTF (skip FPKM/TPM calculation)')
    
    args = parser.parse_args()
    
    # Mode 1: only extract gene lengths
    if args.extract_length:
        extract_gene_lengths_from_gtf(args.extract_length, args.output)
        return
    
    # Mode 2: calculate FPKM and TPM
    if not args.count:
        parser.error("Count matrix file is required (-c/--count)")
    
    # Get gene lengths
    gene_lengths = None
    if args.gtf:
        logger.info("Extracting gene lengths from GTF file...")
        gene_lengths = extract_gene_lengths_from_gtf(args.gtf)
    elif args.length:
        gene_lengths = load_gene_lengths(args.length)
    else:
        parser.error("You must provide either a GTF file (-g) or a gene length file (-l)")
    
    # Calculate FPKM and TPM
    calculate_fpkm_tpm(args.count, gene_lengths, args.output)
    
    logger.info("Calculation completed!")


if __name__ == '__main__':
    main()
