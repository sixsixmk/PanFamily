#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import argparse
from pathlib import Path


def load_id_map(idmap_file):
    """
    Load ID mapping file
    
    Args:
        idmap_file: Path to ID mapping file (two columns: OriginalID, NewID)
    
    Returns:
        Dictionary: {OriginalID: NewID}
    """
    if not os.path.isfile(idmap_file):
        print(f"[ERROR] ID mapping file does not exist: {idmap_file}", file=sys.stderr)
        sys.exit(1)
    
    id_map = {}
    unmapped_count = 0
    
    with open(idmap_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            # Skip empty lines and comment lines
            if not line or line.startswith('#'):
                continue
            
            # Skip header
            if line_num == 1 and '\t' in line:
                parts = line.split('\t')
                if 'ID' in parts[0] or 'Original' in parts[0]:
                    continue
            
            parts = line.split('\t')
            if len(parts) >= 2:
                original_id = parts[0].strip()
                new_id = parts[1].strip()
                
                if original_id and new_id:
                    id_map[original_id] = new_id
            else:
                unmapped_count += 1
    
    print(f"[INFO] Loaded {len(id_map)} ID mappings")
    if unmapped_count > 0:
        print(f"[INFO] Skipped {unmapped_count} invalid mappings")
    
    return id_map


def rename_matrix(input_file, output_file, id_map, strain_name):
    """
    Rename gene IDs and sample names in the expression matrix
    
    Args:
        input_file: Input matrix file
        output_file: Output matrix file
        id_map: ID mapping dictionary
        strain_name: Strain name (used for renaming sample columns)
    """
    print(f"\n[Processing] {os.path.basename(input_file)}")
    
    # Read matrix
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    original_gene_count = len(df)
    original_sample_count = len(df.columns)
    
    print(f"  Original matrix: {original_gene_count} genes × {original_sample_count} samples")
    
    # Rename gene IDs
    rename_dict = {}
    unmapped_genes = []
    
    for gene_id in df.index:
        if gene_id in id_map:
            new_id = id_map[gene_id]
            rename_dict[gene_id] = new_id
        else:
            # Keep unmapped IDs unchanged
            rename_dict[gene_id] = gene_id
            unmapped_genes.append(gene_id)
    
    # Apply renaming
    df.index = df.index.map(rename_dict)
    
    # Handle many-to-one mapping: average expression values
    if len(df.index) != len(df.index.unique()):
        print(f"  Detected many-to-one mapping, calculating average expression...")
        df = df.groupby(df.index).mean()
    
    # Rename sample columns
    if strain_name:
        new_columns = []
        for col in df.columns:
            if not col.startswith(strain_name + "_"):
                new_columns.append(f"{strain_name}_{col}")
            else:
                new_columns.append(col)
        df.columns = new_columns
    
    # Save results
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    df.to_csv(output_file, sep='\t')
    
    print(f"  After renaming: {len(df)} genes × {len(df.columns)} samples")
    print(f"  ✓ Saved: {output_file}")
    
    if unmapped_genes and len(unmapped_genes) <= 10:
        print(f"  [Note] {len(unmapped_genes)} genes were not mapped (kept original IDs)")


def rename_deg_file(input_file, output_file, id_map):
    """
    Rename gene IDs in DEG result file
    
    Args:
        input_file: Input DEG file (CSV format)
        output_file: Output DEG file
        id_map: ID mapping dictionary
    """
    print(f"\n[Processing] {os.path.basename(input_file)}")
    
    # Read DEG file
    df = pd.read_csv(input_file, index_col=0)
    original_count = len(df)
    
    print(f"  Original DEG gene count: {original_count}")
    
    # Rename gene IDs
    rename_dict = {}
    unmapped_genes = []
    
    for gene_id in df.index:
        if gene_id in id_map:
            new_id = id_map[gene_id]
            rename_dict[gene_id] = new_id
        else:
            rename_dict[gene_id] = gene_id
            unmapped_genes.append(gene_id)
    
    # Apply renaming
    df.index = df.index.map(rename_dict)
    
    # Handle many-to-one mapping: keep the first entry
    if len(df.index) != len(df.index.unique()):
        print(f"  Detected many-to-one mapping, keeping the first entry...")
        df = df[~df.index.duplicated(keep='first')]
    
    # Save results
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    df.to_csv(output_file)
    
    print(f"  After renaming: {len(df)} DEG genes")
    print(f"  ✓ Saved: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Gene ID renaming script',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:

1. Rename expression matrix:
   python3 rename_gene_ids.py -i Blh-1_count_matrix.tsv -o renamed_count.tsv -m idmap.txt -s Blh-1 -t matrix

2. Rename DEG results:
   python3 rename_gene_ids.py -i DEG_up.csv -o DEG_up_renamed.csv -m idmap.txt -t deg

ID mapping file format (two columns, tab-separated):
OriginalID    NewID
AT1G01010     TPS1
AT2G02020     TPS2
gene12345     TPS3

Notes:
- When multiple original IDs map to one new ID, the expression matrix will use the average value.
- DEG result files will keep the first entry.
- Unmapped IDs remain unchanged.
'''
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help='Input file path')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file path')
    parser.add_argument('-m', '--idmap', required=True,
                        help='ID mapping file path')
    parser.add_argument('-s', '--strain', default='',
                        help='Strain name (used for renaming sample columns, only valid for matrix type)')
    parser.add_argument('-t', '--type', choices=['matrix', 'deg'], default='matrix',
                        help='File type: matrix (expression matrix) or deg (DEG results)')
    
    args = parser.parse_args()
    
    if not os.path.isfile(args.input):
        print(f"[ERROR] Input file does not exist: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    print(f"{'='*70}")
    print(f"  Gene ID Renaming")
    print(f"{'='*70}")
    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    print(f"ID mapping file: {args.idmap}")
    print(f"File type: {args.type}")
    if args.strain:
        print(f"Strain name: {args.strain}")
    print(f"{'='*70}")
    
    # Load ID map
    id_map = load_id_map(args.idmap)
    
    try:
        if args.type == 'matrix':
            rename_matrix(args.input, args.output, id_map, args.strain)
        elif args.type == 'deg':
            rename_deg_file(args.input, args.output, id_map)
        
        print(f"\n{'='*70}")
        print(f"  Renaming Completed!")
        print(f"{'='*70}\n")
        
    except Exception as e:
        print(f"\n[ERROR] Processing failed: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
