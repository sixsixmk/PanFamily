#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Integrate renamed expression matrices from multiple strains
Used for merging expression data from different strains in pan-transcriptome analysis
"""

import os
import sys
import pandas as pd
from pathlib import Path

def merge_matrices(input_dir, output_file, matrix_type="counts"):
    """
    Merge all TSV expression matrix files in the specified directory
    
    Parameters:
        input_dir: Directory containing all strain matrix files
        output_file: Output file path
        matrix_type: Matrix type (counts/fpkm/tpm), for log display
    """
    
    # Get all .tsv files in the directory
    input_path = Path(input_dir)
    tsv_files = sorted(list(input_path.glob("*.tsv")))
    
    if not tsv_files:
        print(f"[ERROR] No .tsv files found in {input_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"[INFO] Found {len(tsv_files)} strain matrix files")
    
    # Read and merge all matrices
    all_dfs = []
    
    for tsv_file in tsv_files:
        strain_name = tsv_file.stem.replace(f"_{matrix_type}", "")
        print(f"  Processing strain: {strain_name} ({tsv_file.name})")
        
        try:
            # Read matrix file
            df = pd.read_csv(tsv_file, sep='\t', index_col=0)
            
            # Add strain prefix to column names (if not already added)
            new_columns = []
            for col in df.columns:
                if not col.startswith(strain_name + "_"):
                    new_columns.append(f"{strain_name}_{col}")
                else:
                    new_columns.append(col)
            df.columns = new_columns
            
            all_dfs.append(df)
            print(f"    - Genes: {len(df)}, Samples: {len(df.columns)}")
            
        except Exception as e:
            print(f"[WARNING] Failed to read file: {tsv_file.name} - {e}", file=sys.stderr)
            continue
    
    if not all_dfs:
        print(f"[ERROR] No matrix files were successfully read", file=sys.stderr)
        sys.exit(1)
    
    # Merge all dataframes (outer join, fill missing values with 0)
    print(f"\n[INFO] Merging matrices...")
    merged_df = all_dfs[0]
    
    for i, df in enumerate(all_dfs[1:], start=2):
        print(f"  Merging {i}/{len(all_dfs)} strains...")
        merged_df = merged_df.join(df, how='outer')
    
    # Fill missing values with 0
    merged_df = merged_df.fillna(0)
    
    # Ensure output directory exists
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save merged matrix
    merged_df.to_csv(output_file, sep='\t')
    
    print(f"Merged matrix saved: {output_file}")
    print(f"  - Total genes: {len(merged_df)}")
    print(f"  - Total samples: {len(merged_df.columns)}")
    print(f"  - Matrix size: {len(merged_df)} × {len(merged_df.columns)}")
    
    return merged_df


def main():
    """Main program"""
    
    if len(sys.argv) < 3:
        print("Usage: merge_renamed_matrices.py <input_dir> <output_file> [matrix_type]")
        print("")
        print("Arguments:")
        print("  input_dir   - Directory containing all renamed strain matrices")
        print("  output_file - Output path for the merged matrix")
        print("  matrix_type - counts/fpkm/tpm (optional, default: counts)")
        print("")
        print("Example:")
        print("  python3 merge_renamed_matrices.py ./rename_matrix/Counts/ ./merged_counts.tsv counts")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    matrix_type = sys.argv[3] if len(sys.argv) > 3 else "counts"
    
    # Check input directory
    if not os.path.isdir(input_dir):
        print(f"[ERROR] Input directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"{'='*70}")
    print(f"  Merge Strain Expression Matrices - {matrix_type.upper()}")
    print(f"{'='*70}")
    print(f"Input directory: {input_dir}")
    print(f"Output file: {output_file}")
    print(f"Matrix type: {matrix_type}")
    print(f"{'='*70}\n")
    
    # Execute merging
    merge_matrices(input_dir, output_file, matrix_type)
    
    print(f"\n{'='*70}")
    print(f"  Merge completed!")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
