#!/usr/bin/env python3
"""
Fixed version of merge_counts.py
Properly merges count files by gene ID, not by position
"""

import os
import sys
import pandas as pd

if len(sys.argv) > 1:
    input_dir = sys.argv[1]
else:
    input_dir = "./"

if len(sys.argv) > 2:
    output_file = sys.argv[2]
else:
    output_file = os.path.join(input_dir if input_dir != "./" else ".", "gene_count_matrix_fixed.tsv")

files = sorted([f for f in os.listdir(input_dir) if f.endswith("_cut.counts")])

if not files:
    print(f"Error: No _cut.counts files found in {input_dir}")
    sys.exit(1)

print(f"Found {len(files)} sample files")
print("="*70)

dfs = []

for f in files:
    file_path = os.path.join(input_dir, f)
    sample_name = f.replace('_cut.counts', '')
    
    print(f"Processing sample: {sample_name}")
    
    # Read file and set Geneid as index - THIS IS THE KEY FIX!
    df = pd.read_csv(file_path, sep='\t', header=0, index_col=0)
    
    # Rename the count column to sample name
    df.columns = [sample_name]
    
    print(f"  - Genes: {len(df)}")
    print(f"  - First gene: {df.index[0]}")
    print(f"  - Last gene: {df.index[-1]}")
    
    dfs.append(df)

print("\n" + "="*70)
print("Merging all samples...")

# Merge by gene ID (index) - outer join to keep all genes
count_matrix = pd.concat(dfs, axis=1, join='outer', sort=True)

# Fill NaN with 0 (if some genes are missing in some samples)
count_matrix = count_matrix.fillna(0).astype(int)

print(f"Merged matrix:")
print(f"  - Total genes: {len(count_matrix)}")
print(f"  - Total samples: {len(count_matrix.columns)}")

# Check for any genes that were in original files but missing after merge
for i, df in enumerate(dfs):
    missing = set(df.index) - set(count_matrix.index)
    if missing:
        print(f"  ⚠️  WARNING: File {files[i]} has {len(missing)} genes not in merged matrix")

# Save the matrix
count_matrix.to_csv(output_file, sep='\t')
print(f"\n✅ Expression matrix saved: {output_file}")
print(f"   - Number of genes: {len(count_matrix)}")
print(f"   - Number of samples: {len(count_matrix.columns)}")
print(f"   - Samples: {', '.join(count_matrix.columns.tolist())}")

# Verify a specific gene as a sanity check
test_genes = ['ATBlh1-1G40660', 'ATBlh1-1G10010', 'ATBlh1-5G53340']
print(f"\n🔍 Sanity check - looking for test genes:")
for gene in test_genes:
    if gene in count_matrix.index:
        print(f"   ✅ {gene} found (counts: {count_matrix.loc[gene].tolist()})")
    else:
        print(f"   ❌ {gene} NOT FOUND")

