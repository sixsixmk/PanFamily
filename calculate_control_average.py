#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate the average expression of Control group samples
Used for SV pre-analysis
"""

import pandas as pd
import sys
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s', stream=sys.stderr)
logger = logging.getLogger(__name__)

def calculate_control_average(input_matrix, output_file, control_samples_file, strain_name):

    try:
        # 1. Read Control sample list
        logger.info(f"Reading Control sample list: {control_samples_file}")
        with open(control_samples_file, 'r') as f:
            control_samples = [line.strip() for line in f if line.strip()]
        
        if not control_samples:
            logger.error("Control sample list is empty")
            return False
        
        logger.info(f"Number of Control samples: {len(control_samples)}")
        
        # 2. Read input matrix
        logger.info(f"Reading input matrix: {input_matrix}")
        df = pd.read_csv(input_matrix, sep='\t', index_col=0)
        logger.info(f"Matrix dimensions: {df.shape[0]} genes × {df.shape[1]} samples")
        
        # 3. Extract sample names from column names and match Control group
        # Column name format: StrainA_Sample1 -> extract Sample1
        control_cols = []
        for col in df.columns:
            if '_' in col:
                # Split column name to extract sample part (remove strain prefix)
                sample_part = col.split('_', 1)[1]
                if sample_part in control_samples:
                    control_cols.append(col)
        
        if len(control_cols) == 0:
            logger.error(f"No matching Control sample columns found")
            logger.error(f"Expected sample names: {control_samples}")
            logger.error(f"Column names in matrix: {list(df.columns[:5])}... (first 5 columns)")
            return False
        
        logger.info(f"Matched {len(control_cols)} Control sample columns")
        
        # 4. Extract Control group columns and calculate the mean
        df_control = df[control_cols]
        avg_expr = df_control.mean(axis=1)
        
        # 5. Output two-column file
        result = pd.DataFrame({
            'GeneID': avg_expr.index,
            strain_name: avg_expr.values
        })
        
        result.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Success: Output file {output_file}")
        logger.info(f"  Number of genes: {len(result)}")
        logger.info(f"  Number of samples used: {len(control_cols)}")
        
        return True
        
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        return False
    except Exception as e:
        logger.error(f"Error occurred during processing: {e}")
        import traceback
        traceback.print_exc(file=sys.stderr)
        return False

def main():
    parser = argparse.ArgumentParser(
        description='Calculate average expression of Control group samples',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example usage:

  python3 calculate_control_average.py \\
      -i gene_cluster/count/StrainA_count_matrix.tsv \\
      -o sv_pre/count/StrainA_control_count.tsv \\
      -c control_samples.txt \\
      -s StrainA

Format of Control sample list (control_samples.txt):
  Sample1
  Sample2
  Sample3

Input matrix format (TSV):
  GeneID    StrainA_Sample1    StrainA_Sample2    StrainA_Sample3
  Gene1     100                150                120
  Gene2     200                180                210

Output file format (TSV):
  GeneID    StrainA
  Gene1     123.33
  Gene2     196.67
        '''
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help='Input matrix file (TSV format)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file path (two-column TSV)')
    parser.add_argument('-c', '--control_samples', required=True,
                        help='Control sample list file (one sample name per line)')
    parser.add_argument('-s', '--strain', required=True,
                        help='Strain name (used as column name in output file)')
    
    args = parser.parse_args()
    
    # Execute calculation
    success = calculate_control_average(
        args.input,
        args.output,
        args.control_samples,
        args.strain
    )
    
    if success:
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()
