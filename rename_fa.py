#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def load_mapping(mapping_file):
    mapping = {}
    with open(mapping_file) as f:
        for i, line in enumerate(f):
            if i == 0:  
                continue
            parts = line.strip().split('\t') 
            if len(parts) >= 2:
                old_id, new_id = parts[0], parts[1]
                mapping[old_id] = new_id
    return mapping

def replace_ids(input_fasta, mapping_file, output_fasta):
    mapping = load_mapping(mapping_file)
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in mapping:
                new_id = mapping[record.id]
                record.id = new_id
                record.description = new_id 
            SeqIO.write(record, out_f, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Replace FASTA IDs using a mapping file")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-m", "--mapping", required=True, help="Mapping file (tab-delimited, with header)")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    replace_ids(args.input, args.mapping, args.output)

