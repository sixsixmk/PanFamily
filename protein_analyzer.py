import sys
import argparse
from Bio import SeqIO
from Bio.SeqUtils import ProtParam

def aliphatic_index(sequence):
    """计算脂肪族指数"""
    aliphatic_aa = {'A': 0.283, 'V': 0.606, 'I': 0.530, 'L': 0.603}
    return sum(aliphatic_aa.get(aa, 0) for aa in sequence) / len(sequence) * 100

def analyze_proteins(fasta_file, output_file):
    """分析蛋白质序列的理化性质"""
    records = SeqIO.parse(fasta_file, "fasta")

    with open(output_file, "w") as out_f:
        header = "ID,Protein_Length,Molecular_Weight_kDa,Isoelectric_Point,Hydrophilicity,Aliphatic_Index,Instability_Index\n"
        out_f.write(header)

        for record in records:
            protein_length = len(record.seq)
            protein_analyzer = ProtParam.ProteinAnalysis(str(record.seq))
            molecular_weight = protein_analyzer.molecular_weight() / 1000
            isoelectric_point = protein_analyzer.isoelectric_point()
            hydrophilicity = protein_analyzer.gravy()
            aliphatic_idx = aliphatic_index(record.seq)
            instability_index = protein_analyzer.instability_index()

            line = f"{record.id},{protein_length},{molecular_weight:.4f},{isoelectric_point:.4f},{hydrophilicity:.4f},{aliphatic_idx:.4f},{instability_index:.4f}\n"
            out_f.write(line)

def main():
    parser = argparse.ArgumentParser(description="Analyze protein sequences in a FASTA file and save the results in a CSV file.")
    parser.add_argument("--fasta", required=True, help="Input FASTA file containing protein sequences.")
    parser.add_argument("--csv", required=True, help="Output CSV file to store the results.")

    args = parser.parse_args()
    analyze_proteins(args.fasta, args.csv)

if __name__ == "__main__":
    main()
