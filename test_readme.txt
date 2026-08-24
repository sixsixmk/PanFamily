待补充
================================================================================
PanFamily Installation and Usage Guide
================================================================================

PanFamily is an integrated toolkit for gene family multi-omics analysis, supporting systematic analyses across genomics, pan-genomics, transcriptomics, and single-cell transcriptomics.

================================================================================
I. Installation and Environment Setup
================================================================================

1.1 Runtime Environment Requirements
--------------------------------------------------------------------------------

[Python Environment]
- Python version: 3.8 or above
- Dependencies: os, sys, argparse, BioPython (v1.84)

[R Environment]
- R version: 4.3 or above
- Dependencies:
  optparse (v1.7.5), dplyr (v1.1.4), tidyr (v1.3.1), ggplot2 (v3.5.1),
  ComplexHeatmap (v2.18.0), circlize (v0.4.16), ggtree (v3.10.0),
  ggridges (v0.5.6), viridis (v0.6.5), scales (v1.3.0),
  RColorBrewer (v1.1-3), colorspace (v2.1-1), DESeq2 (v1.42.0)

[External Tools]
- Sequence alignment and search: BLAST+ (v2.15.0), HMMER (v3.4), MAFFT (v7.520)
- Phylogenetic analysis: FastTree (v2.1.11), IQ-TREE (v2.3.0), trimAl (v1.4)
- Sequence processing: seqkit (v2.8.0)
- Selection pressure analysis: KaKs_Calculator (v3.0), ParaAT (v2.0), MUSCLE (v5.1)
- Motif analysis: MEME Suite (v5.5.5)
- Collinearity analysis: MCScanX
- Variant annotation: ANNOVAR
- Transcriptome analysis: featureCounts (v2.0.6)

1.2 Installation Methods
--------------------------------------------------------------------------------

Method 1: Add to system PATH (recommended)
# Download and extract the PanFamily package
# Add to PATH
echo 'export PATH="/path/to/PanFamily:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Verify installation
PanFamily -h

Method 2: Create a symbolic link
sudo ln -s /path/to/PanFamily/PanFamily /usr/local/bin/PanFamily

1.3 View Help Information
--------------------------------------------------------------------------------
PanFamily -h                    # View main help
PanFamily --list-modules        # List all available modules
PanFamily GeneDetection --help  # View help for a specific module

================================================================================
II. Genomics Module (Single-Genome Analysis)
================================================================================

The genomics module is designed for gene family analysis at the single reference genome level, including gene family member identification, phylogenetic analysis, selection pressure analysis, conserved motif analysis, protein physicochemical property analysis, and collinearity analysis.

--------------------------------------------------------------------------------
2.1 GeneDetection - Gene Family Member Identification
--------------------------------------------------------------------------------

Function description:
Identify target gene family members in a single reference genome. Both BLASTP and HMMER methods are supported and can be used independently or in combination.

Usage:
# Method 1: Run via configuration file
PanFamily GeneDetection -c gene_detection.yaml

# Method 2: Run via command-line parameters
PanFamily GeneDetection --method both -i proteins.fa -o results --seed seeds.fa --hmm model.hmm

Configuration parameters (gene_detection.yaml):
--------------------------------------------------------------------------------
Parameter name                Description                                   Example
--------------------------------------------------------------------------------
methods.blastp               Enable BLASTP method                            true/false
methods.hmmer                Enable HMMER method                             true/false
files.input_file             Target species protein FASTA file               /path/to/proteins.fa
files.seed                   Reference gene family protein sequences         /path/to/seeds.fa
files.hmm                    HMM model file path (comma-separated)            /path/to/PF00001.hmm
files.output_dir             Output directory                                 /path/to/output
parameters.threads           Number of threads                                12
parameters.evalue_blast      BLASTP E-value threshold                         1e-5
parameters.evalue_hmm        HMMER E-value threshold                          1e-5
parameters.coverage_threshold HMMER domain coverage threshold                 0.9
--------------------------------------------------------------------------------

Input file formats:
- Protein sequence file: standard FASTA format
- HMM model file: HMMER3 format (downloadable from Pfam)
- Seed sequence file: known gene family protein sequences (FASTA)

Output files:
- gene_family_ids.txt: list of identified gene family member IDs
- gene_family_proteins.fa: protein sequences of gene family members
- blastp_results.txt: BLASTP alignment results
- hmmsearch_results.txt: HMMER search results

Module Tools Download Example
Genome and 
Pangenome
HMMER（v3.4）
yq（v4.2.0）、
Bitacora（v1.4）
seqkit（v2.13）
MUSCLE（v5.1）
ParaAT（v2.0）
MAFFT（v7.526）
trimAl（v1.5.1）
FastTree（v2.1.11）
IQ-TREE（v3.1.2）
ANNOVAR
MEME Suite
KaKs_Calculator
（v3.0）
MCScanX
BLAST+（v2.15.0）
ParaAT（v2.0）
conda install -c bioconda hmmer
conda install -c conda-forge yq
https://github.com/molevol-ub/bitacora
conda install -c bioconda seqkit
conda install -c bioconda muscle
conda install -c bioconda paraat
conda install -c bioconda mafft
conda install -c bioconda trimal
conda install -c bioconda fasttree
conda install -c bioconda iqtree
https://www.openbioinformatics.org/annovar/annovar_downlo
ad_form.php
conda install -c bioconda meme
https://github.com/Chenglin20170390/KaKs_Calculator-3.0 
git clone https://github.com/wyp1125/MCScanX.git
conda install -c bioconda blast
ftp://download.big.ac.cn/bigd/tools/ParaAT2.0.tar.gz
PanFamily GeneDetection -c gene_detection.yaml
PanFamily Analysis Phylogenetic -c 
phylogenetic.yaml
PanFamily Analysis KaKs -c kaks.yaml
PanFamily Analysis motif -c motif.yaml
PanFamily Analysis protein -c protein_analysis.yaml
PanFamily Analysis collinearity -c collinearity.yaml
PanFamily PanDetection -c PanDetection.yaml
PanFamily Analysis pangenelist -c 
PanGeneFamily_list.yaml
PanFamily Analysis PAV -c PanPAV.yaml
PanFamily Analysis PanKaKs -c pankaks.yaml
PanFamily Analysis panmotif -c panmotif.yaml
PanFamily Analysis panprotein -c panprotein.yaml
PanFamily Analysis SV -c sv.yaml
PanFamily Analysis MaptoRef -c maptoref.yaml

Module Tools Download Example
(Pan-) 
Transcriptome
samtools（v1.22.1）
fastqc（v0.12.1）
trim_galore
（v0.6.6）
hisat2（v2.2.1）
featureCounts
（v2.1.1）
conda install -c bioconda samtools
conda install -c bioconda fastqc
conda install -c bioconda trim-galore
conda install -c bioconda hisat2
conda install -c bioconda subread
PanFamily RNAseq -c RNAseq.yaml
PanFamily PanRNAseq -c PanRNAseq.yaml
ScRNA-seq
cellranger（v10.0.0）
Seurat（v5.2.1）
scMayoMap（v0.2.0）
Monocle（v2.32.0）
igraph（v1.2.10）
ggsci（v4.1.0）
harmony（v1.2.4）
https://www.10xgenomics.com/support/software/cellranger/downloads
install.packages("Seurat")
https://github.com/chloelulu/scMayoMap
install.packages("igraph")
install.packages("ggsci")
install.packages("harmony")
PanFamily scRNA -c SC.yaml
Visualization
optparse（v1.7.5）
dplyr（v1.1.4）
tidyr（v1.3.1）
ggplot2（v4.0.0）
ComplexHeatmap
（v2.25.2）
circlize（v0.4.16）
ggtree（v3.12.0）
ggridges（v0.5.7）
install.packages(c("optparse", "dplyr", "tidyr", 
"ggplot2","circlize", "ggridges", "viridis", 
"scales","RColorBrewer", "colorspace"))
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")BiocManager::install(c("Com
plexHeatmap", "ggtree"))
PanFamily -V -i input_dir -o output_dir
