# TMF_LuisTorres 
## Iron homeostasis and osmoadaptation: identification of the direct regulon of Fur metalloregulators in Chromohalobacter salexigens
This repository contains the scripts utilized for the ChIP-Seq data analysis of Fur1 and Fur2 mutants under 0.6 and 2.5M NaCl salinity conditions.

## Contents

### 1. *SHELL*
This section includes scripts for read mapping and peak calling for Fur1 and Fur2 mutants under different salinity conditions.

### FULL DATASET ANALYSIS
- FUR1_analysis.sh
- FUR2_analysis.sh
### ANALYISIS WITH SUBSAMPLED READS
- FUR1_subsampled_analysis.sh
- FUR2_subsampled_analysis.sh

**Programs required:** bowtie2, macs2, Seqtk

### 2. *R*
This section includes R scripts for ChIP-Seq peak-gene association, gene-operon association, and correlation of ChIP-Seq data with expression data.

- Peak-gene_association.R
- OPERON.R
- RNASEQ_CHIPSEQ.R

## Contributors
- Luis Torres Ares





