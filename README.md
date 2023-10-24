# TMF_LuisTorres 
## Iron homeostasis and osmoadaptation: identification of the direct regulon of Fur metalloregulators in Chromohalobacter salexigens
Code used for ChIP-Seq analysis of Fur1 and Fur2 mutants during master thesis.

The analyses were divided in two parts

## 1. *SHELL*
Read mapping and peak calling for Fur1 and Fur2 under 0.6M and 2.5 ChIP-Seq reads 

### FULL DATASET ANALYSIS
- FUR1_analysis.sh
- FUR2_analysis.sh
### ANALYISIS WITH SUBSAMPLED READS
- FUR1_subsampled_analysis.sh
- FUR2_subsampled_analysis.sh

Programs required: bowtie2, macs2, Seqtk

## 2. *R-Studio*
### ChIP-Seq peaks association with genes
- Peak-gene_association.R
### Gene-operon association
- OPERON.R
### Correlation of ChIP-seq Data with Expression Data
- RNASEQ_CHIPSEQ.R





