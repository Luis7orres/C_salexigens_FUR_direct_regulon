#!/bin/bash
#SBATCH --job-name=chip06

#### FUR2 SUBSAMPLED ANALYSIS ####

###  SUBSAMPLING ###

## AP 0.6 --> R1 y R2

seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/AP/chip/0.6/SM001943-9_S15_L002_R2_001.fastq.gz 2500000 > /home/omicas/luistorres/TFM/samples/AP/chip/0.6/R2_chip_subsampled.fastq.gz
seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/AP/chip/0.6/SM001943-9_S15_L002_R1_001.fastq.gz 2500000 > /home/omicas/luistorres/TFM/samples/AP/chip/0.6/R1_chip_subsampled.fastq.gz

## INPUT 0.6 --> R1 y R2 

seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/AP/input/0.6/SM001943-3_S9_L002_R1_001.fastq.gz 2500000 > /home/omicas/luistorres/TFM/samples/AP/input/0.6/R1_input_subsampled.fastq.gz
seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/AP/input/0.6/SM001943-3_S9_L002_R2_001.fastq.gz 2500000 > /home/omicas/luistorres/TFM/samples/AP/input/0.6/R2_input_subsampled.fastq.gz

## AP 2.5 --> R1 y R2

seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/AP/chip/2.5/SM001943-12_S18_L002_R1_001.fastq.gz 2500000 > /home/omicas/luistorres/TFM/samples/AP/chip/2.5/R1_chip_subsampled.fastq.gz
seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/AP/chip/2.5/SM001943-12_S18_L002_R2_001.fastq.gz 2500000 > /home/omicas/luistorres/TFM/samples/AP/chip/2.5/R2_chip_subsampled.fastq.gz

## INPUT 2.5 --> R1 y R2

seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/AP/input/2.5/SM001943-6_S12_L002_R1_001.fastq.gz 2500000 > /home/omicas/luistorres/TFM/samples/AP/input/2.5/R1_input_subsampled.fastq.gz
seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/AP/input/2.5/SM001943-6_S12_L002_R2_001.fastq.gz 2500000 > /home/omicas/luistorres/TFM/samples/AP/input/2.5/R2_input_subsampled.fastq.gz


### READ MAPPING 0.6M NaCl ###

## CHIP 

cd /home/omicas/luistorres/TFM/samples/AP/chip/0.6/

bowtie2 -x ../../../../genome/index -1 R1_chip_subsampled.fastq.gz -2 R2_chip_subsampled.fastq.gz -S chip_0.6_AP_subsampled.sam
samtools sort -o chip_0.6_AP_subsampled.bam chip_0.6_AP_subsampled.sam
samtools index chip_0.6_AP_subsampled.bam 

## INPUT

cd /home/omicas/luistorres/TFM/samples/AP/input/0.6/

bowtie2 -x ../../../../genome/index -1 R1_input_subsampled.fastq.gz -2 R2_input_subsampled.fastq.gz -S input_0.6_AP_subsampled.sam
samtools sort -o input_0.6_AP_subsampled.bam input_0.6_AP_subsampled.sam
samtools index input_0.6_AP_subsampled.bam

### READ MAPPING 2.5M NaCl ###

# CHIP

cd /home/omicas/luistorres/TFM/samples/AP/chip/2.5/

bowtie2 -x ../../../../genome/index -1 R1_chip_subsampled.fastq.gz -2 R2_chip_subsampled.fastq.gz -S chip_2.5_AP_subsampled.sam
samtools sort -o chip_2.5_AP_subsampled.bam chip_2.5_AP_subsampled.sam
samtools index chip_2.5_AP_subsampled.bam

# INPUT

cd /home/omicas/luistorres/TFM/samples/AP/input/2.5/

bowtie2 -x ../../../../genome/index -1 R1_input_subsampled.fastq.gz -2 R2_input_subsampled.fastq.gz -S input_2.5_AP_subsampled.sam
samtools sort -o input_2.5_AP_subsampled.bam input_2.5_AP_subsampled.sam
samtools index input_2.5_AP_subsampled.bam


### PEAK CALLING ###

## 0.6M

cd /home/omicas/luistorres/TFM/results/0.6/AP/peaks_submuestreado/

macs2 callpeak -t ../../../../samples/AP/chip/0.6/chip_0.6_AP_subsampled.bam -c ../../../../samples/AP/input/0.6/input_0.6_AP_subsampled.bam -f BAM --outdir . --nomodel --extsize 147 -n AP_peaks_subsampled

## 2.5

cd /home/omicas/luistorres/TFM/results/2.5/AP/peaks_submuestreado/

macs2 callpeak -t ../../../../samples/AP/chip/2.5/chip_2.5_AP_subsampled.bam -c ../../../../samples/AP/input/2.5/input_2.5_AP_subsampled.bam -f BAM --outdir . --nomodel --extsize 147 -n AP_peaks_subsampled_2.5





