#!/bin/bash
#SBATCH --job-name=FUR1sub

#### FUR1 0.6M NaCl SUBSAMPLED ANALYSIS ####

### SUBSAMPLING ###

## AP 0.6 --> R1 AND R2

seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/FUR1/0.6/chip/SM001943-8_S14_L002_R1_001.fastq.gz 500000 > /home/omicas/luistorres/TFM/samples/FUR1/0.6/chip/subsampling/R1_chip_0.6.fastq.gz
seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/FUR1/0.6/chip/SM001943-8_S14_L002_R2_001.fastq.gz 500000 > /home/omicas/luistorres/TFM/samples/FUR1/0.6/chip/subsampling/R2_chip_0.6.fastq.gz 

## INPUT 0.6 --> R1 AND  R2

seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/FUR1/0.6/input/SM001943-2_S8_L002_R1_001.fastq.gz 500000 > /home/omicas/luistorres/TFM/samples/FUR1/0.6/input/subsampling/R1_input_0.6.fastq.gz
seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/FUR1/0.6/input/SM001943-2_S8_L002_R2_001.fastq.gz 500000 > /home/omicas/luistorres/TFM/samples/FUR1/0.6/input/subsampling/R2_input_0.6.fastq.gz

## AP 2.5 --> R1 AND R2

seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/FUR1/2.5/chip/SM001943-11_S17_L002_R1_001.fastq.gz 500000 > /home/omicas/luistorres/TFM/samples/FUR1/2.5/chip/subsampling/R1_chip_2.5.fastq.gz
seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/FUR1/2.5/chip/SM001943-11_S17_L002_R2_001.fastq.gz 500000 > /home/omicas/luistorres/TFM/samples/FUR1/2.5/chip/subsampling/R2_chip_2.5.fastq.gz

## INPUT 2.5 --> R1 ND R2

seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/FUR1/2.5/input/SM001943-5_S11_L002_R1_001.fastq.gz 500000 > /home/omicas/luistorres/TFM/samples/FUR1/2.5/input/subsampling/R1_input_2.5.fastq.gz 
seqtk sample -s 66 /home/omicas/luistorres/TFM/samples/FUR1/2.5/input/SM001943-5_S11_L002_R2_001.fastq.gz 500000 > /home/omicas/luistorres/TFM/samples/FUR1/2.5/input/subsampling/R2_input_2.5.fastq.gz 

### READ MAPPING 0.6M ###

## CHIP

cd /home/omicas/luistorres/TFM/samples/FUR1/0.6/chip/subsampling/

bowtie2 -x ../../../../../genome/index -U R1_chip_0.6.fastq.gz  -S chip_FUR1_0.6_sub.sam
samtools sort -o chip_FUR1_0.6_sub.bam chip_FUR1_0.6_sub.sam
samtools index chip_FUR1_0.6_sub.bam
rm chip_FUR1_0.6_sub.sam

## INPUT

cd /home/omicas/luistorres/TFM/samples/FUR1/0.6/input/subsampling

bowtie2 -x ../../../../../genome/index -U R1_input_0.6.fastq.gz  -S input_FUR1_0.6_sub.sam
samtools sort -o input_FUR1_0.6_sub.bam input_FUR1_0.6_sub.sam
samtools index input_FUR1_0.6_sub.bam
rm input_FUR1_0.6_sub.sam


### READ MAPPING 2.5M ###

## CHIP

cd /home/omicas/luistorres/TFM/samples/FUR1/2.5/chip/subsampling

bowtie2 -x ../../../../../genome/index -U R1_chip_2.5.fastq.gz -S chip_FUR1_2.5_sub.sam
samtools sort -o chip_FUR1_2.5_sub.bam chip_FUR1_2.5_sub.sam
samtools index chip_FUR1_2.5_sub.bam
rm chip_FUR1_2.5_sub.sam 

## INPUT

cd /home/omicas/luistorres/TFM/samples/FUR1/2.5/input/subsampling

bowtie2 -x ../../../../../genome/index -U R1_input_2.5.fastq.gz -S input_FUR1_2.5_sub.sam
samtools sort -o input_FUR1_2.5_sub.bam input_FUR1_2.5_sub.sam
samtools index input_FUR1_2.5_sub.bam
rm input_FUR1_2.5_sub.sam

### PEAK CALLING ###

## 0.6M

cd /home/omicas/luistorres/TFM/results/FUR1/0.6/

macs2 callpeak -t ../../../samples/FUR1/0.6/chip/subsampling/chip_FUR1_0.6_sub.bam -c ../../../samples/FUR1/0.6/input/subsampling/input_FUR1_0.6_sub.bam -f BAM --gsize 3.7e6 --outdir . -n FUR1_Peaks_0.6_sub

## 2.5

cd /home/omicas/luistorres/TFM/results/FUR1/2.5/

macs2 callpeak -t ../../../samples/FUR1/2.5/chip/subsampling/chip_FUR1_2.5_sub.bam -c ../../../samples/FUR1/2.5/input/subsampling/input_FUR1_2.5_sub.bam -f BAM --gsize 3.7e6 --outdir . -n FUR1_Peaks_2.5_sub
