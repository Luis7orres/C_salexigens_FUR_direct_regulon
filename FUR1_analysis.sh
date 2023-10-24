#!/bin/bash
#SBATCH --job-name=FUR1

####  FUR1 ANALYSIS  ####

### READ MAPPING  0.6M ###

## CHIP

cd /home/omicas/luistorres/TFM/samples/FUR1/0.6/chip/

bowtie2 -x ../../../../genome/index -1 SM001943-8_S14_L002_R1_001.fastq.gz -2 SM001943-8_S14_L002_R2_001.fastq.gz -S chip_FUR1_0.6.sam
samtools sort -o chip_FUR1_0.6.bam chip_FUR1_0.6.sam
samtools index chip_FUR1_0.6.bam
rm chip_FUR1_0.6.sam

## INPUT

cd /home/omicas/luistorres/TFM/samples/FUR1/0.6/input/

bowtie2 -x ../../../../genome/index -1 SM001943-2_S8_L002_R1_001.fastq.gz -2 SM001943-2_S8_L002_R2_001.fastq.gz -S input_FUR1_0.6.sam
samtools sort -o input_FUR1_0.6.bam input_FUR1_0.6.sam
samtools index input_FUR1_0.6.bam
rm input_FUR1_0.6.sam


### READ MAPPING 2.5M ###

## CHIP

cd /home/omicas/luistorres/TFM/samples/FUR1/2.5/chip/

bowtie2 -x ../../../../genome/index -1 SM001943-11_S17_L002_R1_001.fastq.gz -2 SM001943-11_S17_L002_R2_001.fastq.gz -S chip_FUR1_2.5.sam
samtools sort -o chip_FUR1_2.5.bam chip_FUR1_2.5.sam
samtools index chip_FUR1_2.5.bam
rm chip_FUR1_2.5.sam 

## INPUT

cd /home/omicas/luistorres/TFM/samples/FUR1/2.5/input/

bowtie2 -x ../../../../genome/index -1 SM001943-5_S11_L002_R1_001.fastq.gz -2 SM001943-5_S11_L002_R2_001.fastq.gz -S input_FUR1_2.5.sam
samtools sort -o input_FUR1_2.5.bam input_FUR1_2.5.sam
samtools index input_FUR1_2.5.bam
rm input_FUR1_2.5.sam

### PEAK CALLING ###

## FUR1 0.6M NaCl

cd /home/omicas/luistorres/TFM/results/FUR1/0.6/

macs2 callpeak -t ../../../samples/FUR1/0.6/chip/chip_FUR1_0.6.bam -c ../../../samples/FUR1/0.6/input/input_FUR1_0.6.bam -f BAM --outdir . --gsize 3.7e6 -n FUR1_Peaks_0.6

## FUR1 2.5M NaCl

cd /home/omicas/luistorres/TFM/results/FUR1/2.5/

macs2 callpeak -t ../../../samples/FUR1/2.5/chip/chip_FUR1_2.5.bam -c ../../../samples/FUR1/2.5/input/input_FUR1_2.5.bam -f BAM --outdir . --gsize 3.7e6 -n FUR1_Peaks_2.5
