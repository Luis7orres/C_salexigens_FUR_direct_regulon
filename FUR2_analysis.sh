#!/bin/bash
#SBATCH --job-name=2.5

#### FUR2 ANALYSIS  ####

### READ MAPPING 0.6M ###

## CHIP

cd /home/omicas/luistorres/TFM/samples/AP/chip/0.6/chip_todas
bowtie2 -x ../../../../../genome/index -1 ../SM001943-9_S15_L002_R1_001.fastq.gz -2 ../SM001943-9_S15_L002_R2_001.fastq.gz -S chip_0.6_AP.sam
samtools sort -o chip_0.6_AP.bam chip_0.6_AP.sam
samtools index chip_0.6_AP.bam
rm chip_0.6_AP.sam

## INPUT

cd /home/omicas/luistorres/TFM/samples/AP/input/0.6/input_todas

bowtie2 -x ../../../../../genome/index -1 ../SM001943-3_S9_L002_R1_001.fastq.gz -2 ../SM001943-3_S9_L002_R2_001.fastq.gz -S input_0.6_AP.sam
samtools sort -o input_0.6_AP.bam input_0.6AP.sam
samtools index input_0.6_AP.bam
rm input_0.6_AP.sam

### READ MAPPING 2.5M ###

## CHIP

cd /home/omicas/luistorres/TFM/samples/AP/chip/2.5/chip_todas
bowtie2 -x ../../../../../genome/index -1 ../SM001943-12_S18_L002_R1_001.fastq.gz -2 ../SM001943-12_S18_L002_R2_001.fastq.gz -S chip_2.5_AP.sam
samtools sort -o chip_2.5_AP.bam chip_2.5_AP.sam
samtools index chip_2.5_AP.bam
rm chip_2.5_AP.sam

## INPUT

cd /home/omicas/luistorres/TFM/samples/AP/input/2.5/input_todas

bowtie2 -x ../../../../../genome/index -1 ../SM001943-6_S12_L002_R1_001.fastq.gz -2 ../SM001943-6_S12_L002_R2_001.fastq.gz -S input_2.5_AP.sam
samtools sort -o input_2.5_AP.bam input_2.5_AP.sam
samtools index input_2.5_AP.bam
rm input_2.5_AP.sam


### PEAK CALLING ###

## FUR2 0.6M NaCl

cd /home/omicas/luistorres/TFM/results/0.6/AP/all
macs2 callpeak -t ../../../samples/AP/chip/0.6/chip_todas/chip_0-6_AP.bam -c ../../../samples/AP/input/0.6/input_todas/input_0.6_AP.bam -f BAM --outdir . --gsize 3.69e7 --nomodel --extsize 147 -n AP_peaks_2.5

## FUR2 2.5M NaCl

cd /home/omicas/luistorres/TFM/results/2.5/AP/all
macs2 callpeak -t ../../../samples/AP/chip/2.5/chip_todas/chip_2.5_AP.bam -c ../../../samples/AP/input/2.5/input_todas/input_2.5_AP.bam -f BAM --outdir . --gsize 3.69e7 --nomodel --extsize 147 -n AP_peaks_2.5


