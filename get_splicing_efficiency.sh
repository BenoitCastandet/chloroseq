#!/bin/bash

#Collect variables from script                
bam_name=$1;
splice=$2;
exon=$3;
intron=$4;

# strategy: run bedtools intersect to extract all the reads intersecting the exon/intron junction, in 3 steps: 
#1 find all the reads that intersect the splice sites, use -s to force strandeness
#2 eliminate purely exonic reads 
#3 eliminate purely intronic reads
                                            
bedtools intersect -s -a $bam_name -b $splice | bedtools intersect -a stdin -b $intron | bedtools intersect -a stdin -b $exon > splice_junctions.bam

#use samtools to only extract the reads that are spliced (exon/exon) from the splice_junctions.bam file and calculate the coverage using bedtools coverage
samtools view -h splice_junctions.bam | awk '$6 ~ /N/ || $1 ~ /^@/' | samtools view -bS - | bedtools coverage -bed -a $intron -b stdin |cut -f 1,2,3,4,5,10 > spliced_coverage.txt

#use samtools to only extract the reads that are NOT spliced (exon/intron) from the splice_junctions.bam file and calculate the coverage using bedtools coverage
samtools view -h splice_junctions.bam | awk '$6 !~ /N/ || $1 ~ /^@/' | samtools view -bS - | bedtools coverage -bed -a $intron -b stdin |cut -f 1,2,3,4,5,10 > unspliced_coverage.txt

