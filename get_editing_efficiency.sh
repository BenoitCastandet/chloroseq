#!/bin/bash

#Collect variables from script
bam_name=$1
fasta=$2

# create a pileup file containing the SNP for the genome of interest using samtools mpileup, -I to not consider indels
samtools mpileup -I -f $fasta $bam_name > editing.pileup 
