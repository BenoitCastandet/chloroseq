#!/bin/bash

#Shell script called in chloroseq.pl that prepares bam file for further processing

#Collect variables from script
bam=$1
name=$2

#index bam file
samtools index $bam

#get chloroplast reads
samtools view -b $bam "$name" > $name.bam
