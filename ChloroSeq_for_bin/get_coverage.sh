#!/bin/bash

#Collect variables from script
name=$1
bam_name=$2
exon=$3
intron=$4

#get counts
counts="$(samtools view -F 256 $bam_name | wc -l)" 
echo "$counts" > counts

echo "$bam_name" > bam_name

#get strand specific, single base coverage using bedtools on the bam file from the tophat alignment. use -s for strandeness. First for plus strand, then minus strand. Normalized with the number of aligned reads.
bedtools genomecov -strand + -d -ibam $bam_name -g genome_ChrC.bed | grep "^$name" | awk -v awk_count=$counts '{print $1 "\t" $2 "\t" ($3*1000000)/(awk_count)}' > coverage_plus.csv
bedtools genomecov -strand - -d -ibam $bam_name -g genome_ChrC.bed | grep "^$name" | awk -v awk_count=$counts '{print $1 "\t" $2 "\t" ($3*1000000)/(awk_count)}'> coverage_minus.csv

#get strand specific coverage on a 100nt sliding window (every 50nt), normalized by the number of aligned reads
bedtools intersect -s -bed -wb -a $bam_name -b 100nt_plus.bed | cut -f 13,14,15 | sort -k2,2n | uniq -c |  awk -v awk_count=$counts '{print $2 "\t" $3 "\t" $4 "\t" ($1*1000000)/(awk_count)}' > window_plus.csv
bedtools intersect -s -bed -wb -a $bam_name -b 100nt_minus.bed | cut -f 13,14,15 | sort -k2,2n | uniq -c | awk -v awk_count=$counts '{print $2 "\t" $3 "\t" $4 "\t" ($1*1000000)/(awk_count)}' > window_minus.csv

#calculate RPKM for all the exons
bedtools intersect -s -bed -wb -a $bam_name -b $exon | cut -f 13,16,17 | sort -k2,2n | uniq -c |sed -r 's/^( *[^ ]+) +/\1\t/' > exp_exon.txt

#calculate RPKM for all the introns
bedtools intersect -s -bed -wb -a $bam_name -b $intron | cut -f 13,16,17 | sort -k2,2n | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' > exp_intron.txt
