# chloroseq
a bioinformatic pipeline to systematically analyse the chloroplast transcriptome using RNA-Seq

This is a perl script that runs a bioinformatic pipeline to:
1) get single nt coverage, window coverage, exons and introns rpkm for chloroplast genes, following tophat alignments of RNA-Seq data (Analysis 1)
2) get the splicing efficiency, following tophat alignement of RNA-Seq data (Analysis 2)
3) look at the editing sites in chloroplast, following tophat alignment of RNA-Seq data (Analysis 3) 

Bedtools v2.25.0 must be installed. If all 3 analyses are run it will produce nine files as listed:

    Analysis 1 files: counts, bam_name, exon_rpkm.txt, intron_rpkm.txt, nt_coverage.txt, and window_coverage.txt
    Analysis 2 file:  splicing_efficiency.txt
    Analysis 3 files: editing.pileup, editing_efficiency.txt 
    
The ChloroSeq_for_bin directory contains the scripts to download if willing to put them in a bin directory.
The ChloroSeq_scripts directory contains the scripts to download if willing to use them in the working directory (./ to call the scripts).

AUTHORS

Benoît Castandet
(benoit.castandet@ips2.universite-paris-saclay.fr)

Susan Strickler
(srs57@cornell.edu)
