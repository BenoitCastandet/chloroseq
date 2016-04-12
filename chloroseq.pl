#!/usr/bin/perl

=head1 NAME

chloroseq.pl

perl script that runs pipeline to:
1) get single nt coverage, window coverage, exons and introns rpkm for chloroplast genes, following tophat alignments of RNA-Seq data (Analysis 1)
2) get the splicing efficiency, following tophat alignement of RNA-Seq data (Analysis 2)
3) look at the editing sites in chloroplast, following tophat alignment of RNA-Seq data (Analysis 3)  

=cut

=head1 SYPNOSIS

chloroseq.pl [-h] -a <analysis> -b <accepted_hits.bam> -e <exon.gff3> -i <intron.gff3> -g <genome_size>
-n <name> -v <editing_sites.gff3> -f <fasta> -s <splice_sites.gff3> -k <keep_files> 

=head2 I<Flags:>

=over

=item -a

B<analysis>                    analysis to run; options: 1,2,3,all (mandatory)

=item -b

B<accepted_hits.bam>           aligned reads output from tophat in bam format (mandatory)

=item -e

B<exon.gff3>                   gff3 file of exon coordinates (mandatory)

=item -i

B<intron.gff3>                 gff3 file of intron coordinates (mandatory)

=item -g

B<genome_size>                 size of chloroplast genome in base pairs (mandatory)

=item -n

B<name>                        exact name of the chloroplast chromosome, ex: ChrC for Arabidopsis (mandatory)

=item -f 

B<fasta>                       genome file in fasta format (mandatory for analysis 3) 

=item -v 

B<editing_sites.gff3>	       file of known editing sites (optional for analysis 3)

=item -s

B<splice_sites.gff3>           gff3 file of splice sites (mandatory for analysis 2)

=item -k

B<keep_files>                  keep intermediate files

=item -h

B<help>                        print the help

=back

=cut

=head1 DESCRIPTION

This is a perl script to get single nt coverage, window coverage, exon rpkm, intron rpkm, and splicing and editing efficiencies for chloroplast genes following tophat alignments.

Bedtools v2.25.0 must be installed. If all 3 analyses are run it will produce nine files as listed:

    Analysis 1 files: counts, bam_name, exon_rpkm.txt, intron_rpkm.txt, nt_coverage.txt, and window_coverage.txt
    Analysis 2 file:  splicing_efficiency.txt
    Analysis 3 files: editing.pileup, editing_efficiency.txt 

=cut

=head1 AUTHORS

Benoit Castandet
(bc467@cornell.edu)

Susan Strickler
(srs57@cornell.edu)

=cut

=head1 METHODS

chloroseq.pl

=cut

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

#Get commandline options and check mandatory arguments are supplied
our ($opt_a, $opt_b, $opt_e, $opt_i, $opt_g, $opt_n, $opt_s, $opt_f, $opt_v, $opt_k, $opt_h);
getopts("a:b:e:i:g:c:n:s:f:v:kh");

if (!$opt_a && !$opt_b && !$opt_e && !$opt_i && !$opt_g && $opt_n && !$opt_s && !$opt_f && !$opt_v && !$opt_k && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

if ($opt_h) {
    help();
}

my $analysis = $opt_a
    || die("MANDATORY ARGUMENT ERROR: -a <analysis> was not supplied.\n");
my $bam = $opt_b
    || die("MANDATORY ARGUMENT ERROR: -b <accepted_hits.bam> was not supplied.\n");
my $exon = $opt_e
    || die("MANDATORY ARGUMENT ERROR: -e <exon.gff3> was not supplied.\n");
my $intron = $opt_i
    || die("MANDATORY ARGUMENT ERROR: -i <intron.gff3> was not supplied.\n");
my $name = $opt_n
    || die("MANDATORY ARGUMENT ERROR: -n <name> was not supplied.\n");
my $fasta  = $opt_f; 
my $genome = $opt_g;
my $snp    = $opt_v;
my $splice = $opt_s;
my $keep   = $opt_k;


#Check if snp file is provided, warn if ncessary
if ((!$opt_v && $opt_a =~ /all/) || (!$opt_v && $opt_a =~ /3/)){
    print "Warning: No SNP file provided, base count is not strand-specific in analysis 3.\n";
}
#Check for bedtools installation
my $bedtoolsv = `bedtools --version`;  

if ($bedtoolsv !~ "bedtools v2.25.0"){
    die("Version $bedtoolsv is wrong version.  Need v2.25.0. \n\n");
}

#run prep_bam.sh shell script to prep bam file
my $bam_name = $name . ".bam";

if (-e $bam_name) { 
    print "$bam_name file exists, moving to next step.\n"; 
}

else {
print STDERR "Starting prep_bam.sh.\n";
system("./prep_bam.sh $bam $name");
}

#Run Analysis 1 if chosen
if ($analysis =~ /1/ || $analysis =~ /all/){
    print STDERR "Starting coverage analysis (Analysis 1).\n";
    &coverage($name, $bam_name, $exon, $intron, $genome);
}

#Run Analysis 2 if chosen
if ($analysis =~ /2/ || $analysis =~ /all/){
    print STDERR "Starting splicing analysis (Analysis 2).\n";
    unless (defined $splice) {die "Need to supply splice file\n";}
    &splicing($splice, $exon, $intron);
}

#Run Analysis 3 if chosen
if ($analysis =~ /3/ || $analysis =~ /all/){
    print STDERR "Starting editing analysis (Analysis 3).\n";
    unless (defined $fasta) {die "Need to supply the fasta file\n";}
    &editing($fasta, $snp); 
}

##############################################################################################################
##########################################SUBROUTINES#########################################################
##############################################################################################################

#### ANALYSIS 1: coverage subroutine ####
sub coverage {
    #pass arguments
    my $name = $_[0];
    my $bam_name = $_[1];
    my $exon = $_[2];
    my $intron = $_[3];
    my $genome = $_[4];
    
    #generate sliding window files 100 nt sliding window every 50nt
    open (OUTPUT , ">100nt_plus.bed") || die "Could not open file '100nt_plus.bed\n";
    my $i = 0;
    
    while ($i < $genome - 100){ 
	my $end = $i + 100;
	print OUTPUT $name . "\t" . $i . "\t" . $end . "\tsliding_window\t\.\t+\n";
	$i += 50;
    }

    print OUTPUT $name . "\t" . $i . "\t" . $genome . "\tsliding_window\t\.\t+\n";
    close OUTPUT;
    
    open (OUTPUT , ">100nt_minus.bed") || die "Could not open file '100nt_minus.bed\n";

    $i =0;

    while ($i < $genome - 100){
	my $end = $i + 100;
	print OUTPUT $name . "\t" . $i . "\t" . $end . "\tsliding_window\t\.\t-\n";
	$i += 50;
    }
    
    print OUTPUT $name . "\t" . $i . "\t" . $genome . "\tsliding_window\t\.\t-\n";
    close OUTPUT;
    
    #run get_coverage.sh shell script
    print STDERR "Running get_coverage.sh (Analysis 1).\n";
    system("./get_coverage.sh $name $bam_name $exon $intron");
    
    ##join coverage files##
    #Open input file
    open (COVPLUS, "<coverage_plus.csv") || die "Cannot open plus coverage file.\n";
    open (COVMINUS, "<coverage_minus.csv") || die "Cannot open minus coverage file.\n";
    
    my %ids;
    while (<COVMINUS>) {
	chomp($_);
	my @feature = split (/\t/,$_);
	my $location = $feature[0] . "\t" . $feature[1];
	$ids{$location} = $feature[2];
    }
    
    open (OUTPUT , ">nt_coverage.txt") || die "Could not open file nt_coverage.txt";
    print OUTPUT "CHROMOSOME\tPOSITION\tPLUS\tMINUS\n";  
    
    # Select only if it is in the list
    while (<COVPLUS>){
	chomp ($_);
	my $line = $_;
	my @feature2 = split (/\t/,$line);
	my $location2 = $feature2[0] . "\t" . $feature2[1];
	
	if (exists $ids{$location2}) {
	    print OUTPUT $line . "\t" . $ids{$location2} . "\n";
	}
	
	else {
	    next;
	}
    }
     
     close OUTPUT;
    
     #Open input file
     open (WINPLUS, "<window_plus.csv") || die "Cannot open plus coverage file.\n";
     open (WINMINUS, "<window_minus.csv") || die "Cannot open minus coverage file.\n";
     open (ALLWIN, "<100nt_plus.bed") ||die "Cannot open 100nt_plus.bed file.\n";

     my %ids2;
     my %ids3;

     while (<WINPLUS>){
	 chomp($_);
	 my @feature = split (/\t/,$_);
	 my $location = $feature[0] . "\t" . $feature[1];
	 $ids2{$location} = $feature[3];
     }

     while (<WINMINUS>){
	 chomp($_);
	 my @feature = split (/\t/,$_);
	 my $location = $feature[0] . "\t" . $feature[1];
	 $ids3{$location} = $feature[3];
     }
     
     open (OUTPUT , ">window_coverage.txt") || die "Could not open file window_coverage.txt";
     print OUTPUT "CHROMOSOME\tSTART\tEND\tPLUS\tMINUS\n";
     
     ## Select only if it is in the list 
     while (<ALLWIN>) {
	chomp($_);
	my $line = $_;
	my @feature2 = split (/\t/,$line);
	my $location2 = $feature2[0] . "\t" . $feature2[1];      

	print OUTPUT $location2 . "\t" . $feature2[2] . "\t";
	
	if (exists $ids2{$location2}) {                                                                                                                               print OUTPUT $ids2{$location2} . "\t";                                                                                                 } 

	else {print OUTPUT "\t0\t";}

	if (exists $ids3{$location2}) {                                                                                                         
	    print OUTPUT $ids3{$location2} . "\n";                                                                                         
	}

	else {print OUTPUT "0\n";}
	
     }
    
     close OUTPUT;

     ##Join exon file with rpkm files
     #Open input file
     open (EXON, "<$exon") || die "Cannot open exon file.\n";
     open (ERPKM, "<exp_exon.txt") || die "Cannot open exp_exon.txt file.\n";
     open (COUNTS, "<counts") || die "Cannot open counts.\n";

     my $counts = <COUNTS>;
	   
     my %rpkm_ids;

     while (<ERPKM>) {
        chomp($_);
        my @feature = split (/\t/,$_);
        my $location = $feature[1] . "\t" . $feature[2];
        $rpkm_ids{$location} = $feature[0];
     }
     
     open (OUTPUT , ">exon_rpkm.txt") || die "Could not output file exon_rpkm.txt";
     print OUTPUT "CHROMOSOME\tSTART\tEND\tRPKM\n";
     
     ## Select only if it is in the list                                                                               
     while (<EXON>){
        chomp ($_);
        my $line = $_;
        my @feature = split (/\t/,$line);
        my $location = $feature[0] . "\t" . $feature[3];

        if (exists $rpkm_ids{$location}) {
	    my $rpkms = ($rpkm_ids{$location} * 1000000000)/($counts * ($feature[4] - $feature[3]));
	    print OUTPUT $location . "\t" . $feature[4] . "\t" . $rpkms . "\n";
        }

        else {
	    print OUTPUT $location . "\t" . $feature[4] . "\t0\n";
        }
    }

     close OUTPUT;

     ##Join intron file with rpkm files
     #Open input file
     open (INTRON, "<$intron") || die "Cannot open intron file.\n";
     open (EIRPKM, "<exp_intron.txt") || die "Cannot open exp_intron.txt file.\n";

     my %rpkm_ids2;

     while (<EIRPKM>) {
	 chomp($_);
	 my @feature = split (/\t/,$_);
	 my $location = $feature[1] . "\t" . $feature[2];
	 $rpkm_ids2{$location} = $feature[0];
     }

     open (OUTPUT , ">intron_rpkm.txt") || die "Could not output file intron_rpkm.txt";
     print OUTPUT "CHROMOSOME\tSTART\tEND\tRPKM\n";

     ## Select only if it is in the list
     while (<INTRON>){
	 chomp ($_);
	 my $line = $_;
	 my @feature = split (/\t/,$line);
	 my $location = $feature[0] . "\t" . $feature[3];

	 if (exists $rpkm_ids2{$location}) {
	     my $rpkms = ($rpkm_ids2{$location} * 1000000000)/($counts * ($feature[4] - $feature[3]));
	     print OUTPUT $location . "\t" . $feature[4] . "\t" . $rpkms . "\n";
	 }

	 else {
	     print OUTPUT $location . "\t" . $feature[4] . "\t0\n";
	 }
     }
     
     #remove the intermediate files if no -k options
     if (!$keep){
         system('rm coverage_plus.csv coverage_minus.csv 100nt_minus.bed 100nt_plus.bed window_plus.csv window_minus.csv exp_exon.txt exp_intron.txt');
     }
}

#### ANALYSIS 2: splicing subroutine ####
sub splicing {
    #pass arguments
    my $splice = $_[0];
    my $exon = $_[1];
    my $intron = $_[2];

    #run get_splicing_efficiency.sh shell script
    system("./get_splicing_efficiency.sh $bam_name $splice $exon $intron");

    ##join output files##
    #Open input file
    open (SPLICED, "<spliced_coverage.txt") || die "Cannot open spliced coverage file.\n";
    open (UNSPLICED, "<unspliced_coverage.txt") || die "Cannot open unspliced coverage file.\n";

    my %ids3;

    while (<SPLICED>) {
	chomp($_);
	my @feature = split (/\t/,$_);
	my $location = $feature[0] . "\t" . $feature[3] . "\t" . $feature[4];
	$ids3{$location} = $feature[5];
    }

    open (OUTPUT , ">splicing_efficiency.txt") || die "Could not open file window_coverage.txt";
    print OUTPUT "CHROMOSOME\tSTART\tSTOP\tSPLICED\tUNSPLICED\tEFFICIENCY\n";

    # Select only if it is in the list
    while (<UNSPLICED>){
	chomp ($_);
	my $line = $_;
	my @feature2 = split (/\t/,$line);
	my $location2 = $feature2[0] . "\t" . $feature2[3] . "\t" . $feature2[4];

	if (exists $ids3{$location2}) {
	    my $speff =  ($ids3{$location2} / ( $ids3{$location2} + ($feature2[5]/2))) * 100;
	    print OUTPUT $location2 . "\t" . $ids3{$location2} . "\t" . $feature2[5] . "\t" .  $speff . "\n";
	}

	else {
	    next;
	}
    }

    close OUTPUT;

    #remove the intermediate files if no -k options
    if (!$keep){
	system('rm splice_junctions.bam spliced_coverage.txt unspliced_coverage.txt');
    }
}
    
#### ANALYSIS 3: editing subroutine ####
sub editing {
    #pass arguments
    my $fasta = $_[0];
    my $snp = $_[1];

    if (-e "editing.pileup") {
	print "Editing pileup file exists, moving to next step.\n";
    }

    else{
	#run get_editing_efficiency.sh shell script
	system("./get_editing_efficiency.sh $bam_name $fasta");
    }
    
    ##extract snp frequencies##
    #Open input file
    open (PILEUP, "<editing.pileup") || die "Cannot open pileup file.\n";

    my %ids;
    my $feature;
    my $edit = 0;

    #store sites of interest in hash
    if ($snp){
	open (GFF, "<$snp");
	while (<GFF>) {
	    chomp($_);
	    my @sites = split(/\t/,$_);
	    my $feature = $sites[0] . "_" . $sites[3];
	    $ids{$feature} =1;
	}
    }
	
    #find sites of interest (if required) and calculate base frequency at each site of pileup file
    open (OUTPUT , ">editing.txt") || die "Could not open file 'editing.txt\n";
    if (defined $snp){
	print OUTPUT "Chromosome\tPosition\tReference\tTotal\tA\tT\tG\tC\tEfficiency\n";
    }

    else {print OUTPUT "Chromosome\tPosition\tReference\tTotal\tA\tT\tG\tC\n";}
    
    while (<PILEUP>){
	chomp ($_);
	my $line = $_;
	my @cells= split (/\t/,$line);

	#make sure all 5 elements of array are initialized (element 4 and 5 won't be assigned if there is a gap)
	if (scalar @cells == 4){
	    push @cells, "0", "0";
	}
	
	#get rid of indels
	if ($cells[4] =~ /(-1.)/){
	    $cells[4] =~ s/-1.//g;
	}
	
	if ($cells[4] =~ s/-2..//g){
	    $cells[4] =~ s/-2..//g;
	}
	
	if ($cells[4] =~ s/\+1.//g){
	    $cells[4] =~ s/\+1.//g;
	}
	
	if ($cells[4] =~ s/\+2..//g){
	    $cells[4] =~ s/\+2..//g;
	}
	
	#make sure all indels are removed
	if ($cells[4] =~ /\+/ || $cells[4] =~ /-/){
	    die "Indel size not accounted for\n";
	}
	
	#count bases
	my $name = $cells[0] . "_" . $cells[1];
	
	if ($cells[4] =~ /[AaTtGgCc\.\,]/){
	    my $ca_count = $cells[4] =~ tr/[A]//;
	    my $a_count = $cells[4] =~ tr/[a]//;
	    my $ct_count = $cells[4] =~ tr/[T]//;
	    my $t_count = $cells[4] =~ tr/[t]//;
	    my $cg_count = $cells[4] =~ tr/[G]//;
	    my $g_count = $cells[4] =~ tr/[g]//;
	    my $cc_count = $cells[4] =~ tr/[C]//;
	    my $c_count = $cells[4] =~ tr/[c]//;
	    my $ref1 = $cells[4] =~ tr/\.//;
	    my $ref2 = $cells[4] =~ tr/\,//;
	    
	    #Only choose selected if a snp gff3 file is provided
	    if (defined $snp){

		#Add strand-specific snps
		if ($cells[2] =~ /C/){
		    $c_count = $ref1;
		}

		if ($cells[2] =~ /G/){
		    $g_count = $ref2;
		}
		
		if (exists $ids{$name}) {
		    if ($cells[2] =~ /G/){
			if ($a_count != 0 || $ref2 != 0){
			    $edit = ($a_count / ($a_count + $ref2)) * 100;
			}
			else {$edit = 0;}
		    }
		    
		    elsif ($cells[2] =~ /C/){
			if ($ct_count != 0 || $ref1 != 0){
			    $edit = ($ct_count / ($ct_count + $ref1)) * 100;
			}
			
			else {$edit = 0;}
		    }
		    
		    else {$edit = "NA";}

		    if ($cells[2] =~ /C/){
			print OUTPUT $cells[0] . "\t" . $cells[1] . "\t" . $cells[2] . "\t" . $cells[3] . "\t" . $ca_count . "\t" . $ct_count . "\t". $cg_count . "\t" . $ref1 . "\t" . $edit . "\n";
		    }

		    elsif ($cells[2] =~ /G/){
			print OUTPUT $cells[0] . "\t" . $cells[1] . "\t" . $cells[2] . "\t" . $cells[3] . "\t" . $a_count . "\t" . $t_count . "\t". $ref2 . "\t" . $c_count . "\t" . $edit . "\n";
		    }
		
		    else {
			next;
		    }
		}
	    }
	    
	    #In no snp gff3 file provided, print all sites
	    else {
		my $a_count2 = 0;
		my $t_count2 = 0;
		my $g_count2 = 0;
		my $c_count2 = 0; 

		#Counting bases
		if ($cells[2] =~ /[Aa]/){
		    $a_count2 = $ref1 + $ref2;
		}

		elsif ($cells[2] =~ /[Tt]/){
		    $t_count2 = $ref1 + $ref2;
		}

		elsif ($cells[2] =~ /[Gg]/){
		    $g_count2 = $ref1 + $ref2;
		}

		elsif ($cells[2] =~ /[Cc]/){
		    $c_count2 = $ref1 + $ref2;
		}

		#Get base frequency
		my $a_freq = $ca_count + $a_count + $a_count2;
		my $t_freq = $ct_count + $t_count + $t_count2;
		my $g_freq = $cg_count + $g_count + $g_count2;
		my $c_freq = $cc_count + $c_count + $c_count2;                            
		
		print OUTPUT $cells[0] . "\t" . $cells[1] . "\t" . $cells[2] . "\t" . $cells[3] . "\t" . $a_freq . "\t" . $t_freq . "\t". $g_freq . "\t" . $c_freq . "\n";
	    }
		
	}
	
	else {
	    next;
	}
    }
    
    close OUTPUT;
}

=head2 help

  Usage: help()
  Desc: print help of this script
  Ret: none
  Args: none
  Side_Effects: exit of the script
  Example: if (!@ARGV) {
               help();
           }

=cut

sub help {
    print STDERR <<EOF;
    $0:
    
  Description:

    perl script that runs pipeline to:                                                                         
    - Analysis 1: get single nt coverage, window coverage, exons and introns rpkm for chloroplast genes, following tophat alignment of RNA-Seq data
    - Analysis 2: get the splicing efficiency, following tophat alignment of RNA-Seq data
    - Analysis 3: get the editing efficiency following tophat alignment of RNA-Seq data
	
	For example:     chloroseq.pl -a all -b accepted_hits.bam -e TAIR10_ChrC_files/TAIR10_ChrC_exon.gff3 
                         -i TAIR10_ChrC_files/TAIR10_ChrC_introns.gff3 -g 158457 -n ChrC -f TAIR10_ChrC_files/TAIR10_ChrC_bowtie2_index/TAIR10_ChrC.fa
			 -v TAIR10_ChrC_files/TAIR10_ChrC_editing_sites.gff3 -s TAIR10_ChrC_files/TAIR10_ChrC_splice_sites_sort.gff3
	
      Usage:
	
          chloroseq.pl [-h] -a <analysis> -b <accepted_hits.bam> -e <exon.gff3> -i <intron.gff3> -g <genome_size>
          -n <name> -v <editing_sites.gff3> -f <fasta> -v <snps> -s <splice_sites.gff3 -k <keep_files> > output.tab   
	
      Flags:
	
	-a <analysis>                    analysis to run; options: 1,2,3,all (mandatory)
        -b <accepted_hits>               aligned reads output from tophat in bam format (mandatory)
        -e <exons>                       exon coordinates in gff3 format (mandatory)
        -i <intron>                      intron coordinates in gff3 format (mandatory)
        -g <genome_size>                 size of chloroplast genome in base pairs (mandatory)
        -n <name>                        exact name of the chloroplast, ex: ChrC for Arabidopsis (mandatory)
        -f <fasta>                       genome files in fasta format (mandatory for Analysis 2)
        -v <snps>                        file of known editing sites (optional for Analysis 2)
        -s <splice_sites>                splice sites in gff3 format (mandatory for Analysis 3)
        -k <keep_files>                  keep intermediate files (optional)
	-h <help>                        print the help
	
EOF
exit (1);
}
