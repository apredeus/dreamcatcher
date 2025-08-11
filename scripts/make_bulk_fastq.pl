#!/usr/bin/env perl 

## prepare unmapped reads for bulk RNA-seq. If BAMs are used, any mapper would do. 
## only STARsolo can be used for the unmapped reads 
## will apply similar homopolymer filtering as make_bulk_fastq.pl; also, all reads with at least 1 mate shorter than 40 bp are discarded

use strict; 
use warnings;

if (scalar @ARGV != 1) { 
    print STDERR "USAGE: ./make_bulk_fastq.pl <output directory of bulk RNA-seq mapped with STAR or other tool>\n"; 
    exit 1; 
} 

my $dir = shift @ARGV;
my $out = shift @ARGV;
my $cmd = $ENV{'CMD'};
my $len_cutoff = 40; 

if (-s $dir."/Log.final.out") { 
    print STDOUT "make_bulk_fastq.pl: sample is provided in the form of STAR output! Will attempt to find Unmapped.out.mate* files..\n";
   
    ## we store these reads archived so let's account for this  
    my $paired = 0;
    my $R1 = `find $dir/* | grep Unmapped.out.mate1`; 
    my $R2 = `find $dir/* | grep Unmapped.out.mate2`;
    chomp $R1; 
    chomp $R2; 
    if ($R1 eq "") { 
        print STDERR "make_bulk_fastq.pl: WARNING: cannot find unmapped read files Unmapped.out.mate(1/2) in the STAR directory! Will try the BAM file ..\n";
        goto STARBAM; 
    } elsif ($R2 eq "") { 
        print STDOUT "make_bulk_fastq.pl: found unmapped reads in fastq format; sample was determined to be a bulk single-end experiment!\n";
    } else { 
        print STDOUT "make_bulk_fastq.pl: found unmapped reads in fastq format; sample was determined to be a bulk paired-end experiment!\n";
        $paired = 1; 
    }

    if ($paired) { 
        print STDOUT "make_bulk_fastq.pl: using the following input files:\n\tR1 = $R1\n\tR2 = $R2\n";     
        my $star_log = $dir."/Log.final.out"; 
        my $total_reads = `grep "Number of input reads" $star_log | awk '{print \$6}' | tr -d '\n'`; 
        my $total_unmapped = `grep "Number of reads unmapped" $star_log | awk -F "|" '{print \$2}' | awk '{sum+=\$1} END {print sum}' | tr -d '\n'`; 
        print STDOUT "make_bulk_fastq.pl: using STARsolo input: $total_reads total reads, $total_unmapped reads not mapped to host genome..\n"; 

        open OUT1,">","Unmapped_filt.R1.fastq" or die "$!"; 
        open OUT2,">","Unmapped_filt.R2.fastq" or die "$!"; 
        ## this occasionally stupid concoction reads 2 fastq files synchronously.
        if ($R1 =~ m/gz$/ && $R2 =~ m/gz$/) { 
            open READ1,"$cmd pigz -cd $R1 |" or die "$!"; 
            open READ2,"$cmd pigz -cd $R2 |" or die "$!"; 
        } elsif ($R1 =~ m/bz2$/ && $R2 =~ m/bz2$/) { 
            open READ1,"$cmd pbzip2 -cd $R1 |" or die "$!"; 
            open READ2,"$cmd pbzip2 -cd $R2 |" or die "$!";
        } else { 
            open READ1,"$cmd cat $R1 |" or die "$!"; 
            open READ2,"$cmd cat $R2 |" or die "$!";
        }

        my $nr = 1; 
        my ($rname1,$seq1,$len1,$qual1,$rname2,$seq2,$len2,$qual2) = ('') x 8;
        my ($passed,$discarded) = (0) x 2; 
        while (<READ1>) {
            my $r1 = $_; 
            chomp $r1;
            my $r2 = <READ2>;
            chomp $r2;
            if ($nr % 4 == 1) {
                $rname1 = (split /\s/,$r1)[0];
                $rname2 = (split /\s/,$r2)[0];
                print STDERR "make_bulk_fastq.pl: WARNING - read names do not match: $rname1 =/= $rname2!\n" if ($rname1 ne $rname2);
            } elsif ($nr % 4 == 2) { 
                $len1 = length $r1; 
                $len2 = length $r2; 
                $seq1 = $r1;
                $seq2 = $r2;
            } elsif ($nr % 4 == 0) { 
                $qual1 = $r1;
                $qual2 = $r2;
                if (! low_complexity($seq1) && ! low_complexity($seq2) && $len1 >= $len_cutoff && $len2 >= $len_cutoff) { 
                    printf OUT1 "%s\n%s\n+\n%s\n",$rname1,$seq1,$qual1;
                    printf OUT2 "%s\n%s\n+\n%s\n",$rname2,$seq2,$qual2;
                    $passed++; 
                } else {
                    $discarded++; 
                }
            } 
            $nr++; 
        }
        close READ1;
        close READ2;
        close OUT1; 
        close OUT2;
        print STDOUT "make_bulk_fastq.pl: writing filtered unmapped reads: outputted $passed reads to Unmapped_filt.R1/2.fastq, discarded $discarded reads as low-complexity/short (<40bp) sequences..\n";
    } else { 
        ## same but for single-end reads 
        print STDOUT "make_bulk_fastq.pl: using the following input files:\n\tR1 = $R1\n";     
        my $star_log = $dir."/Log.final.out"; 
        my $total_reads = `grep "Number of input reads" $star_log | awk '{print \$6}' | tr -d '\n'`; 
        my $total_unmapped = `grep "Number of reads unmapped" $star_log | awk -F "|" '{print \$2}' | awk '{sum+=\$1} END {print sum}' | tr -d '\n'`; 
        print STDOUT "make_bulk_fastq.pl: using STARsolo input: $total_reads total reads, $total_unmapped reads not mapped to host genome..\n"; 
        
        open OUT1,">","Unmapped_filt.R1.fastq" or die "$!"; 
        if ($R1 =~ m/gz$/) { 
            open READ1,"$cmd pigz -cd $R1 |" or die "$!"; 
        } elsif ($R1 =~ m/bz2$/) { 
            open READ1,"$cmd pbzip2 -cd $R1 |" or die "$!"; 
        } else { 
            open READ1,"$cmd cat $R1 |" or die "$!"; 
        }

        my $nr = 1; 
        my ($rname1,$seq1,$len1,$qual1) = ('') x 4;
        my ($passed,$discarded) = (0) x 2; 
        while (<READ1>) {
            my $r1 = $_; 
            chomp $r1;
            if ($nr % 4 == 1) {
                $rname1 = (split /\s/,$r1)[0];
            } elsif ($nr % 4 == 2) { 
                $len1 = length $r1; 
                $seq1 = $r1;
            } elsif ($nr % 4 == 0) { 
                $qual1 = $r1;
                if (! low_complexity($seq1) && $len1 >= $len_cutoff) { 
                    printf OUT1 "%s\n%s\n+\n%s\n",$rname1,$seq1,$qual1;
                    $passed++; 
                } else {
                    $discarded++; 
                }
            } 
            $nr++; 
        }
        close READ1;
        close OUT1; 
        print STDOUT "make_bulk_fastq.pl: writing filtered unmapped reads: outputted $passed reads to Unmapped_filt.R1.fastq, discarded $discarded reads as low-complexity/short (<40bp) sequences..\n";
    } 
} else {
    ## any type of BAM file that has unmapped reads in it will do. 
    ## we will take only unmapped reads (-f4) for single-end, and separate R1/R2 for paired-end (-f76=4+8+64/-f140=4+8+128)
    ## change the heuristics to find the right bam - this is just a guess based on my STAR output dirs..
    STARBAM:
    my $bam = `find $dir/* | grep -iv transcriptome | grep \"\.bam\$\"`; 
    chomp $bam; 
    if (length $bam) { 
        print STDOUT "make_bulk_fastq.pl: sample is assumed to be bulk RNA-seq BAM with unmapped reads!..\n";
    } else { 
        print STDERR "make_bulk_fastq.pl: ERROR - failed to find a BAM file in the directory provided! Exiting .. \n"; 
        exit 1; 
    }
    print STDOUT "make_bulk_fastq.pl: using the following input files: \n\tBAM = $bam\n";

    ## we need to check if the file was mapped as paired-end; for bulk, this changes some stats and fastq extraction samtools command.
    my $preads = `grep \"properly paired\" host.bam.flagstat | awk '{printf "%d",\$1+\$3}'`; 
    my $paired = ($preads > 0) ? 1 : 0; 
    if ($paired) { 
        my $total_reads  = `grep "primary\$"      host.bam.flagstat | awk '{printf "%d",\$1/2+\$3/2}'`; 
        my $total_mapped = `grep "primary mapped" host.bam.flagstat | awk '{printf "%d",\$1/2+\$3/2}'`;
        my $total_unmapped = $total_reads - $total_mapped; 
        print STDOUT "make_bulk_fastq.pl: using BAM input: $total_reads total reads, $total_unmapped reads not mapped to host genome..\n"; 
        print STDOUT "make_bulk_fastq.pl: extracting unmapped paired-end reads from the BAM file!\n";
        system "$cmd samtools fastq -\@4 -f4 -1 Unmapped_unfilt.R1.fastq -2 Unmapped_unfilt.R2.fastq -s Unmapped_unfilt.S.fastq $bam"; 
    } else { 
        my $total_reads  = `grep "primary\$"      host.bam.flagstat | awk '{printf "%d",\$1+\$3}'`; 
        my $total_mapped = `grep "primary mapped" host.bam.flagstat | awk '{printf "%d",\$1+\$3}'`;
        my $total_unmapped = $total_reads - $total_mapped; 
        print STDOUT "make_umi_fastq.pl: using BAM input: $total_reads total reads, $total_unmapped reads not mapped to host genome..\n"; 
        print STDOUT "make_bulk_fastq.pl: extracting unmapped single-end reads from the BAM file!\n";
        system "$cmd samtools fastq -\@4 -f4 -0 Unmapped_unfilt.R1.fastq $bam";
    } 
    
    ## now we do the same thing as with the unmapped reads above (sync reading) 
    if ($paired) { 
        open OUT1,">","Unmapped_filt.R1.fastq" or die "$!"; 
        open OUT2,">","Unmapped_filt.R2.fastq" or die "$!"; 
        open READ1,"$cmd cat Unmapped_unfilt.R1.fastq |" or die "$!";
        open READ2,"$cmd cat Unmapped_unfilt.R2.fastq |" or die "$!";

        my $nr = 1; 
        my ($rname1,$seq1,$len1,$qual1,$rname2,$seq2,$len2,$qual2) = ('') x 8;
        my ($passed,$discarded) = (0) x 2; 
        while (<READ1>) {
            my $r1 = $_; 
            chomp $r1;
            my $r2 = <READ2>;
            chomp $r2;
            if ($nr % 4 == 1) {
                $rname1 = (split /\s/,$r1)[0];
                $rname2 = (split /\s/,$r2)[0];
                print STDERR "make_bulk_fastq.pl: WARNING - read names do not match: $rname1 =/= $rname2!\n" if ($rname1 ne $rname2);
            } elsif ($nr % 4 == 2) { 
                $len1 = length $r1; 
                $len2 = length $r2; 
                $seq1 = $r1;
                $seq2 = $r2;
            } elsif ($nr % 4 == 0) { 
                $qual1 = $r1;
                $qual2 = $r2;
                if (! low_complexity($seq1) && ! low_complexity($seq2) && $len1 >= $len_cutoff && $len2 >= $len_cutoff) { 
                    printf OUT1 "%s\n%s\n+\n%s\n",$rname1,$seq1,$qual1;
                    printf OUT2 "%s\n%s\n+\n%s\n",$rname2,$seq2,$qual2;
                    $passed++; 
                } else {
                    $discarded++; 
                }
            } 
            $nr++; 
        }
        close READ1;
        close READ2;
        close OUT1; 
        close OUT2;
        print STDOUT "make_bulk_fastq.pl: writing filtered unmapped reads: outputted $passed reads to Unmapped_filt.R1/2.fastq, discarded $discarded reads as low-complexity/short (<40bp) sequences..\n";
    } else { 
        ## same but for single-end reads 
        open OUT1,">","Unmapped_filt.R1.fastq" or die "$!"; 
        open READ1,"$cmd cat Unmapped_unfilt.R1.fastq |" or die "$!"; 

        my $nr = 1; 
        my ($rname1,$seq1,$len1,$qual1) = ('') x 4;
        my ($passed,$discarded) = (0) x 2; 
        while (<READ1>) {
            my $r1 = $_; 
            chomp $r1;
            if ($nr % 4 == 1) {
                $rname1 = (split /\s/,$r1)[0];
            } elsif ($nr % 4 == 2) { 
                $len1 = length $r1; 
                $seq1 = $r1;
            } elsif ($nr % 4 == 0) { 
                $qual1 = $r1;
                if (! low_complexity($seq1) && $len1 >= $len_cutoff) { 
                    printf OUT1 "%s\n%s\n+\n%s\n",$rname1,$seq1,$qual1;
                    $passed++; 
                } else {
                    $discarded++; 
                }
            } 
            $nr++; 
        }
        close READ1;
        close OUT1; 
        print STDOUT "make_bulk_fastq.pl: writing filtered unmapped reads: outputted $passed reads to Unmapped_filt.R1.fastq, discarded $discarded reads as low-complexity/short (<40bp) sequences..\n";
    }
}

sub low_complexity {
  my $low_comp = 0; 
  ## identify low complexity sequences of few kinds:
  ## (1) over 80% (changed from 90% because of empiric evidence!) of the same letter overall;
  ## (2) a run spanning 60%+ of the read of the most frequent letter with 1 mismatch allowed
  my $str = shift @_;
  chomp $str;
  my $C = {};
  ## find the most common letter in the string
  my $top_ltr;
  my $top_count = 0;
  my @ltrs = split //,$str;
  foreach my $ltr (@ltrs) {
    if (! defined $C->{$ltr}) {
      $C->{$ltr} = 1;
    } else {
      $C->{$ltr}++;
    }
    if (! defined $top_ltr) {
      $top_ltr = $ltr;
      $top_count = 1;
    } elsif ($C->{$ltr} >= $top_count) {
      $top_ltr = $ltr;
      $top_count = $C->{$ltr};
    }
  }

  ## now calculate fraction of the top letter & set low_comp to 1 if it's 80%+ 
  my $str_len = length $str;
  $low_comp = 1 if ($top_count/$str_len >= 0.8);

  ## find longest substring of $top_ltr, allowing for 1 mismatch
  my $mismatch = 0;
  my $max_substr = 0;
  my $cur_substr = 0;
  foreach my $ltr (@ltrs) {
    if ($ltr eq $top_ltr) {
      $cur_substr++;
      $max_substr = ($cur_substr > $max_substr) ? $cur_substr : $max_substr;
    } elsif ($ltr ne $top_ltr && $cur_substr == 0) {
      next; ## don't count mismatches if we are not in the homopolymer yet
    } elsif ($ltr ne $top_ltr && $mismatch == 0) {
      $mismatch = 1;
      $cur_substr++;
      $max_substr = ($cur_substr > $max_substr) ? $cur_substr : $max_substr;
    } elsif ($ltr ne $top_ltr && $mismatch == 1) {
      $mismatch = 0;
      $cur_substr = 0;
    }
  }
  $low_comp = 1 if ($max_substr/$str_len >= 0.6); 
  return $low_comp;
}
