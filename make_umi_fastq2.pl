#!/usr/bin/env perl 

## make UMI-style single-end fastq file from STARsolo or CellRanger output. 
## this does rely on some systems calls, so have gzip/samtools/etc installed
## (need to think about how this would work with Singularity)

use strict; 
use warnings;

if (scalar @ARGV != 2) { 
	print STDERR "USAGE: ./make_umi_fastq.pl <STARsolo or Cell/Space Ranger directory> <output file>\n"; 
	exit 1; 
} 

my $dir = shift @ARGV;
my $out = shift @ARGV; 

open OUT,">",$out or die "$!"; 

if (-d $dir."/output" && -e $dir."/Log.final.out") { 
	print STDOUT "make_umi_fastq.pl: sample was determined to be STARsolo output! \nMaking a single-end, UMI-tools-formatted fastq file..\n";
	
	## we need to check if the file was mapped as paired-end; this changes what R1 and R2 actually are.
	## if $paired == 1, mate1 contains the barcode + UMI; if $paired == 0, it's mate2.  
	my $paired = 0;
	open LOG,"<",$dir."/Log.out" or die "$!"; 
	while (<LOG>) { 
		$paired = 1 if (m/clip5pNbases 39 0/); 
	}
	close LOG;
	if ($paired) { 
		print STDOUT "make_umi_fastq.pl: sample was determined to be 5' 10x paired-end experiment!\n";
	} else { 
		print STDOUT "make_umi_fastq.pl: sample was determined to be 3'/5' 10x single-end experiment!\n"; 
	} 

	## now do the processing. We keep only the BC+UMI part of the barcode read
	## this would lose some biological info in case of PE 5' experiments, but I think it should not matter this much for bacterial reads
	my ($R1,$R2) = ('') x 2;
	if (-e $dir."/Unmapped.out.mate1.gz" && -e $dir."/Unmapped.out.mate2.gz") { 
	  $R1 = ($paired) ? $dir."/Unmapped.out.mate1.gz" : $dir."/Unmapped.out.mate2.gz"; 
	  $R2 = ($paired) ? $dir."/Unmapped.out.mate2.gz" : $dir."/Unmapped.out.mate1.gz";
	} elsif (-e $dir."/Unmapped.out.mate1.bz2" && -e $dir."/Unmapped.out.mate2.bz2") {
	  $R1 = ($paired) ? $dir."/Unmapped.out.mate1.bz2" : $dir."/Unmapped.out.mate2.bz2"; 
	  $R2 = ($paired) ? $dir."/Unmapped.out.mate2.bz2" : $dir."/Unmapped.out.mate1.bz2";
	} else { 
		print STDERR "make_umi_fastq.pl: cannot find unmapped read files Unmapped.out.mate(1/2).<gz/bz2> in the STARsolo directory! Exiting .."; 
		exit 1; 
	}

  print STDOUT "Using unmapped STARsolo files:\nR1 = $R1\nR2 = $R2\n"; 	
  ## this occasionally stupid concoction reads 2 fastq files synchronously.
	if ($R1 =~ m/gz$/ && $R2 =~ m/gz$/) { 
		open READ1,"pigz -cd $R1 |" or die "$!"; 
		open READ2,"pigz -cd $R2 |" or die "$!"; 
	} else { 
		open READ1,"lbzip2 -cd $R1 |" or die "$!"; 
		open READ2,"lbzip2 -cd $R2 |" or die "$!";
	}

	my $nr = 1; 
  my ($rname,$bc,$umi,$seq,$qual) = ('') x 5;
	my ($passed,$discarded) = (0) x 2; 
	
	while (<READ1>) {
		my $r1 = $_; 
		chomp $r1;
		my $r2 = <READ2>;
		chomp $r2;

		if ($nr % 4 == 1) {
		  $r1 =~ m/(.*?)\s/;                                     ## shortest string until space 
			$rname = $1;
			$r2 =~ m/(.*?)\s/; 
			my $rname2 = $1; 
			print STDERR "make_umi_fastq.pl: WARNING - read names do not match: $rname =/= $rname2!\n" if ($rname ne $rname2);
		} elsif ($nr % 4 == 2) { 
			my $len = length $r1; 
			my $bclen = ($len > 24) ? 16 : 14; 
			my $umilen = ($len > 28) ? 10 : $len - $bclen;         ## it's 10 for 5' PE, 10 for 3' v1(24-14) or v2(26-16), and 12 for 3'v3 (28-16)
			print STDERR "make_umi_fastq.pl: WARNING - UMI length is not equal to 10 or 12!\n" if ($umilen != 10 && $umilen != 12); 
			$bc = substr $r1,0,$bclen; 
			$umi = substr $r1,$bclen,$umilen;
			$seq = $r2;
		} elsif ($nr % 4 == 0) { 
			$qual = $r2;
			if (! low_complexity($seq)) { 
			  printf OUT "%s_%s_%s\n%s\n+\n%s\n",$rname,$bc,$umi,$seq,$qual;
			  $passed++; 
			} else {
			  #printf STDOUT "DEBUG: filtered low complexity sequence $seq\n"; 
		      $discarded++; 
			}
		} 
		$nr++; 
	}
  close READ1;
  close READ2;
	print STDOUT "make_umi_fastq.pl: making UMI-tools formatted reads: outputted $passed reads, discarded $discarded reads because of 50nt+ homopolymers ..\n"; 
} else { 
	## TODO: add read counter, add homopolymer filter for Cell Ranger BAM input

	## here we assume the presence of a Cell/Space Ranger/STARsolo-formatted BAM file: 
	## single- or paired-end (meaning biological reads), with CR and UR set for all reads, and CB/UB for reads that mapped successfully
	## we will take only unmapped reads (-f4), and only R2 for paired-end reads (-f128)
	## since all the BAM files are named differently, we would need to do a little sing-and-a-dance first: 
	my $bam = `find $dir/* | grep -iv atac | grep \"\.bam\$\"`; 
	chomp $bam; 
	if (length $bam) { 
		print STDOUT "make_umi_fastq.pl: sample is assumed to be Cell/Space Ranger/STARsolo BAM! \nMaking a single-end, UMI-tools-formatted fastq file..\n";
	} else { 
		print STDERR "make_umi_fastq.pl: ERROR - failed to find a BAM file in the directory provided! Exiting .. \n"; 
		exit 1; 
	}
	print STDOUT "Using files: \nBAM = $bam\n";

	## we need to check if the file was mapped as paired-end; this changes what R1 and R2 actually are.
	## if $paired == 1, flags 64/128 will be set (R1 will be trimmed for BC+UMI).  
	my $preads = `samtools view -h $bam | head -10000 | samtools flagstat - | grep \"properly paired\" | cut -d' ' -f1`; 
	my $paired = ($preads > 0) ? 1 : 0; 
	if ($paired) { 
		print STDOUT "make_umi_fastq.pl: sample was determined to be 5' 10x paired-end experiment!\n";
	} else { 
		print STDOUT "make_umi_fastq.pl: sample was determined to be 3'/5' 10x single-end experiment!\n"; 
	} 

	## now do the processing. We keep only the BC+UMI part of the barcode read
	## this would lose some biological info in case of PE 5' experiments, but I think it should not matter this much for bacterial reads
	if ($paired) { 
		open BAM,"samtools view -\@4 -f132 $bam |" or die "ERROR: failed to open the bam file using samtools!";
	} else { 
		open BAM,"samtools view -\@4 -f4 $bam |" or die "ERROR: failed to open the bam file using samtools!";
	}

	while (<BAM>) { 
    my ($rname,$bc,$umi,$seq,$qual) = ('') x 5;
    $rname = (split /\t/)[0];
    $seq   = (split /\t/)[9];
    $qual  = (split /\t/)[10];
    m/CR:\w:(\w+)\t.*UR:\w:(\w+)\t/;
    $bc = $1;
    $umi = $2;
    printf OUT "@%s_%s_%s\n%s\n+\n%s\n",$rname,$bc,$umi,$seq,$qual;
	} 
	
	close BAM;
} 

close OUT; 

sub low_complexity {
  my $low_comp = 0; 
  ## identify low complexity sequences of few kinds:
  ## (1) over 90% of the same letter overall;
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

  ## now calculate fraction of the top letter & set low_comp to 1 if it's 90%+ 
  my $str_len = length $str;
  $low_comp = 1 if ($top_count/$str_len >= 0.9);

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
