#!/usr/bin/env perl 

## make UMI-style single-end fastq file from STARsolo or CellRanger output. 
## this does rely on some systems calls, so have gzip/samtools/etc installed
## (need to think about how this would work with Singularity)

use strict; 
use warnings;

if (scalar @ARGV != 1) { 
	print STDERR "USAGE: ./make_umi_fastq.pl <STARsolo or Cell/Space Ranger directory>\n"; 
	exit 1; 
} 

my $dir = shift @ARGV; 

if (-d $dir."/output" && -e $dir."/Log.final.out" && -e $dir."/Unmapped.out.mate1.gz" && $dir."//Unmapped.out.mate2.gz" ) { 
	print STDERR "Sample was determined to be STARsolo output! \nMaking a single-end, UMI-tools-formatted fastq file..\n";
	
	## we need to check if the file was mapped as paired-end; this changes what R1 and R2 actually are.
	## if $paired == 1, mate1 contains the barcode + UMI; if $paired == 0, it's mate2.  
	my $paired = 0;
	open LOG,"<",$dir."/Log.out" or die "$!"; 
	while (<LOG>) { 
		$paired = 1 if (m/clip5pNbases 39 0/); 
	}
	close LOG;
	if ($paired) { 
		print STDERR "Sample was determined to be 5' 10x paired-end experiment!\n";
	} else { 
		print STDERR "Sample was determined to be 3'/5' 10x single-end experiment!\n"; 
	} 

	## now do the processing. We keep only the BC+UMI part of the barcode read
	## this would lose some biological info in case of PE 5' experiments, but I think it should not matter this much for bacterial reads
	my $R1 = ($paired) ? $dir."/Unmapped.out.mate1.gz" : $dir."/Unmapped.out.mate2.gz"; 
	my $R2 = ($paired) ? $dir."/Unmapped.out.mate2.gz" : $dir."/Unmapped.out.mate1.gz"; 
  print STDERR "Using files:\nR1 = $R1\nR2 = $R2\n"; 	
  ## this occasionally stupid concoction reads 2 fastq files synchronously. 
	open R1,"zcat $R1 |" or die "$!"; 
	open R2,"zcat $R2 |" or die "$!"; 
	my $nr = 1; 
  my ($rname,$bc,$umi,$seq,$qual) = ('') x 5;
	
	while (<R1>) {
		my $r1 = $_; 
		chomp $r1;
		my $r2 = <R2>;
		chomp $r2;
		#print STDERR "Line number $nr, r1 = $r1, r2 = $r2\n";

		if ($nr % 4 == 1) {
		  $r1 =~ m/(.*?)\s/;                                     ## shortest string until space 
			$rname = $1;
			$r2 =~ m/(.*?)\s/; 
			my $rname2 = $1; 
			print STDERR "WARNING: Read names do not match: $rname =/= $rname2!\n" if ($rname ne $rname2);
		} elsif ($nr % 4 == 2) { 
			my $len = length $r1; 
			my $bclen = ($len > 24) ? 16 : 14; 
			my $umilen = ($len > 28) ? 10 : $len - $bclen;         ## it's 10 for 5' PE, 10 for 3' v1(24-14) or v2(26-16), and 12 for 3'v3 (28-16)
			$bc = substr $r1,0,$bclen; 
			$umi = substr $r1,$bclen,$umilen;
			$seq = $r2;
		} elsif ($nr % 4 == 0) { 
			$qual = $r2;
			printf "%s_%s_%s\n%s\n+\n%s\n",$rname,$bc,$umi,$seq,$qual;
		} 
		$nr++; 
	}
  close R1;
  close R2;
} else { 
  ## here we assume the presence of a Cell/Space Ranger/STARsolo-formatted BAM file: 
	## single- or paired-end (meaning biological reads), with CR and UR set for all reads, and CB/UB for reads that mapped successfully
	## we will take only unmapped reads (-f4), and only R2 for paired-end reads (-f128)
	## since all the BAM files are named differently, we would need to do a little sing-and-a-dance first: 
	my $bam = `find $dir/* | grep -iv atac | grep \"\.bam\$\"`; 
	chomp $bam; 
	if (length $bam) { 
		print STDERR "Sample is assumed to be Cell/Space Ranger/STARsolo BAM! \nMaking a single-end, UMI-tools-formatted fastq file..\n";
	} else { 
		print STDERR "ERROR: failed to find a BAM file in the directory provided!\n"; 
		exit 1; 
	}
	print STDERR "Using files: \nBAM = $bam\n";

	## we need to check if the file was mapped as paired-end; this changes what R1 and R2 actually are.
	## if $paired == 1, flags 64/128 will be set (R1 will be trimmed for BC+UMI).  
	my $preads = `samtools view -h $bam | head -10000 | samtools flagstat - | grep \"properly paired\" | cut -d' ' -f1`; 
	my $paired = ($preads > 0) ? 1 : 0; 
	if ($paired) { 
		print STDERR "Sample was determined to be 5' 10x paired-end experiment!\n";
	} else { 
		print STDERR "Sample was determined to be 3'/5' 10x single-end experiment!\n"; 
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
    printf "@%s_%s_%s\n%s\n+\n%s\n",$rname,$bc,$umi,$seq,$qual;
	} 
	
	close BAM;
} 
