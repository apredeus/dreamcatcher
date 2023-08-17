#!/usr/bin/env perl 

## run this on a name-sorted SAM file to get the best alignment score (AS), 
## number of reads with the best score, and the overall number of alignments

use strict; 
use warnings; 
#use Data::Dumper; 

if (scalar @ARGV != 1) {
	print "USAGE: ./count_read_occ.pl <name-sorted sam>\n"; 
	exit 1;
}

my $sam = shift @ARGV; 
open SAM,"<",$sam or die "$!"; 

my $prev_read;
my $curr_read;
my $R = {}; ## hash to store read 
my $count = 1; 

while (<SAM>) { 
	next if (m/^@/); 
	my $curr_read = (split /\t/)[0];  
	$prev_read = $curr_read if (! defined $prev_read);
  
	if ($prev_read eq $curr_read) { 
		## still counting reads
		m/\tAS:i:(\d+)\t/; 
		my $as = $1;
		$R->{$count} = $as;
		$count++; 
	} else {                        
		## time to print all the info for the previous read
	  my ($max_as,$max_count) = (0) x 2;
		## find the highest AS for this read
		foreach my $c (keys %{$R}) {
			$max_as = $R->{$c} if ($R->{$c} > $max_as); 
		}

		## count how many SAM occurrences of this read have the max AS
		foreach my $c (keys %{$R}) {
			$max_count++ if ($R->{$c} == $max_as); 
		}
		printf "%s\t%d\t%d\t%d\n",$prev_read,$max_as,$max_count,scalar(keys %{$R}); 

	  ## reset the counter and the hash
		$count = 1; 
		$R = {}; 
		
		## now, populate things for the new read
		m/\tAS:i:(\d+)\t/; 
		my $as = $1;
		$R->{$count} = $as;
		$count++;
	}
	$prev_read = $curr_read; 
} 

## after the loop is done, we need to print one final record
my ($max_as,$max_count) = (0) x 2;
## find the highest AS for this read
foreach my $count (keys %{$R}) {
	$max_as = $R->{$count} if ($R->{$count} > $max_as); 
} 
    
## count how many SAM occurrences of this read have the max AS  
foreach my $count (keys %{$R}) {
  $max_count++ if ($R->{$count} == $max_as);
} 

printf "%s\t%d\t%d\t%d\n",$prev_read,$max_as,$max_count,scalar(keys %{$R});

close SAM; 

