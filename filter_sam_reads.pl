#!/usr/bin/env perl 

## use the table of mapped reads and max AS to output only the best reads with appropriate NH tag

use strict;
use warnings;
use Data::Dumper; 

if (scalar @ARGV != 2) {
  print "USAGE: ./filter_sam_reads.pl <name-sorted sam> <count table>\n"; 
  exit 1;
}

my $sam = shift @ARGV; 
my $counts = shift @ARGV; 
open SAM,"<",$sam or die "$!"; 
open COUNTS,"<",$counts or die "$!"; 

my $R = {}; 

while (<COUNTS>) { 
	my $rname = (split /\t/)[0];
	my $as = (split /\t/)[1];
	my $nh = (split /\t/)[2];
	$R->{$rname}->{AS} = $as; 
	$R->{$rname}->{NH} = $nh; 
}

print STDERR Dumper($R); 

while (<SAM>) {
	if (m/^@/) {
		print; 
	} else {
		chomp; 
	  my $rname = (split /\t/)[0]; 
	  m/\tAS:i:(\d+)\t/; 
	  my $as = $1;
		printf "%s\tNH:i:%d\n",$_,$R->{$rname}->{NH} if ($as == $R->{$rname}->{AS}); 
	} 
} 

close COUNTS;
close SAM; 
		


