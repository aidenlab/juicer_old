#!/usr/bin/perl 

# Perl script to convert to fragment map from infile. The infile should be the
# "intermediate" form: no duplicates, 10 fields, laid out as:
#
# strand1 chr1 pos1 strand2 chr2 pos2 mapq1 sequence1 mapq2 sequence2 
#
# The script also requires a restriction site file, which lists on 
# each line, the sorted locations of the enzyme restriction sites.
#
# Usage:  fragment.pl <infile>

use POSIX;

$site_file = "/broad/aidenlab/restriction_sites/hg19_DpnII.txt";
# Check arguments
if (scalar(@ARGV) == 1) {
  ($infile) = @ARGV;
}
elsif (scalar(@ARGV) == 2) {
  ($infile,$site_file) = @ARGV;
}
else {
  print "Usage: fragment.pl <infile> [site file]\n";
  print " <infile>: file in intermediate format to calculate statistics on\n";
  print " [site file]: list of restriction sites, one line per chromosome (default DpnII hg19)\n";
  exit;
}

# Global variables for calculating statistics
my %chromosomes;
my %hindIII;

# read in restriction site file and store as multidimensional array
open FILE, $site_file or die $!;
while (<FILE>) {
  my @locs = split;
  my $key = shift(@locs);
  my $ref = \@locs;
  $chromosomes{$key} = $ref;
	if ($key == "14") {
		$chromosomes{$key."m"} = $ref;  
		$chromosomes{$key."p"} = $ref;
	}
}
close(FILE);

# read in infile and calculate statistics
open FILE, $infile or die $!;
while (<FILE>) {
  my @record = split;

  # find upper index of position in sites array via binary search
  my $index1 = &bsearch($record[2],$chromosomes{$record[1]});
  my $index2 = &bsearch($record[5],$chromosomes{$record[4]});
  print $record[0] . " " . $record[1] . " " . $record[2] . " " . $index1 . " " . $record[3] . " " . $record[4] . " " . $record[5] . " " . $index2 . " ";
	for (my $i=6; $i < scalar(@record); $i++){
		print $record[$i] . " ";
	}
	print "\n";
}
close(FILE);

# Binary search, array passed by reference
# search array of integers a for given integer x
# return index where found or upper index if not found
sub bsearch {
  my ($x, $a) = @_;            # search for x in array a
  my ($l, $u) = (0, @$a - 1);  # lower, upper end of search interval
  my $i;                       # index of probe
  while ($l <= $u) {
    $i = int(($l + $u)/2);
    if ($a->[$i] < $x) {
      $l = $i+1;
    }
    elsif ($a->[$i] > $x) {
      $u = $i-1;
    } 
    else {
      return $i+1; # found
    }
  }
  return $l;         
}
