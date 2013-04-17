#!/usr/bin/perl

# Perl script to calculate statistics on the infile. The infile should be the
# "intermediate" form: no duplicates, 14 fields, laid out as:
#
# str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 cigar1 seq1 mapq2 cigar2 seq2
#
# (Old form no fragments.)
# (Old form contained 10 fields with the cigar fields eliminated.)
#
# This intermediate form makes it easy to sort (to remove duplicates) and to
# grab just the first six fields to create the HiC map.
#
# The last six fields are needed for this statistics scripts, in order to
# determine the histogram of mapQ values, the number of dangling ends, and the
# number of ligation junctions.  Cigar is needed to properly calculate dups.
#
# The script also requires a HindIII restriction site file, which lists on
# each line, the sorted locations of the HindIII restriction sites.
#
# The script will print out the total reads, # dangling ends, # ligation
# junctions, #intra vs inter chromosomal reads, # inner/outer/right/left
# read pairs, a histogram of the MAPQ values, a histogram of the distance
# from the closest HindIII site, and a count of which end was closest to the
# HindIII site.
#
# Usage:	statistics.pl <infile>

use File::Basename;
use POSIX;
use List::Util qw[min max];
use Getopt::Std;
use vars qw/ $opt_s $opt_l $opt_d $opt_o /;

# Check arguments
getopts('s:l:d:o:');

my $site_file = "/broad/aidenlab/restriction_sites/hg19_HindIII.txt";
my $ligation_junction = "AAGCTAGCTT";
my $outdir;

if (scalar(@ARGV) == 1) {
	($infile) = @ARGV;
	($fname,$outdir,$suffix) = fileparse($infile);
	$stats_file = $outdir . "/stats.txt";
}
else {
	print "Usage: statistics.pl -s[site file] -d[outdir] -l[ligation] -o[stats file] <infile>\n";
	print " <infile>: file in intermediate format to calculate statistics on\n";
	print " [site file]: list of HindIII restriction sites, one line per chromosome (default HindIII hg19)\n";
	print " [outdir]: output directory for figures (default same directory as infile)\n";
	print " [ligation]: ligation junction (default $ligation_junction)\n";
	print " [stats file]: output file containing total reads, for library complexity (default infile directory/stats.txt)\n";
	exit;
}


if ($opt_s) {
  $site_file = $opt_s;
}
if ($opt_l) {
  $ligation_junction = $opt_l;
}
if ($opt_d) {
  $outdir = $opt_d;
}
if ($opt_o) {
	$stats_file = $opt_o;
}

my $dangling_junction = substr $ligation_junction, length($ligation_junction)/2;

# Global variables for calculating statistics
my %chromosomes;
my %hindIII;
my %mapQ;
my %mapQ_inter;
my %mapQ_intra;
my %innerM;
my %outerM;
my %rightM;
my %leftM;
my $three_prime_end=0;
my $five_prime_end=0;
my $total = 0;
my $dangling = 0;
my $ligation = 0;
my $inner = 0;
my $outer = 0;
my $left = 0;
my $right = 0;
my $inter = 0;
my $intra = 0;
my $small = 0;
my $large = 0;
my $very_small = 0;
my $very_small_dangling = 0;
my $small_dangling = 0;
my $large_dangling = 0;
my $inter_dangling = 0;
my $true_dangling_intra_small = 0;
my $true_dangling_intra_large = 0;
my $true_dangling_inter = 0;
my $total_current = 0;
# logspace bins
my @bins = (10,12,15,19,23,28,35,43,53,66,81,100,123,152,187,231,285,351,433,534,658,811,1000,1233,1520,1874,2310,2848,3511,4329,5337,6579,8111,10000,12328,15199,18738,23101,28480,35112,43288,53367,65793,81113,100000,123285,151991,187382,231013,284804,351119,432876,533670,657933,811131,1000000,1232847,1519911,1873817,2310130,2848036,3511192,4328761,5336699,6579332,8111308,10000000,12328467,15199111,18738174,23101297,28480359,35111917,43287613,53366992,65793322,81113083,100000000,123284674,151991108,187381742,231012970,284803587,351119173,432876128,533669923,657933225,811130831,1000000000,1232846739,1519911083,1873817423,2310129700,2848035868,3511191734,4328761281,5336699231,6579332247,8111308308,10000000000);

# read in restriction site file and store as multidimensional array
open FILE, $site_file or die $!;
while (<FILE>) {
	my @locs = split;
	my $key = shift(@locs);
	my $ref = \@locs;
	$chromosomes{$key} = $ref;
}
close(FILE);

# read in infile and calculate statistics
open FILE, $infile or die $!;
while (<FILE>) {
	$total_current++;
	my @record = split;

	my $num_records = scalar(@record);
	# position distance
	my $pos_dist = abs($record[2] - $record[6]);
	
	my $hist_dist = &bsearch($pos_dist,\@bins);	 
	
	my $is_dangling = 0;
	# one part of read pair has unligated end
	if ($num_records > 8 && ($record[10] =~ m/^$dangling_junction/ || $record[13] =~ m/^$dangling_junction/)) {
		$dangling++;
		$is_dangling=1;
	}
	
	# look at chromosomes
	if ($record[1] eq $record[5]) {
		$intra++;
		# determine right/left/inner/outer ordering of chromosomes/strands
		if ($record[0] == $record[4]) {
	    if ($record[0] == 0) {
				if ($pos_dist >= 20000) {
					$right++;
				}
				$rightM{$hist_dist}++;
	    }
	    else {
				if ($pos_dist >= 20000) {
					$left++;
				}
				$leftM{$hist_dist}++;
	    }
		}
		else {
	    if ($record[0] == 0) {
				if ($record[2] < $record[6]) {
					if ($pos_dist >= 20000) {
						$inner++;
					}
					$innerM{$hist_dist}++;
				}
				else {
					if ($pos_dist >= 20000) {
						$outer++;
					}
					$outerM{$hist_dist}++;
				}
	    }
	    else {
				if ($record[2] < $record[6]) {
					if ($pos_dist >= 20000) {
						$outer++;
					}
					$outerM{$hist_dist}++;
				}
				else {
					if ($pos_dist >= 20000) {
						$inner++;
					}
					$innerM{$hist_dist}++;
				}
	    }
		}
		# intra reads less than 20KB apart
		if ($pos_dist < 10) {
	    $very_small++;
	    if ($is_dangling) {
				$very_small_dangling++;
	    }
		}
		elsif ($pos_dist < 20000) {
	    $small++;
	    if ($is_dangling) {
				$small_dangling++;
	    }
		}
		else {
	    $large++;
	    if ($is_dangling) {
				$large_dangling++;
	    }
		}
	}
	else {
		$inter++;
		if ($is_dangling) {
	    $inter_dangling++;
		}
	}
	if ($num_records > 8) {
		my $mapq_val = min($record[8],$record[11]);
		if ($mapq_val <= 200) {
			$mapQ{$mapq_val}++;
			if ($record[1] eq $record[5]) {
				$mapQ_intra{$mapq_val}++;
			}
			else {
				$mapQ_inter{$mapq_val}++;
			}
		}
		# read pair contains ligation junction
		if ($record[10] =~ m/$ligation_junction/ || $record[13] =~ m/$ligation_junction/) {
			$ligation++;
		}
	}
	# determine distance from nearest HindIII site, add to histogram
	my $report = ($record[1] != $record[5]) || ($pos_dist >= 20000);
	my $dist = &distHindIII($record[0], $record[1], $record[2], $record[3], $report);
	if ($dist <= 2000) {
		$hindIII{$dist}++;
	}
	
	$dist = &distHindIII($record[4], $record[5], $record[6], $record[7], $report);
	if ($dist <= 2000) {
		$hindIII{$dist}++;
	}
	
	if ($is_dangling) {
		if ($record[10] =~ m/^$dangling_junction/) {
	    $dist = &distHindIII($record[0], $record[1], $record[2], $record[3], 1);
		}
		else { #	$record[13] =~ m/^$dangling_junction/) 
	    $dist = &distHindIII($record[4], $record[5], $record[6], $record[7], 1);
		}
		if ($dist == 1) {
	    if ($record[1] == $record[5]) {
				if ($pos_dist < 20000) {
					$true_dangling_intra_small++;
				}
				else {
					$true_dangling_intra_large++;
				}
	    }
	    else {
				$true_dangling_inter++;
	    }
		}
	}
}
close(FILE);

my($statsfilename, $directories, $suffix)= fileparse($stats_file, qr/\.[^.]*/);
my $histsfile = $directories . $statsfilename . "_hists.m";
my $statssupfile = $directories . $statsfilename . "_supp.txt";

open FILE, " >> $stats_file", or die $!;
print FILE "Total reads in current file: " . commify($total_current) . "\n";
if ($total_current==0) {
	$total_current=1;
}
printf FILE "Ligations: %s (%0.2f\%)\n", commify($ligation), $ligation*100/$total_current;
if ($five_prime_end + $three_prime_end > 0) {
	printf FILE "Five prime: %s (%0.2f\%)\n", commify($five_prime_end), $five_prime_end*100/($five_prime_end + $three_prime_end);
	printf FILE "Three prime: %s (%0.2f\%)\n", commify($three_prime_end), $three_prime_end*100/($five_prime_end + $three_prime_end);
}
else {
	printf FILE "Five prime: $five_prime_end (%0.2f\%)\n", 0;
	printf FILE "Three prime: $three_prime_end (%0.2f\%)\n", 0;
}
printf FILE "Inter: %s (%0.2f\%)\n", commify($inter), $inter*100/$total_current;
printf FILE "Intra: %s (%0.2f\%)\n", commify($intra), $intra*100/$total_current;
printf FILE "Small: %s (%0.2f\%)\n", commify($small), $small*100/$total_current;
printf FILE "Large: %s (%0.2f\%)\n", commify($large), $large*100/$total_current;
printf FILE "Very small: %s (%0.2f\%)\n", commify($very_small), $very_small*100/$total_current;
if ($large > 0) {
	printf FILE "Inner: %s (%0.2f\%) \n", commify($inner), $inner*100/$large;
	printf FILE "Outer: %s (%0.2f\%) \n", commify($outer), $outer*100/$large;
	printf FILE "Left: %s (%0.2f\%) \n", commify($left), $left*100/$large;
	printf FILE "Right: %s (%0.2f\%) \n", commify($right), $right*100/$large;
}
close FILE;

open FILE, " > $statssupfile", or die $!;
printf FILE "Dangling: %s (%0.2f\%)\n", commify($dangling), $dangling*100/$total_current;
printf FILE "  Very small: %s (%0.2f\%)\n", commify($very_small_dangling), $very_small_dangling*100/$total_current;
printf FILE "  Small: %s (%0.2f\%)\n", commify($small_dangling),$small_dangling*100/$total_current;
printf FILE "  Large: %s (%0.2f\%)\n", commify($large_dangling),$large_dangling*100/$total_current;
printf FILE "  Inter: %s (%0.2f\%)\n", commify($inter_dangling),$inter_dangling*100/$total_current;
printf FILE "  True Small: %s (%0.2f\%)\n", commify($true_dangling_intra_small),$true_dangling_intra_small*100/$total_current;
printf FILE "  True Large: %s (%0.2f\%)\n", commify($true_dangling_intra_large),$true_dangling_intra_large*100/$total_current;
printf FILE "  True Inter: %s (%0.2f\%)\n", commify($true_dangling_inter),$true_dangling_inter*100/$total_current;
close FILE;

open FILE, "> $histsfile", or die $!; 
print FILE "A = [\n";
for (my $i=1; $i <= 2000; $i++) {
	my $tmp =	 $hindIII{$i} || 0;
	print FILE "$tmp ";
}
print FILE "\n];\n";
print FILE "B = [\n";
for (my $i=0; $i <= 200; $i++) {
	my $tmp = $mapQ{$i} || 0;
	my $tmp2 = $mapQ_intra{$i} || 0;
	my $tmp3 = $mapQ_inter{$i} || 0;
	print FILE "$tmp $tmp2 $tmp3\n ";
}
print FILE "\n];\n";
print FILE "D = [\n";
for (my $i=0; $i < scalar(@bins); $i++) {
	my $tmp = $innerM{$i} || 0;
	print FILE "$tmp ";
	$tmp = $outerM{$i} || 0;
	print FILE "$tmp ";
	$tmp = $rightM{$i} || 0;
	print FILE "$tmp ";
	$tmp = $leftM{$i} || 0;
	print FILE "$tmp\n";
}

print FILE "\n];";
print FILE "x = [\n";
for (my $i=0; $i < scalar(@bins); $i++) {
	print FILE "$bins[$i] ";
}
print FILE "\n];\n";
close FILE;

# Find distance to nearest HindIII restriction site
sub distHindIII {
	# find upper index of position in sites array via binary search
	#my $index = &bsearch($_[2],$chromosomes{$_[1]});
	my $index = $_[3];
	# get distance to each end of HindIII fragment
	my $dist1;
	if ($index == 0) {
		# first fragment, distance is position
		$dist1 =	$_[2];	
	}
	else {
		$dist1 = abs($_[2] - $chromosomes{$_[1]}[$index-1]);
	}
	my $dist2 = abs($_[2] - $chromosomes{$_[1]}[$index]);
	
#	if ($index == 0) {
#		if (!($_[2] >= 0 && $_[2] <= $chromosomes{$_[1]}[$index])) {
#	    print(STDERR "Problem not " . $_[2] . " >= 0 && <= " . $chromosomes{$_[1]}[$index] . " " . $_[1] . " index = " . $index . " chr = " . $_[1] . "\n");
#		}
#	}
#	else {
#		if (!($_[2] >= $chromosomes{$_[1]}[$index-1] && $_[2] <= $chromosomes{$_[1]}[$index])) {
#	    print(STDERR "Problem not " . $_[2] . " >= " . $chromosomes{$_[1]}[$index-1] . " && <= " . $chromosomes{$_[1]}[$index] . " index = " . $index . " chr = " . $_[1] . "\n");
#		}
#	}
	# get minimum value -- if (dist1 <= dist2), it's dist1, else dist2
	my $retval = $dist1 <= $dist2 ? $dist1 : $dist2; 
	# get which end of the fragment this is, 3' or 5' (depends on strand)
	if ($retval == $dist1 && $_[4]) {
		$_[0] == 0 ? $five_prime_end++ : $three_prime_end++;
	}
	elsif ($retval == $dist2 && $_[4]) {
		$_[0] == 16 ? $five_prime_end++ : $three_prime_end++;
	}
	return $retval;
}

# Binary search, array passed by reference
# search array of integers a for given integer x
# return index where found or upper index if not found
sub bsearch {
	my ($x, $a) = @_;		 # search for x in array a
	my ($l, $u) = (0, @$a - 1);	 # lower, upper end of search interval
	my $i;	      		 # index of probe
	while ($l <= $u) {
		$i = int(($l + $u)/2);
		if ($a->[$i] < $x) {
	    $l = $i+1;
		}
		elsif ($a->[$i] > $x) {
	    $u = $i-1;
		}
		else {
	    return $i; # found
		}
	}
	return $l;				 # not found, return upper
}
sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}
