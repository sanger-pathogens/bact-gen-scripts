#! /usr/bin/perl -w
use strict;
use warnings;

my @features;
my @strain_1;
my @strain_2;
my @strain_3;
my @strain_4;
my $n = 0;
my $n_max;
my @patterns;
my %table;
my %hash;
my $x;
my $y;
my $i;
my $test;
my $strainsum;
my $totalsize;
my $numloci;
my %strains;
my $pres;
my $in;
my $abs;
my $out;

open IN, $ARGV[0] or die print "mauve2artemis.pl <Mauve backbone file> <ordered list of strains>\nConverts a Mauve backbone file to an Artemis readable tab file for each strain indicating conservation of LCB\n";
open STR, $ARGV[1], or die;

$i = 0;
while (<STR>) {
	chomp $_;  
	$strains{$i} = $_;
	$i++;
}

my $i_max = $i;
$i = 0;


while (<IN>) {
	unless (substr($_,0,1) =~ /s/) {
		my @loci = split(/\t/,$_);
		$numloci=scalar(@loci);
		$strainsum=0;
		$totalsize=0;
		for ($x =0; $x < $numloci; $x += 2){
			$y=$x/2;
			if ($loci[$x+1]-$loci[$x] == 0) {
				$hash{$y}[$n] = 0;
			} else {
				$hash{$y}[$n] = 1;
			}
			$strainsum+=$hash{$y}[$n];
			$totalsize+=$loci[$x+1];
			$totalsize-=$loci[$x];
		}
		$y = 0;
		while ($y <= $numloci/2) {
			if ($hash{$y}[$n] == 1) {
				open OUT, ">> $strains{$y}_mauve.tab";
				if ($loci[$y*2] >= 0) {
					print OUT "FT   misc_feature    $loci[$y*2]..$loci[($y*2)+1]\n";
				} else {
					my $num_a = abs($loci[$y*2]);
					my $num_b = abs($loci[(2*$y)+1]);
					print OUT "FT   misc_feature    $num_a..$num_b\n";
				}
				while ($i <= $i_max) {
					if ($hash{$i}[$n] == 1) {
						$pres++;
						$in = "$in $strains{$i}";
					} elsif ($hash{$i}[$n] == 0) {
						$abs++;
						$out = "$out $strains{$i}";
					}
					$i++;
				}
				if ($pres == $numloci/2) {
					print OUT 'FT                   /note='.'"'."Core genomic region".'"'."\n";
				} else {
					print OUT 'FT                   /note='.'"'."Present in$in".'"'."\n";
					print OUT 'FT                   /note='.'"'."Absent from$out".'"'."\n";
				}
				if ($pres <= 13) {
					print OUT "FT                   /colour=$pres\n";
				} else {
					print OUT "FT                   /colour=0\n";
				}
				$pres = 0;
				$abs = 0;
				$in = "";
				$out = "";
				$y++;
				$i = 0;
			} elsif ($hash{$y}[$n] == 0) {
				$y++;
				$i = 0;
			}
		}
	}
	$n++;
}
