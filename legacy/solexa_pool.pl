#!/usr/local/bin/perl

=head1 NAME

solexa_pool_CI.pl

=head1 SYNOPSIS

Wrapper to run Maq, extract pileup for SNPs, calculate weighted and
unweighted 95% confidence interval for the number of strains in which 
each potential SNP is present.

=head1 AUTHORS

Kathryn Holt <kh2(at)sanger.ac.uk>

=head1 COPYRIGHT

=head1 BUGS

if you witness any bug please contact the authors

=head1 LICENSE

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

=cut




use Getopt::Long;



sub show_help {
    print STDERR <<EOF;

$0:

 	Wrapper to run Maq, extract pileup for SNPs, calculate weighted 
	and unweighted 95% confidence interval for the number of strains 
	in which each potential SNP is present.

Usage: $0 -i reads -r reference -n numstrains -q maxqual -p prefix -id poolid -opt maqoptions [-rep 1]

	-i   Reads file (fastq)
	-r   Reference sequence file (fasta)
	-n   Number of strains in the pool
	-q   Maximum quality score for calibrated qualities
	-p   Prefix for output files on this set of pools
	-id  Unique identifier for this pool
	-opt Strain of additional options to pass to Maq 
	-rep 0 to switch off running the repeat filter [optional]

EOF
}

	&GetOptions(	"i=s"=>\$reads_file,
			"r=s"=>\$reference,
			"n=s"=>\$numstrains,
			"q=s"=>\$maxqual,
			"p=s"=>\$prefix,
			"id=s"=>\$poolid,
			"opt=s"=>\$maq_options,
			"rep=s"=>\$filter,
			"h"=>\$help,
			"help"=>\$help);
	
	if( $help ) {
    		&show_help();
    		exit(1);
	}
	
	if (!$filter) {$filter=1;}
	
	#### run Maq programs
	
	# SETTING: location of easy run
	$easyrun = "maq.pl easyrun";
	$maq = "maq";
	
	## run easy run
	
	$dir = $prefix.$poolid;
	print "running maq easyrun:\n  $easyrun -N $numstrains -E 0 -1 35 -d $dir $reference $reads_file\n";
	system `$easyrun -N $numstrains -E 0 -d $dir $reference $reads_file`;

	# filter SNP calls from non-unique regions

	$rawsnp = "$dir/cns.snp";
    	$filtersnp = "$dir/cns.snp.filtered";

    if ($filter>0) {
	print "filtering SNP calls from non-unique regions\n";
	open(IN,$rawsnp) or die("Couldn't open SNP file $rawsnp\n");
	open(OUT,">$filtersnp");
	while ($line=<IN>) {
		my @fields = split /\t/, $line;
		if ($fields[6]<=1) {
			print OUT $line;
		}
	}
	close OUT;
	close IN;
    } # end filter
    
    else {system `cp $rawsnp $filtersnp`;}
	
	# generate pileup
	
	$pileup = "$dir/cns.snp.filtered.pileup";
	print "running maq pileup\n  $maq pileup -v -l $filtersnp $dir/ref.bfa $dir/all.map > $pileup\n";
	system `$maq pileup -v -l $filtersnp $dir/ref.bfa $dir/all.map > $pileup`;
	
	### generate estimates
	
	print "generating strain frequency estimates and confidence intervals\n";
	
	@base_list = ("A","C","G","T");
	
	$outfile=$prefix.$poolid."frequencies.txt";
	open (OUT, ">$outfile");
	print OUT "Pool\tPosition\tReferenceAllele\tReads\tSameSNP\tSNPAllele\tProportion\tStDev\tNumStrains\tNumStrainsInt\tLowerCI\tUpperCI\tSNPAlleleQ\tProportionQ\tStDevQ\tNumStrainsQ\tNumStrainsIntQ\tLowerCIQ\tUpperCIQ\n";
	open (IN,$pileup) or die ("Couldn't open pileup file $pileup");
	while ($line=<IN>) {
		chomp $line;
		my @fields = split /\t/, $line;
		my $pos = $fields[1];
		my $ref = $fields[2];
		my $N = $fields[3];
		my @bases = split//,$fields[4]; # 0 position contains @ symbol
		my @quals = split//,$fields[5]; # 0 position contains @ symbol
		my $num_bases = $#bases;
		if ($N ne $num_bases) {
			print "uh oh. depth should be $N but only found $num_bases reads!\n";
		}
		elsif ($N>$mindepth) {
			print OUT "$prefix$poolid\t$pos\t$ref\t$N\t";
			my @weights; #0=a, 1=c, 2=g, 3=t, 4=ref
			my @weights_sq; #0=a, 1=c, 2=g, 3=t, 4=ref
			my @hits; #0=a, 1=c, 2=g, 3=t, 4=ref
			foreach $b (1..$num_bases) {
				my $base = $bases[$b];
				my $qual = ord($quals[$b])-33; # convert ASCII to number
				my $weight = $qual/$maxqual;#(($qual-20)^2)/(($maxqual-20)^2); # weight of quality score
				#if ($qual<20) {$weight=0;}
				#elsif ($qual>37) {$weight=1;}
				my $i;
				if ($base =~ /[aA]/) {
					$i=0
				}
				elsif ($base =~ /[cC]/) {
					$i=1;
				}
				elsif ($base =~ /[gG]/) {
					$i=2;
				}
				elsif ($base =~ /[tT]/) {
					$i=3;
				}
				else { $i=4; }
				$weights[$i] += $weight;
				$weights_sq[$i] += ($weight*$weight);
				$hits[$i]++;
			}
			## calculate proportion of qualities for each potential SNP
			my @p; #0=a, 1=c, 2=g, 3=t
			my @p_weighted; #0=a, 1=c, 2=g, 3=t
			my $sum_of_weights = $weights[0] + $weights[1] + $weights[2] + $weights[3] + $weights[4];
			my $sum_of_sq_weights = $weights_sq[0] + $weights_sq[1] + $weights_sq[2] + $weights_sq[3] + $weights_sq[4];
			my $total_hits = $hits[0] + $hits[1] + $hits[2] + $hits[3] + $hits[4];
			foreach $i (0..3) {
				if ($hits[$i]+$hits[4]>0) {
					$p[$i] = sprintf("%.3f",$hits[$i]/$total_hits);
				}
				else { $p[$i] = 0; }
				if ($weights[$i]+$weights[4]>0) {
					$p_weighted[$i] = sprintf("%.3f",$weights[$i]/$sum_of_weights);
				}
				else { $p_weighted[$i] = 0;}
			}
			## determine which base is the most likely SNP, using raw proportion
			if ($p[0] > $p[1] && $p[0] > $p[2] && $p[0] > $p[3] ) { $snp=0; }
			elsif ($p[1] > $p[0] && $p[1] > $p[2] && $p[1] > $p[3] ) { $snp=1; }
			elsif ($p[2] > $p[0] && $p[2] > $p[1] && $p[2] > $p[3] ) { $snp=2; }
			else { $snp=3; }
			## determine which base is the most likely SNP, using weighted proportion
			if ($p_weighted[0] > $p_weighted[1] && $p_weighted[0] > $p_weighted[2] && $p_weighted[0] > $p_weighted[3] ) { $snp_weighted=0; }
			elsif ($p_weighted[1] > $p_weighted[0] && $p_weighted[1] > $p_weighted[2] && $p_weighted[1] > $p_weighted[3] ) { $snp_weighted=1; }
			elsif ($p_weighted[2] > $p_weighted[0] && $p_weighted[2] > $p_weighted[1] && $p_weighted[2] > $p_weighted[3] ) { $snp_weighted=2; }
			else { $snp_weighted=3; }
			## determine if they predict the same SNP
			if ($snp eq $snp_weighted) {$match=1;}
			else {$match=0;}
			print OUT $match;
			## calculate variance
			my $var = sprintf("%.3f",sqrt(($p[$snp]*(1-$p[$snp]))/$total_hits));
			my $lower = sprintf("%.3f",($p[$snp]-1.96*$var)*$numstrains);
			my $upper = sprintf("%.3f",($p[$snp]+1.96*$var)*$numstrains);
			my $numstrains_p = sprintf("%.3f",$p[$snp]*$numstrains);
			my $numstrains_p_int = sprintf("%.0f",$p[$snp]*$numstrains);
			print OUT "\t$base_list[$snp]\t$p[$snp]\t$var\t$numstrains_p\t$numstrains_p_int\t$lower\t$upper\t";
			## calculate weighted variance
			my $numerator = ($p_weighted[$snp_weighted]*(1-$p_weighted[$snp_weighted])*$sum_of_sq_weights);
			my $denominator = $sum_of_weights*$sum_of_weights;
			my $var_weighted;
			if ($denominator > 0) {
				$var_weighted = sprintf("%.3f",sqrt($numerator/$denominator));
			}
			else {
				$var_weighted = 0;
			}
			my $lower_weighted = sprintf("%.3f",($p_weighted[$snp_weighted]-1.96*$var_weighted)*$numstrains);
			my $upper_weighted = sprintf("%.3f",($p_weighted[$snp_weighted]+1.96*$var_weighted)*$numstrains);
			my $numstrains_p_weighted = sprintf("%.3f",$p_weighted[$snp_weighted]*$numstrains);
			my $numstrains_p_weighted_int = sprintf("%.0f",$p_weighted[$snp_weighted]*$numstrains);
			print OUT "$base_list[$snp_weighted]\t$p_weighted[$snp_weighted]\t$var_weighted\t$numstrains_p_weighted\t$numstrains_p_weighted_int\t$lower_weighted\t$upper_weighted\t\n";
		}
	}
	close IN;
