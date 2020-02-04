#!/usr/local/bin/perl

=head1 NAME

run_mummer_contigs.pl

=head1 SYNOPSIS

- get list of regions that don't map in one or more sequences
- get list of regions within X bp of contig ends in any sequence
- output table of condensed regions to exclude (use generate_embl.pl to view in artemis)

- call show-snps, write out alignments (20bp either side) & table for SNP calls that are
	- > X [20] bp from nearest contig end
	- maximal quality of basecall
	
- write alignments & table for SNP calls that are
	- <= X [2] bp from next nearest SNP/indel called
	- 20bp either side contain gaps

* Could look for regions that align with multiple sequences in the ref as well...

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

$0: Run MUMmer on published genome (ie no quality data, no contigs, just a finished sequence)

Usage: $0 -r reference_sequence -s query_sequence -o output [-d dist]

-d dist is the minimum number of bases to the next SNP/indel call, default is 2

Note: you must be logged into holly to use MUMmer at Sanger

EOF
}

	&GetOptions(	"s=s"=>\$sequence,
			"r=s"=>\$refseq,
			"d=s"=>\$distance_to_nearest,
			"o=s"=>\$output,
			"h"=>\$help,
			"help"=>\$help);
	
	if( $help ) {
    		&show_help();
    		exit(1);
	}
	
	
	# set defaults
	
	if (!$distance_to_nearest) {$distance_to_nearest = 2;}
	
	### CONTIG FILTERING
	
	# RUN MUMMER
	
	print "\n";
	
	# run nucmer
	$nucmer_command = "/nfs/pathsoft/external/linux/MUMmer3.19/scripts/nucmer";
        $nucmer_command .= " -prefix=$output";
        $nucmer_command .= " $refseq";
        $nucmer_command .= " $sequence";
        print $nucmer_command."\n\n";
        system $nucmer_command;
        $nucmer_output = $output.".delta";
		
	# run filter
	$filter_output = $output.".filter";
	$filter_command = "/nfs/pathsoft/external/linux/MUMmer3.19/delta-filter -r -q ";
	$filter_command .= " $nucmer_output > $filter_output";
	print "\n".$filter_command."\n\n";
	system $filter_command;
	
	# show snps
	
	print "  calling SNPs...\n\n";
	$snps_output = $output.".snps";
	$snps_command = "/nfs/pathsoft/external/linux/MUMmer3.19/show-snps";
	$snps_command .= " -Clr -T -H -x 20 $filter_output > $snps_output"; 
	print $snps_command."\n";
	system $snps_command;
	
	print "\nFiltering SNPs...\n  accepted SNPs are:\n      - no gaps within 10bp, and > $distance_to_nearest from nearest SNP/indel\n";
	print "  alignments are also printed for SNPs that are close to other SNPs or indels, so they may be checked manually\n\n";

	open (SNPS, $snps_output);
	while (<SNPS>) {
		chomp $_;
		$row = $_;
		@fields = split /\t/,$row;
		# read relevant fields
		$refpos = $fields[0];
		$refnt = $fields[1];
		$qpos = $fields[3];
		$qnt = $fields[2];
		if ($refnt =~ /[ACGT]/) { # make sure it is a SNP not a gap we are recording
			if ($qnt =~ /[ACGT]/) { # make sure it is a SNP not a gap we are recording
				$buff = $fields[4]; # distance to nearest SNP/indel
				$dist = $fields[5]; # distance to nearest sequence end
				$ref_context = $fields[8]; # 20bp either side of SNP from reference sequence
				$query_context = $fields[9]; # 20bp either side of SNP from reference sequence
				$ref_seq_name = $fields[12]; # contig query matched to
				$query_contig = $fields[13]; # contig query matched to
					if ( ($query_context =~ /[.]/) || ($ref_context =~ /[.]/) || ($buff <= $distance_to_nearest)) {
						$possible_snpcount++;
						$possibles .= $row."\n"; # add SNP data to possibles table
						# add the alignment to alignment file for checking
						$possible_alignments .= ">$ref_seq_name $refpos $refnt\n    $ref_context\n";
						$possible_alignments .= "    ";
						@ref_aligns = split //,$ref_context;
						@query_aligns = split //,$query_context;
						foreach $base (0..40) {
							if ($ref_aligns[$base] eq $query_aligns[$base]) {$possible_alignments .= "|";}
							else {$possible_alignments .= " "};
						}
						$possible_alignments .= "\n    $query_context\n>$query_contig $qpos $qnt    $buff  $dist\n\n";
					}
					else { # accept the SNPs
						$snpcount++;
						$table .= $row."\n"; # add SNP data to table
						# add the alignment to file
						$alignments .= ">$ref_seq_name $refpos $refnt\n    $ref_context\n";
						$alignments .= "    ";
						@ref_aligns = split //,$ref_context;
						@query_aligns = split //,$query_context;
						foreach $base (0..40) {
							if ($ref_aligns[$base] eq $query_aligns[$base]) {$alignments .= "|";}
							else {$alignments .= " "};
						}
						$alignments .= "\n    $query_context\n>$query_contig $qpos $qnt    $buff  $dist\n\n";
					}
			}
		}
	}
	close SNPS;
	
	# print alignments
	$alignments_output = $output.".aligns";
	print "  printing $snpcount alignments to $alignments_output\n";
	open (ALIGNS, "> $alignments_output");
	print ALIGNS  $alignments;
	close ALIGNS;
	
	# print table;
	$snp_table_output = $output.".txt";
	print "  printing $snpcount SNPs to $snp_table_output\n\n";
	open (SNPTABLE, "> $snp_table_output");
	print SNPTABLE  $table;
	close SNPTABLE;
	
	if ($possible_snpcount) {
	
	# print possible alignments for checking
	$alignments_output = $output."_check.aligns";
	print "  printing $possible_snpcount alignments to $alignments_output for checking\n";
	open (ALIGNS, "> $alignments_output");
	print ALIGNS  $possible_alignments;
	close ALIGNS;
	
	# print table;
	$snp_table_output = $output."_check.txt";
	print "  printing $possible_snpcount SNPs to $snp_table_output for checking\n\n";
	open (SNPTABLE, "> $snp_table_output");
	print SNPTABLE  $possibles;
	close SNPTABLE;
	
	}

END:

