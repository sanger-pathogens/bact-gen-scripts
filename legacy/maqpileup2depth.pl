#!/usr/bin/perl -w

=head1 NAME

maqpileup2depth.pl

=head1 SYNOPSIS

Takes a maq pileup file as input, get depth information out of it.
The output can be loaded as user plot in Artemis, output contains read depths for different strands.

=head1 AUTHORS

Miao He  <mh10(at)sanger.ac.uk>

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


use strict;

sub show_help {
	die <<EOF;

	$0: take maq pileup file as input, get depth information out of it.
	The output can be loaded as user plot in artemis. 
	Output contains depths for different strands. 

	Usage: $0 maq_pileup_file > output


EOF
}

my $infile = shift;
my $prefix = shift;

unless ($infile) {
	&show_help;
	exit;
}

warn "\n\nreading input file...\n\n";


open(IN, $infile) or die;
while (<IN>) {
	chomp;
	my @fields = split /[\t]/, $_;
	my $pos = $fields[1];
	my $depth = $fields[3];
	my $read_bases = substr($fields[4],1);
	my @read_bases_array = split //, $read_bases;
	my $depth_fwd = ''; 
	my $depth_rev = '';
	foreach my $character (@read_bases_array) {
		if ($character =~ /[ACGT\,]/) {
			$depth_fwd++;
		}
		elsif ($character =~ /[acgt\.]/) {
			$depth_rev++;
		}
	}
	if (!$depth_rev) {
		$depth_rev = 0;
	}
	if (!$depth_fwd) {
		$depth_fwd = 0;
	}
	print "$depth_fwd $depth_rev\n";

}
close IN;


exit;
