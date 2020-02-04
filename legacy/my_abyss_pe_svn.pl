#!/software/perl-5.8.8/bin/perl -w

eval 'exec /software/perl-5.8.8/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;

=pod

=head1 my_abyss_pe_svn.pl

my_abyss_pe.pl -i <path_to/*.fast(a|q) -o <output_prefix> -k <kmer size>

=head1 SYNOPSIS

Usage: my_abyss_pe_svn.pl -i <path_to/*.fast(a|q) -o <output_prefix> -k <kmer size>

=head1 DESCRIPTION

This is a wrapper to run ABySS in parallel version on farm2 cluster.
Put your fasta or fastq files for each lane in a single directory.
They could be soft links to the files or copy them to the selected location.

The pipeline will bsub your jobs, so there's no need to bsub the script itself.
Also, it will choose the amount of CPUs and memory needed to run the assembly,
although I don't recommend you running more than 3lanes using this pipeline.

=head1 ARGUMENTS

<path_to/*.fast(a|q)>: Path to FASTA or FASTQ files. You can leave separated pairs 1 and 2 
<output_prefix>: Output prefix. Default: your_abyss
<kmer size>: Kmer size used for the assembly. It can't be longer than the actual read length
		
=head1 SEE ALSO

See the team wiki page http://scratchy.internal.sanger.ac.uk/wiki/index.php/Team_133

=head1 AUTHOR

Alejandro Sanchez, <as9@sanger.ac.uk>

=head1 TIME AND DATE

04-08-2009

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 Genome Research Limited. All Rights Reserved.              
                                                                               
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

#code goes here

use strict;
use Getopt::Std;

my $checkver = `ls -l /software/pathogen/external/applications/ABySS | grep ABySS`;
chomp $checkver;
my $version = (split /\s+/, $checkver)[$#_]; 

my %opts;
unless (@ARGV) {
        &USAGE;
}
getopts('hqvi:o:k:n:', \%opts);
&USAGE if $opts{h};

sub USAGE {

die 'Usage: my_abyss_pe.pl -i <path_to_fast(a|q) -o <output_prefix> -k <kmer size> [-nqv]

THIS SCRIPT WOULD BSUB YOUR JOB AUTOMATICALLY, DO NOT BSUB IT PLEASE!!!

	-i: directory to fasta or fastq files
	-o: output prefix (your_abyss by default)
	-k: kmer size
	-n: ABySS option; is the minimum number of pairs 
	    needed to consider joining two contigs.
	    Determined by trail. (10 by default)
	-q: trimming on. The program will trim by quality the end of the reads.
	-v: No run (verbose). Will print only the job call
	
'
}


### TEST NEW ABYSS##
#my $path = "/software/pathogen/external/applications/ABySS/1.1.0/bin";
#my $abyssp = "$path/ABYSS-P";
#my $pe = "$path/abyss-pe";

my $abyssp = "ABYSS-P";
my $pe = "abyss-pe";
my $mpi = "/nfs/users/nfs_a/as9/bin/mpiSub";
my $in = "";
my $out = "your_abyss";
my $k = "";
my $n = 10; #Default
my $c = 4; #CPUS

#print "in $opts{i} kmer $opts{k}"; <STDIN>;

if ($opts{i} and $opts{k}) {
	$in = $opts{i};
	$k = $opts{k};
	#print "in $in kmer $k"; <STDIN>;
} else {
	print STDERR "No input/kmer size specified...\n";
	&USAGE;
}

my $numfiles = `ls $in/ | grep -c fast`;
chomp $numfiles;
#print "$numfiles"; <STDIN>;

if ($opts{o}) {
	$out = $opts{o};
} else {
	print STDERR "No output prefix specified... Your result would be named 'your_abyss.*'\n";
}

if ($opts{n}) {
	$n = $opts{n};
}
my $abyssopt = "-k $k -g -s ";

if ($opts{q}) {
	$abyssopt .= "--standard-quality";
}

my $host = `hostname`;
chomp $host;
my $q = "parallel";
$c = $numfiles;

my $cpus = $c*4;
if ($cpus > 24) {
	$cpus = 24;
}

my $random = int(rand(100000));
my $call = "team133-bsub.pl -extra '-n $cpus' $q 15 $out.o $out.e $out.abyss.$random $mpi $abyssp $abyssopt \"$in/*.fast*\"";
my $call2 = "team133-bsub.pl -extra \"-w 'done($out.abyss.$random)'\" normal 15 $out-pe.o $out-pe.e $out\_pe '\"mv contigs.fa $out-3.fa && $pe k=$k j=4 n=$n in='$in/*.fast*' name=$out\"'";

if ($host =~ /pcs4/) {
	$q = "basement";
	$call = "team133-bsub.pl -extra '-n 2' $q 30 $out.o $out.e $out.abyss.$random $mpi $abyssp $abyssopt \"$in/*.fast*\"";
} elsif ($host =~ /seq1/) {
	$q = "phrapq";
	$call = "team133-bsub.pl -extra '-n 12' $q 30 $out.o $out.e $out.abyss.$random $mpi $abyssp $abyssopt \"$in/*.fast*\"";
	$call2 = "team133-bsub.pl -extra \"-w 'done($out.abyss.$random)'\" $q 15 $out-pe.o $out-pe.e $out\_pe '\"mv contigs.fa $out-3.fa && $pe k=$k j=4 n=$n in='$in/*.fast*' name=$out\"'";

}

if ($opts{v}) {
	print STDOUT "Using ABySS $version\n";
	print ("ABYSS-P job:\n$call\n\n");
	print ("abyss-pe job:\n$call2\n");
} else {
	print STDOUT "Using ABySS $version\n";
	system ("$call");
	system ("$call2");
}



