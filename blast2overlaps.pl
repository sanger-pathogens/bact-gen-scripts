#!/usr/bin/perl

# Martin Hunt, 26/09/2008

# parses blast output to find overlaps between pairs of reads, looking something
# like this:
#
# read1: --------+----------------+--
#                |     match      |
# read2:      ---+----------------+-----------------
#
# The overhang outside the match must be less than $dist_to_end, and the
# length of the match must be at least $min_length to count as an overlap.

use strict;
use warnings;
use List::Util qw(min max);

if ($#ARGV < 4){
    print "usage: blast2overlaps.pl <fasta file> <blast file> <minimum hit length> <max dist to ends> <min % identity>\n";
    exit;
}

my $path_to_bin = "~/bin/ReadSimulation";  # where this perl script is stored

my $fasta_file  = $ARGV[0];
my $blast_file  = $ARGV[1];
my $min_length  = $ARGV[2];  # minimum length to count as an overlap
my $dist_to_end = $ARGV[3]; # cutoff value of distance between start/end of
                            # read and start/end of overlap
my $min_pc      = $ARGV[4];

#my @overlap_lengths;   # to store the length of each overlap > $min_length
my %readlengths;  # to store read name->read length
#my $current_readname = "none";
#my @current_percents;
#my $per_read_output = "readname " . (join " ", @ranges) . "\n";
#my $freq_dist_output = "";



# put the fasta sequences into an array, one fasta file per element of array
open F, $fasta_file or die("Error opening $fasta_file\n");
my $whole_file = join ("", <F>);    
chomp $whole_file;
$whole_file = substr ($whole_file,1);  # removes the > at the start of the file
my @sequences = split("\n>", $whole_file);
close F;

chomp @sequences;

# fill hash of sequence name -> length
foreach (@sequences){
    my @fasta = split ("\n", $_);
    
    # if the fasta header has a space in it then only use the part before the space
    if ($fasta[0] =~ /(\S+)\s+(\S+)/){$fasta[0] = $1;}

    $readlengths{$fasta[0]} = length $fasta[1];
}




open F, $blast_file or die ("Error opening $blast_file");

while (<F>){
    if (/(^[^\#]\S+)\s\S+/){
#        if ($current_readname eq "none"){  # if 1st time round, set current_readname
#            $current_readname = $1;
#        }
#        else {
#            if ($current_readname ne $1) { # if next readname is not the same as the current readname
#                # reset the counts array
#                my @count = ();
#                for my $i (0 .. $#ranges){$count[$i] = 0;}
#
#                # count up the hits in each % range
#                foreach my $pc (@current_percents){
#                    for my $i (0 .. $#ranges){
#                        if ($pc >= $ranges[$#ranges - $i]) {
#                            $count[$#ranges - $i] += 1;
#                            last;
#                        }
#                    }
#                }              
#
#                $per_read_output .= "$current_readname @count\n";  # ...output results
#                $current_readname = $1;      # reset current readname
#                @current_percents = ();       # reset percentage array
#            }
#        }

        my $percent = isOverlap($_, $min_length, $dist_to_end, \%readlengths);
        if ($percent > $min_pc) {print $_}
    }
}

close F;







#______________________________________________________________________________

# takes line of output from blast, returns the %identity if it is an overlap of two reads,
# returns 0 otherwise
sub isOverlap {
    my $inputline       = shift;
    my $length_cutoff   = shift;
    my $end_range       = shift;
    my $length_hash_ref = shift;
#    my $overlap_ref     = shift;

    my @data = split(/\s+/, $inputline);   # split the input line into its columns

    # if the hit is a read to itself,then don't count it as an overlap
#    if ($data[3] < $length_cutoff || $data[0] eq $data[1]) {return 0;}
    if ($data[0] eq $data[1]) {return 0;}

    my $s1 = $data[6];
    my $e1 = $data[7];
    my $s2 = $data[8];
    my $e2 = $data[9];
  
 
    if ($s1 > $e1){
        $s1 = $length_hash_ref->{$data[0]} - $data[6] + 1;
        $e1 = $length_hash_ref->{$data[0]} - $data[7] + 1;
    }

    if ($s2 > $e2){
        $s2 = $length_hash_ref->{$data[1]} - $data[8] + 1;
        $e2 = $length_hash_ref->{$data[1]} - $data[9] + 1;
    }

    # if the start/end of an overlap is within $end_range of the end of the sequence,
    # then assume the overlap starts at the start/end
    if ($s1 < $end_range) {$s1 = 1}
    if ($s2 < $end_range) {$s2 = 1}
    if ($e1 > $length_hash_ref->{$data[0]} - $end_range) {$e1 = $length_hash_ref->{$data[0]}}
    if ($e2 > $length_hash_ref->{$data[1]} - $end_range) {$e2 = $length_hash_ref->{$data[1]}}

    # if the length of the overlap is < cutoff, don't count it
    if ($e1 - $s1 < $length_cutoff){return 0;}

    my $test_s1 = ($s1 == 1);
    my $test_e1 = ($e1 == $length_hash_ref->{$data[0]});
    my $test_s2 = ($s2 == 1);
    my $test_e2 = ($e2 == $length_hash_ref->{$data[1]});

    # if it is an overlap, return the %identity
    if ( ($test_s1 && $test_e1) || ($test_s1 && $test_e2) || ($test_s2 && $test_e1) || ($test_s2 && $test_e2)){
        my $a = $e1 - $s1 + 1;
 #       push @$overlap_ref, $a;
        return $data[2];
    }
    else { # otherwise, not an overlap, so return 0
        my $a = $e1 - $s1 + 1;
        return 0;
    }
}
