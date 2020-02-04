#!/usr/local/bin/perl -w

# change log:
# tdo 23.10. include 3 times view / crunch file, how all
# promer/nucmer hits / if comparison exist, do not redo / make it more
# modular.
# 

############## this version moves all subroutines to either contigO.pm and primerP.pm

# Copyright (C) 2008 Genome Research Limited. All Rights Reserved.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

=head1 Name

=head2 abacas.pl: Algorithmic Based Automatic Contiguation of Assembled Shotgun sequences

Copyright (C) 2008 The Wellcome Trust Sanger Institute, Cambridge, UK. This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public Licenseas published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

=head1 Synopsis

=head2 Usage:  abacas.pl -r <reference file: single fasta> -q <query sequence file: fasta> -p <nucmer/promer> [OPTIONS]

Assumes a UNIX like system

	-r	reference sequence should be a single fasta file
	-q	contigs in multi-fasta format
        -p	MUMmer program to use: 'nucmer' or 'promer'
	[OPTIONS]
        -h              print usage
	-d	0/1	use default nucmer/promer parameters [default 0]
	-s	int	minimum length of exact matching word (nucmer default = 12, promer default = 4)
	-m	0/1	print ordered contigs to file in multifasta format [default 0]
	-b	0/1	print contigs in bin to file [default 0]
	-N	0/1	print a pseudomolecule without "N"s [default 0]
	-i 	int 	mimimum percent identity [default 40]
	-v	int	mimimum contig coverage [default 40]
	-V	int	minimum contig coverage difference [default 1]
	-l	int	minimum contig length [default 1]
	-t	0/1	run tblastx on contigs that are not mapped [default 0]
	-g 	string (file name)	print gaps on reference to file name
	-a	0/1	append contigs in bin to the pseudomolecule
	-o	prefix  output files will have this prefix
	-P	0/1	pick primer sets to close gaps
	-f	int	number of flanking bases on either side of a gap for primer design (default 350)
    -R 0/1 default 1

=head1 Authors

Sammy Assefa (email: sa4@sanger.ac.uk), Thomas Keane (tk2) & Thomas Otto (tdo)

=head1 Overview

ABACAS is intended to rapidly contiguate (align, order, orientate) , visualize and design primers to close gaps on shotgun assembled contigs based on a reference sequence. It used MUMmer to find alignment positions and identify syntenies of assembly contigs against the reference. The output is then processed to generate a pseudomolecule taking overlaping contigs and gaps in to account. MUMmer's alignment generating programs, Nucmer and Promer are used followed by the 'delta-filter' utility function. Users could also run tblastx on contigs that are not used to generate the pseudomolecule. If the blast search results in mapping of extra contigs, finishers can use the visualization tool to easily modify ordering of contigs on the pseudomolecule. Gaps  in  the pseudololecule  are represented by "N"s. ABACAS could automatically extract gaps on the pseudomolecule and generate primer oligos for gap closure using Primer3. Uniqueness of primer sets is checked by running a sensitive NUCmer alignment.

=head1 Requirements

ABACAS requires MUMmer to be installed in the working directory for ordering, orienting of contigs. The Arthemis Comparision Tool (ACT) should be downloaded for visualizing scaffolding of contigs. Primer design part of the programme requires  Primer3. Optinally, BLASTALL is required in order to run tblastx on the contigs  that are not mapped using  Nucmer or Promer.

=head1 Input:

Two fasta files containing the reference and query (contigs) sequences are required. The reference file should be in a single fasta format for speedy contig ordering and orientation.

=head1 Outout:

=head2 Running the script with default options will generate the following files:

    1.  ordered and orientated sequence file (reference_query.fasta) 
    2.  a feature file (reference_query.tab) 
    3.  a bin file thatcontains contigs that are not used (reference_query.bin) 
    4.  a comparision file (reference_query.crunch) 
    5.  a file with gap information (reference_query.gaps) 
    6.  a file that contains information on contigs that have a mapping information but could not be used in the ordering (unused_contigs.out)
    7.  a feature file to view contigs with ambiguous mapping (reference.notMapped.contigs.tab). This file should be uploaded on the reference side of ACT view.
    8.  a file that shows how repetitive the reference genome is (reference.Repeats.plot).
    Files 7 & 8 should be  uploaded on the reference side of ACT view. 


Please note that contigs in the '.fasta' file will be reverse complemented if they are found to map on the reverse strand. However, the ACT view shows the initial orientation of these contigs i.e. will be shown on the reverse strand. If you write a fasta file of the pseudomolecule from ACT, the resulting sequence will be a set of ordered contigs (the orientation will not change). It is therefore recommended to use the '.fasta' pseudomolecule file automatically generated for further investigation.

=head2 Optional Output files:

It is also possible to generate other files including:
    1.  A list of ordered and orientated contigs in a multi-fasta format (-m 1)
    2.A pseudomolecule with all unmapepd contigs appended to the end for reordering (-a 1) 
    3.A pseudomolecule where the gaps are not padded with N (-N 1) 
    4.A multi-fasta file of all unmapped contigs (-b 1) 
    5.A multi-fasta file of regions on the reference that correspond to gaps on pseudomolecule (-g file_name) 
    6.A list of sense and antisense primer sets in separate files 
    7.A list of locations where sense and antisense primers are found in two separate files 
    8.A standard primer3 output summary file with a detailed information on oligos.

=head1 Options:

-d 	 default 0 i.e. increase mapping sensitivity by using 'all anchor matches regardless of their uniqueness during mapping' ( i.e. --maxmatch) 

       -d  1 runs default NUCmer or PROmer parameters

-m     default -m 0 
  
       *this option is helpful if the user wants to further investigate the ordering using other alignment algorithms such as blast.
 
       -m 1 to print ordered and orientated contigs to file

-b    default -b 0
  
       *contigs that are not used in generating the pseudomolecule will be placed  in a '.bin' file. This file only contains contig names.
  
      -b 1 will print those contigs in the bin to multi-fasta file. This option must be used with -t 1 in order to run 'tblastx' on contigs in bin

-N    default -N 0
 
      *ABACAS produces a pseudomolecule ('.fasta' file) and fills gaps with 'N's. This option will remove "N"s.

       -N 1 will generate another pseudomolecule without 'N's

-i    default 40
  
       *minimum percent identity could vary from 0 to 100 depending on the closeness of the two genomes. Choosing a smaller value will pull in more       	contigs and vice versa

-v    default 40
       
       *minimum contig coverage: set a value between 0 and 100
 
-V    default 1
	
        *minimum contig coverage difference. Use -V 0 to place contigs randomly to one of the positions (in cases where a contig maps to multiple places)

-l     default 100
	 * contigs below this cutoff will not be used

-t     default 0
	 *runs tblastx on contigs that are not used to generate the pseudomolecule.
	 -t 1 will run blastall on contigs in the .bin file

-g     file_name 

       * will print sequences of the reference that correspond with gaps on the pseudomolecule in a multi-fasta format

-a     default -a 0
   
       * this option will append contigs that re not used to the end of the pseudomolecule. Contigs will then be easily manipulated and re-ordered using 	ACT's graphical interface

-o     prefix (string)
 
       *output files will have this prefix

-P    default 0
 
      * -P 1 will pick primer oligos to close gaps

-f	 default 350 

      *number of flanking bases on either side of a gap for primer design (default 350)

=head1 Colour code:

The feature (file 2 from default output section) file has the following colour codes: 
    Dark blue (4):     contigs with forward orientation
    Dark green (3):    contigs with reverse orientations 
    Sky blue (5):      contigs that overlap with the next contig
    Yellow (7):        contigs that have no hit (only added to the pseudomolecule if '-a 1' is used)

=cut
# load modules
use strict;
use warnings;
use POSIX qw(ceil floor);
use Getopt::Std;
use Getopt::Long;
use Data::Dumper;
use lib "/nfs/team81/sa4/abacs/unstable/";
use contigO;
use primerP;

my %options;
getopts ('r:q:p:h:d:s:m:b:N:i:v:V:l:c:t:g:a:o:P:f:R:u:e:k:z', \%options);

#my $ts = ABACAS::test;

if (defined $options{h} && $options{h} eq ""){ die contigO::usage();}
unless (defined $options{r} && defined $options{q}) {die contigO::usage();}
my $trac =0;
my $escapeToPrimers = 0;

if (defined $options{e})
{
    contigO::print_header(); 
    print "Primer design selected,... escaping contig ordering\n";
    $escapeToPrimers =$options{e};
    $options{p} = "nucmer";
    #unless (defined $options{r} && defined $options{q}) {die contigO::usage();}
    #primerP::pickPrimers2 ();
}
else
{
    unless (defined $options{r} && defined $options{q} && defined $options{p} ){ die contigO::usage();}
}
unless (defined $options {d}){$options{d} = 0;}
unless (defined $options{s}) {$options{s} = 12; $trac =1;}
unless (defined  $options{m}){$options{m} =0;}
unless (defined $options {b}){$options{b} =0;}
unless (defined $options {N}){$options{N} = 0;}
unless (defined $options{t}){$options{t} =0;}
unless (defined $options{i}){$options{i} =40;}
unless (defined $options{v}){$options{v} = 40;}
unless (defined $options {V}){$options{V} = 1;}
unless (defined $options {l}){$options{l} = 1;}
unless (defined $options {a}){$options {a} =0;}
unless (defined $options {P}){$options {P} =0;}
unless (defined $options {f}){$options {f} =350;}
unless (defined $options {u}){$options {u} ="nucmer";}
unless (defined $options {R}){$options {R} ="1";}
unless (defined $options {k}){$options {k} = 0;}
#assign options to variables
my ($reference, $query_file, $choice, $sen, $seed, $mlfas, $fasta_bin, $avoid_Ns, $tbx, $min_id, $min_cov, $diff_cov, $min_len,$add_bin_2ps, $pick_primer, $flank, $chk_uniq,$redoMummer, $is_circular)=
($options{r}, $options{q}, "$options{p}","$options{d}",$options{s},$options{m}, $options{b}, $options{N}, $options{t}, $options{i}, $options{v}, $options{V}, $options{l}, $options{a}, $options{P}, $options{f}, $options{u}, $options{R}, $options{k});
#check for user inputs to run tblastx
if ($tbx ==1 && $fasta_bin !=1)
{
	print "ERROR:  Please use  -t 1 -b 1  if you want to run tblastx on contigs in bin\n";
	exit;
}

if ($escapeToPrimers ==1)
{
    #ontigO::print_header();
    primerP::pickPrimers ($reference, $query_file, $flank, $chk_uniq);
    exit;
}
#BEGIN
contigO::print_header(); 
#define vars. for summary stats file
my $num_contigs =0; my $num_inbincontigs=0;
my $avg_cont_size =0; my $num_overlaps =0; my $num_gaps =0; my $num_mapped =0;
my $total_bases_mpd =0; my $p_ref_covered =0; my $num_ambigus =0;
my $ref_bp=10000; my $num_inComparisoncontigs=0; my $num_fortillingcontigs=0;
my $num_notsetTilling=0;
my $debug=1;

#define files to store values for options g (name of gap file), o (prefix)
my ($gaps_2file, $ref_inline, $prefix);
#open file to write gaps regions on reference to file
if (defined $options {g})
{
    $gaps_2file = "$options{g}";
    open (refgapsFH, '>', $gaps_2file) or die "Could not open file $gaps_2file for write:$!\n";
    $ref_inline = contigO::Ref_Inline($reference);
    #print $ref_inline; ## delete    
}

#Get length of the reference sequence

#This lines are used to get the length of the reference sequence. 
my $lin = contigO::Ref_Inline($reference);
my $ref_len =  length($lin);

if (($choice eq "nucmer")||($choice eq "promer"))
{
	print "PREPARING DATA FOR   $choice \n";
}
else
{
	die ("$choice is not a recognised MUMmer function\n");
}
if ($trac ==1 && $choice eq "promer")
{
    $seed =4; 
}
if ($choice eq "promer" && $seed >6)
{
	print "WARNING: minimum value of exact words is greater than the default.\n ";
}
my $df = 'delta-filter';
my $st = 'show-tiling';
my $ask = 'which';
my ($path_dir, $run_mum);
my $path_toPass;
###################
# Running MUMmer   #
###################

if ($debug)
{
    print "the seed is  $seed \n";
    print "RedoMummer= ",$redoMummer."\n";
}



my @do_mum_return;
if ($redoMummer)
{
    @do_mum_return = contigO::doMummer($reference, $query_file, $choice,$sen,$seed,$min_id, $min_cov, $diff_cov,$min_len, $debug, $is_circular);
}

####################################
# Processing tiling output         #
####################################
if ($debug) {
  print "Do tiling...\n";
}

#### ???????????????????? progress ????????
#doTilling();
###########################################

#get nucmer or promer tiling file into variable 

my $mummer_tiling = $do_mum_return[0];
$path_dir = $do_mum_return[1];
$path_toPass = $do_mum_return[2];

if ($debug ){
    print "Tiling file is in ", $mummer_tiling, "\n";
    print "Hash the contigs...\n";
    }
#get nucmer or promer tiling file into variable 

my ($href, $ref_contigs_qual) = contigO::hash_contigs2($query_file);
my $qualAvailable=0;		## flag to know, that the qual exists
if ($debug) {
  print "Hash_ed the contigs...\n";
}


my %contigs_hash = %{$href};
my @c_keys = keys %contigs_hash;
$num_contigs = scalar(@c_keys);
my @cont_lens;

#define variables
my (@ids ,$id, $id_len);
my (@Rs, @Re, @G, @Len, @cov, @pid, @orient, @Cid);
my (@Ps, @Pe);

######## To tdo...... 
if ($debug) {
  print "Get comparison\n";
}

my $file_to_comp = $choice.".filtered.delta";
my ($ref_coords,$ref_rank)= contigO::getMummerComparison($file_to_comp);#### move 

if ($debug) {
  print "Got comparison\n load the tilling graph\n";
}

open (TIL, $mummer_tiling) or die "Could not open $mummer_tiling: $!";
while (<TIL>)
{
	chomp;
	if ($_ =~/^>/)
	{
		my $line = substr $_, 1;
		my @splits = split /\s+/, $line;
		$id = $splits[0];
		push @ids, $id;
		$id_len= $splits[1];
	}
	else 
	{
		my @splits = split /\s+/, $_;
		push @Rs, $splits[0];
		push @Re, $splits[1];
		push @G, $splits[2];
		push @Len, $splits[3];
		push @cov, $splits[4];
		push @pid, $splits[5];
		push @orient, $splits[6];
		push @Cid, $splits[7];
	}
	
}
close (TIL);
my $total;	
if (scalar(@Rs) != scalar(@Re))
{
	print "ERROR: unequal array size\n";
	exit;
}
else
{
	$total = scalar(@Rs);
	$num_mapped = scalar(@Rs);
	
}
my $ref_loc = $reference;  # get locations of reference and query files
my $qry_loc = $query_file;
#define file handles for output files
my ($seqFH,$tabFH,$binFH,$crunchFH, $gapFH, $mlFH, $dbinFH, $avoidNFH);
my $dif_dir =0; 	#assume query and reference are in the working directory
my ($new_query_file, $new_reference);

if ($debug) {
  print "done tilling\n";
  
}

#tdo

########################################
# map unused contigs to refenrece      #
########################################

if (1) {
  
  my @splits_ref = split (/\//, $reference);
  my $new_reference_file = $splits_ref[(scalar(@splits_ref)-1)];
  
  $num_notsetTilling=contigO::writeBinContigs2Ref("unused_contigs.out",$new_reference_file);
  
}


########################################
# find repeats in the reference seq    # 
########################################

#findRepeats($reference,$new_reference_file);

#tdo end



########################################
# Opening files to write final output  #
#######################################

if (defined $options {o}) #check is a prefix is provided by user
{
	$prefix  = "$options{o}";
	open ($seqFH, '>', $prefix . '.fasta') or die "Could not open file $prefix.fasta for write: $!\n";
	open ($tabFH, '>', $prefix . '.tab') or die "Could not open file $prefix.tab for write: $!\n";
	open ($binFH, '>', $prefix . '.bin') or die "Could not open file  $prefix.bin for write: $!\n";
	open ($crunchFH, '>', $prefix . '.crunch') or die "Could not open file $prefix.crunch for write: $!\n";
	open ($gapFH, '>', $prefix . '.gaps') or die "Could not open file $prefix.gaps for write: $!\n";
	if ($mlfas ==1)
	{
		open ($mlFH, '>', $prefix . '.contigs.fas') or die "Could not open file $prefix.contigs.fas for write: $!\n";	
	}
	if ($fasta_bin ==1)
	{
		open ($dbinFH, '>', $prefix . '.contigsInbin.fas') or die "Could not open file $prefix.contigsInbin.fas for write: $!\n";
	}
	if ($avoid_Ns ==1)
	{
		open ($avoidNFH, '>', $prefix .'.NoNs.fasta') or die "Could not open file $prefix.NoNs.fasta for write: $!\n";
	}
}
else
{
	#check if the query and ref are in the working directory and open files for writing 
	if ($query_file =~ /^\../ || $reference =~ /^\../ || $query_file =~ /^\// || $reference =~ /^\// ) 
	{
		$dif_dir =1; 
		my @splits_q = split (/\//, $query_file);
		my @splits_r = split (/\//, $reference);	#split the file name by '/'
		my $num_q = scalar (@splits_q);			#counts the number of directories to the file
		my $num_r = scalar(@splits_r);
		$new_query_file = $splits_q[$num_q -1];
		$new_reference=$splits_r[$num_r -1];
		open ($seqFH, '>', $new_query_file . "_" . $new_reference . '.fasta') or die "Could not open file $query_file\_$reference.fasta for write: $!\n";
		open ($tabFH, '>', $new_query_file . "_" . $new_reference . '.tab') or die "Could not open file $query_file\_$reference.tab for write: $!\n";
		open ($binFH, '>', $new_query_file . "_" . $new_reference . '.bin') or die "Could not open file  $query_file\_$reference.bin for write: $!\n";
		open ($crunchFH, '>', $new_query_file . "_" . $new_reference . '.crunch') or die "Could not open file $query_file\_$reference.crunch for write: $!\n";
		open ($gapFH, '>', $new_query_file . "_" . $new_reference . '.gaps') or die "Could not open file $query_file\_$reference.gaps for write: $!\n";
		if ($mlfas ==1)
		{
			open ($mlFH, '>', $new_query_file . "_" . $new_reference . '.contigs.fas') or die "Could not open file $query_file\_$reference.contigs.fas for write: $!\n";
		}
		if ($fasta_bin ==1)
		{
			open ($dbinFH, '>', $new_query_file . "_" . $new_reference . '.contigsInbin.fas') or die "Could not open file $query_file\_$reference.contigsInbin.fas for write: $!\n";
		}
		if ($avoid_Ns ==1)
		{
			open ($avoidNFH, '>', $new_query_file . "_" . $new_reference . '.NoNs.fasa') or die "Could not open file $query_file\_$reference.NoNs.fasta for write: $!\n";
		}
	}
	else
	{
		open ($seqFH, '>', $query_file . "_" . $reference . '.fasta') or die "Could not open file $query_file\_$reference.fasta for write: $!\n";
		open ($tabFH, '>', $query_file . "_" . $reference . '.tab') or die "Could not open file $query_file\_$reference.tab for write: $!\n";
		open ($binFH, '>', $query_file . "_" . $reference . '.bin') or die "Could not open file  $query_file\_$reference.bin for write: $!\n";
		open ($crunchFH, '>', $query_file . "_" . $reference . '.crunch') or die "Could not open file $query_file\_$reference.crunch for write: $!\n";
		open ($gapFH, '>', $query_file . "_" . $reference . '.gaps') or die "Could not open file $query_file\_$reference.gaps for write: $!\n";
		if ($mlfas ==1)
		{
			open ($mlFH, '>', $query_file . "_" . $reference . '.contigs.fas') or die "Could not open file $query_file\_$reference.contigs.fas for write: $!\n";	
		}
		if ($fasta_bin ==1)
		{
			open ($dbinFH, '>', $query_file . "_" . $reference . '.contigsInbin.fas') or die "Could not open file $query_file\_$reference.contigsInbin.fas for write: $!\n";
		}
		if ($avoid_Ns ==1)
		{
			open ($avoidNFH, '>', $query_file . "_" . $reference . '.NoNs.fasta') or die "Could not open file $query_file\_$reference.NoNs.fasta for write: $!\n";
		}
	}
}

### new crunch file...
my $resCrunch;

#define start of a pseudomolecule
my $ps_start =1;
$Ps[0] = 1; 
$Pe[0] = $Ps[0] + $Len[0] -1;

print $tabFH "ID   ",$id, "\n";
print $seqFH ">", "ordered_", $id, "\n";

my $tmp_qual;
my $tmp_nqual;
my $tmp_seq ="";
my $tmp_nseq ="";
#$num_mapped =0; my

print "total $total \n";


for (my $i=1; $i <= $total; $i+=1)
{
        my $g;
        if ($G[$i -1] <= 0)
	{
            $g =1;
	  
            if (defined($Len[$i]))
	    {
		  $Ps[$i] =  $Pe[$i-1] +$g;
		  $Pe[$i] = $Ps[$i] + $Len[$i] -1;
		  $total_bases_mpd+=$Len[$i];
	    }
	}
        else
	{
	  $g =$G[$i -1]; 
	  
	  if (defined($Len[$i]))
		{
		  $Ps[$i] =  $Pe[$i-1] +$g +1;
		  $Pe[$i] = $Ps[$i] + $Len[$i] -1;
		  $total_bases_mpd+=$Len[$i];
		}
	}
        if ($Rs[$i -1] <0) #check if a reference starting position is less than 0
	{
            $Rs[$i -1] =1;
	}
        my $covv =sprintf("%.0f",$cov[$i -1]); #ROUNDING
        my $pidd = sprintf("%.0f", $pid[$i -1]);
  
        if ($mlfas ==1)
	{
            print $mlFH ">",$Cid[$i -1], "\n";
	}
        my $cont_cord;
        my $col;
        if($orient[$i-1] eq "+")
        {
            $cont_cord = $Ps[$i -1]."..".$Pe[$i-1];
            $col = 4;
            my $c_seq = $contigs_hash{$Cid[$i-1]};
            push (@cont_lens, length($c_seq));
            # tdo
            my $c_qual='';
	  
            if (defined($$ref_contigs_qual{$Cid[$i-1]})) {
                  ## flag to know, that the qual exists
                  $qualAvailable=1;
                  
                  $c_qual = $$ref_contigs_qual{$Cid[$i-1]};
            }
            $tmp_qual .= $c_qual;
  
            $tmp_seq = $tmp_seq.$c_seq;
	  
            if ($avoid_Ns ==1)
              {
                $tmp_nseq = $tmp_nseq.$c_seq;
                #tdo
                $tmp_nqual .= $c_qual;
              }
            if ($mlfas ==1)
            {
                my $len_mlf = length($c_seq);
                if ($len_mlf <= 60)
                      {
                        print $mlFH $c_seq, "\n"; 
                      }
                elsif ($len_mlf> 60 )
                      {
                        for (my $i =0; $i < $len_mlf; $i+=60)
                          {
                                my $tmp_s = substr $c_seq, $i, 60;
                                print $mlFH $tmp_s, "\n"; 
                          }
                      } 
            }
	}
        else
	{
		$cont_cord = "complement(".$Ps[$i -1]."..".$Pe[$i-1].")";
		$col =3;
		my $c_seq = contigO::revComp($contigs_hash{$Cid[$i-1]}); #REVERSE COMPLEMENT A SEQUENCE
                push @cont_lens, length($c_seq);
		$tmp_seq = $tmp_seq.$c_seq;

		# tdo
		my $c_qual='';
                if (defined($$ref_contigs_qual{$Cid[$i-1]})) {
		  $c_qual = contigO::reverseQual($$ref_contigs_qual{$Cid[$i-1]});
		}
		$tmp_qual .= $c_qual;
		
		if ($avoid_Ns ==1)
		{
			$tmp_nseq = $tmp_nseq.$c_seq;
			# tdo
			$tmp_nqual .= $c_qual;
			
		}
		if ($mlfas ==1)
		{
			#print $mlFH $c_seq, "\n";
			my $len_mlf = length($c_seq);
			if ($len_mlf <= 60)
			{
			    print $mlFH $c_seq, "\n"; 
			}
			elsif ($len_mlf> 60 )
			{
			    for (my $i =0; $i < $len_mlf; $i+=60)
			    {
				my $tmp_s = substr $c_seq, $i, 60;
				print $mlFH $tmp_s, "\n"; 
			    }
			} 	
		}
       	}
	if ($Re[$i -1] > $ref_len)
	{
		$Re[$i -1] = $ref_len -1;
	}
	if ($Pe[$i -1] > length($tmp_seq))
	{
		$Pe[$i -1] = length($tmp_seq);
	}


#### sa4 commented the following two lines...
  #print "here $i \n";
  
  #$resCrunch .= getPosCoords($ref_coords,$Cid[$i-1], $Ps[$i -1]);
  
  
	print $crunchFH $covv, " ", $pidd, " ", $Ps[$i -1], " ", $Pe[$i -1], " ", $Cid[$i -1], " ", $Rs[$i -1], " ", $Re[$i-1], " ", "unknown NONE\n";
	
        #WRITE FEATURE FILE
	print $tabFH "FT   contig          ",$cont_cord, "\n";
	print $tabFH "FT                   ", "/systematic_id=\"", $Cid[$i-1],"\"","\n";
	print $tabFH "FT                   ", "/method=\"", "mummer\"", "\n";
        print $tabFH "FT                   ", "/Contig_coverage=\"",$cov[$i -1], "\"", "\n";
        print $tabFH "FT                   ", "/Percent_identity=\"",$pid[$i -1], "\"", "\n";
               
	
        if ($g>1)       #WRITE GAP LOCATIONS AND SIZE TO FILE
	{
		my $gs = $Pe[$i -1] +1;
		my $ge = "";
		if (defined $Ps[$i])
		{
			$ge = $Ps[$i] -1;
		}
                my $rs = $Re[$i -1] +1;
                my $re;
                if (defined $Rs[$i])
                {
                         $re =$Rs[$i]-1;
                }
                else
                {
                        $re = "END";
                }
		print $gapFH "Gap\t",$g, "\t", $gs, "\t", $ge, "\t", $rs, "\t", $re,"\n";
		print $tabFH "FT                   ", "/colour=\"",$col, "\"", "\n";
		if (defined $options{g} && $rs < $ref_len)
		{
			my $gapOnref = substr ($ref_inline, $rs, $g);
			print refgapsFH ">",$g,"_",$rs, "\n";
			my $lenG = length($gapOnref);
			if ($lenG <= 60)
			{
			    print refgapsFH $gapOnref, "\n"; 
			}
			elsif ($lenG> 60 )
			{
			    for (my $i =0; $i < $lenG; $i+=60)
			    {
				my $tmp_s = substr $gapOnref, $i, 60;
				print refgapsFH $tmp_s, "\n"; 
			    }
			}
			
		}
		
		my $ns = contigO::makeN($g);
		$tmp_seq = $tmp_seq.$ns;
		#tdo
		for (1..length($ns)) {
		  $tmp_qual .= "0 ";
		}
	  }
	else
	{
		$col = 5;
		print $tabFH "FT                   ", "/Overlapping=\"", "YES\"", "\n";
		print $tabFH "FT                   ", "/colour=\"",$col, "\"", "\n";
			
	}
        
        #print "done...", $i, "\t", $Cid[$i -1], "\n";
}

#causing errors...
        #open (F,"> test.crunch2") or die "an";
        #print F $resCrunch;
        #close(F);


#tdo
my @Quality_Array;
if ($qualAvailable) {
  @Quality_Array = split(/\s/,$tmp_qual);
  my $res;
  foreach (@Quality_Array) {
	$res .= "$_\n";
	
  }
  ## get name
  my @splits_query = split (/\//, $query_file);
  $new_query_file = $splits_query[(scalar(@splits_query)-1)];
  open (F,"> $new_query_file.qual.plot") or die "problems\n";
  print F $res;
  close(F);
}

##WRITE PSEUDOMOLECULE WITHOUT 'N's
if ($avoid_Ns ==1)
{
	print $avoidNFH ">", "ordered_", $id, "without 'N's","\n";
	my $len2 = length($tmp_nseq);
	if ($len2<=60)
	{
		print $avoidNFH $tmp_nseq, "\n"; 
	}
	elsif ($len2>60)
	{
		for (my $i =0; $i < $len2; $i+=60)
		{
			my $tmp_s = substr $tmp_nseq, $i, 60;
			print $avoidNFH $tmp_s, "\n"; 
		}
	}
}

####################################
#WRITE CONTIGS WITH NO HIT TO FILE #
#################################
my %Cids;

foreach(@Cid)
{
	chomp;
	$Cids{$_} = 1;
}
my @contigs_2bin = ();
my %h_contigs_2bin;

foreach (@c_keys)
{
    push(@contigs_2bin, $_) unless exists $Cids{$_};
	
}
foreach(@contigs_2bin)
{
  $h_contigs_2bin{$_}=1;
  
  print $binFH "$_ \n";
	
}
$num_inbincontigs= scalar(@contigs_2bin);

###############################
#WRITE PSEUDOMOLECULE TO FILE #
###############################
if (1) #appending unmapped contigs to pseudomolecule
{
  my $newseq=">rubishishish\n";
  my $crunch="";
  my $length=1;
  my $tab="";
  my $empty_bin=0;
  
    #print    $contigs_hash{$contigs_2bin[0]};
  foreach my $pos (sort {$a <=> $b} keys %$ref_rank)
  {
	my $contig =$$ref_rank{$pos};
	if ($debug ) {
	  print "contig $contig\n";
	  
	}
	if (defined($h_contigs_2bin{$contig})) {
	  
	  #  for (my $i =0; $i < scalar(@contigs_2bin); $i+=1)
	  #	{
	  my $binseq = $contigs_hash{$contig};
	  #	  print length ($binseq);
	  if ($binseq ne "") {$empty_bin +=1;}
	  $newseq .= $binseq."\n";
	  $crunch .= contigO::getPosCoordsTurn($ref_coords,$contig, $length);
	  
	  $tab .= "FT   contig          $length..".(length($binseq)+$length-1)."\n";
	  $tab .= "FT                   /systematic_id=\"". $contig."\"\n";
	  $tab .= "FT                   /method=\"mummer\"\n";
	  $tab .= "FT                   /colour=\"10\"\n";
	  #		$tab .= "FT                   ", "/", $note, "\n";
	  $length += length ($binseq);
	  
	}
        
        
 }
  ### print
  
  if ($empty_bin != 0)
  {
    open (F,"> head1.fasta") or die "p \n";
    print F $newseq;
    close(F);
  }
  
if ($crunch ne "")
  {
  open (F,"> head1.crunch") or die "p \n";
  print F $crunch;
  close(F);
  }
  
  if ($tab ne "")
  {
    open (F,"> head1.tab") or die "p \n";
    print F $tab;
    close(F);
  }
}


if ($add_bin_2ps ==1) #appending unmapped contigs to pseudomolecule
{
    my $new_seq = $tmp_seq;
    #print "length before bin\t", $new_seq, "\n";
    
    if ($debug) {print scalar(@contigs_2bin), " contigs in bin\n";}
    for (my $i =0; $i < scalar(@contigs_2bin); $i+=1)
    {
            my $binseq = $contigs_hash{$contigs_2bin[$i]};
            $new_seq = $new_seq.$binseq;
    }
    my $prev_len = length($tmp_seq);
    print "Prev_len..", $prev_len, "\n";
    my $len=length($new_seq);
    print "New Len= ", $len, "\n";
    for (my $i =0; $i < scalar(@contigs_2bin); $i+=1)
    {
            my $binseq = $contigs_hash{$contigs_2bin[$i]};
            my $len_current_contig = length($binseq);
            my $start = $prev_len +1;
            my $end = $start + $len_current_contig -1; 
            my $col = 7;
            if ($start > $len)
            {
                    $start = $len;
            }
            if ($end >$len)
            {
                    $end = $len;
            }
            my $co_cord = $start."..".$end;
            my $note = "NO_HIT";
            print $tabFH "FT   contig          ",$co_cord, "\n";
            print $tabFH "FT                   ", "/systematic_id=\"", $contigs_2bin[$i],"\"","\n";
            print $tabFH "FT                   ", "/method=\"", "mummer\"", "\n";
            print $tabFH "FT                   ", "/colour=\"",$col, "\"", "\n";
            print $tabFH "FT                   ", "/", $note, "\n";
            $prev_len= $end;
    }
    my $num_lines = ceil ($len/60);
    if ($len <= 60)
    {
        print $seqFH $new_seq, "\n"; 
    }
    elsif ($len> 60 )
    {
        for (my $i =0; $i < $len; $i+=60)
        {
            my $tmp_s = substr $new_seq, $i, 60;
            print $seqFH $tmp_s, "\n"; 
        }
    }
}
else
{
    my $len=length($tmp_seq);  #print "len ps= ", $len, "\n";
    my $num_lines = ceil ($len/60);
    if ($len <= 60)
    {
        print $seqFH $tmp_seq, "\n"; 
    }
    elsif ($len> 60 )
    {
        for (my $i =0; $i < $len; $i+=60)
        {
            my $tmp_s = substr $tmp_seq, $i, 60;
            print $seqFH $tmp_s, "\n"; 
        }
    }
}
##################################
#WRITE	CONTIGS IN BIN TO FILE   #
#################################
if ($fasta_bin ==1)
{
    foreach(@contigs_2bin)
    {
        print $dbinFH ">", $_, "\n";
        my $s = $contigs_hash{$_};
        my $len3 = length($s);
        if ($len3<=60)
        {
            print $dbinFH $s, "\n"; 
        }
        elsif ($len3>60)
        {
            for (my $i =0; $i < $len3; $i+=60)
            {
                    my $tmp_s = substr $s, $i, 60;
                    print $dbinFH $tmp_s, "\n"; 
            }
        }
    }
}

###################
#Final message    #
##################
if (0)
{
unlink ("$choice.delta");
unlink ("$choice.filtered.delta");
unlink ("$choice.cluster");
unlink ("$choice.tiling");
}

my $sumFH;
print " FINISHED CONTIG ORDERING\n";

############## commented by sa4
# print summary stats:
#contigO::printStats($num_fortillingcontigs, $num_notsetTilling, $num_mapped, $num_contigs,$num_inComparisoncontigs, $ref_len, $total_bases_mpd); ### sub taken to contigO

###################################
# Running blast & picking primers #
##################################

if (defined $options {o})
{
	$prefix  = "$options{o}";
		print STDERR "\nTo view your results in ACT\n\t\t Sequence file 1: $reference\n\t\t Comparison file 1: $prefix.crunch\n\t\t Sequence file 2: $prefix.fasta\n
\t\tACT feature file is: $prefix.tab\n
\t\tContigs bin file is: $prefix.bin\n
\t\tGaps in pseudomolecule are in: $prefix.gaps\n\n";
	if ($tbx ==1)
	  {
	  print "Running tblastx on contigs in bin...\nThis may take sevral minutes ...\n";
		my $formatdb = 'formatdb -p F -i';
		#		my @formating = `
		#		system("$formatdb $reference"); or die "ERROR: Could not find 'formatdb' for blast\n";
		my $blast_opt = 'blastall  -m 9 -p tblastx -d ';
		my $contigs_inBin = $prefix.'.contigsInbin.fas';
		
		my $call = "$blast_opt $reference -i $contigs_inBin -o blast.out";
	#	!system($call) or die "ERROR: Could not find 'blastall' , please install blast in your working path (in picking primers)\n\n $blast_opt $reference -i $contigs_inBin -o blast.out\n" ;
	}
	if ($pick_primer == 1)
	{
		print " DESIGNING PRIMERS FOR GAP CLOSURE...\n";
		my $qq = "$prefix.fasta";
                if ($debug) {print "Passing inputs ...\n";}
		primerP::pickPrimers($qq, $reference, $flank, $path_toPass, $chk_uniq,$qualAvailable,@Quality_Array);
	}
}
else
{
	if ($dif_dir ==0)
	{
		print STDERR "To view your results in act use:   act $reference $query_file\_$reference.crunch $query_file\_$reference.fasta\n
		ACT feature file is: $query_file\_$reference.tab\n
		Contigs not used in creating a pseudomolecult : $query_file\_$reference.bin\n
		Gaps in pseudomolecule are in: $query_file\_$reference.gaps\n\n";
		if ($tbx ==1)
		{
			print "Running tblastx on contigs in bin...\nThis may take sevral minutes ...\n";
			my $formatdb = 'formatdb -p F -i' ;
	#		my @formating = `
			!system("$formatdb $reference") or die "ERROR: Could not find 'formatdb' for blast\n";
			my $blast_opt = 'blastall  -m 9 -p tblastx -d ';
			my $contigs_inBin = $query_file . "_" . $reference . '.contigsInbin.fas';
#			my @bigger_b = `
			!system("$blast_opt $reference -i $contigs_inBin -o blast.out") or die "ERROR: Could not find 'blastall' , please install blast in your working path (other dir==0)\n$blast_opt $reference -i $contigs_inBin -o blast.out\n \n";
		}
		if ($pick_primer == 1)
		{
			print " DESIGNING PRIMERS FOR GAP CLOSURE...\n";
			my $qq = "$query_file\_$reference.fasta";
                        if ($debug) {print "Passing inputs ...\n";}
                        print $qq, " ",$reference, " ",$flank, " ",$path_toPass, " ",$chk_uniq, " ", $qualAvailable, "\n";
                        print " qual arr: ", @Quality_Array;
			primerP::pickPrimers($qq, $reference, $flank, $path_toPass, $chk_uniq,$qualAvailable,@Quality_Array);
                        
		}
	}
	else
	{
		
		print STDERR "To view your results in ACT use:   act $new_reference $new_query_file\_$new_reference.crunch $new_query_file\_$new_reference.fasta\n
		ACT feature file is: $new_query_file\_$new_reference.tab\n
		Contig bin file is: $new_query_file\_$new_reference.bin\n
		Gaps in pseudomolecule: $new_query_file\_$new_reference.gaps\n\n";
		if ($tbx ==1)
		{
			print "Running tblastx on contigs in bin...\nThis may take sevral minutes ...\n";
			my $formatdb = 'formatdb -p F -i';
#			my @formating = `
			!system("$formatdb $new_reference") or die "ERROR: Could not find 'formatdb' for blast\n";
			my $blast_opt = 'blastall  -m 9 -p tblastx -d ';
			my $contigs_inBin = $new_query_file . "_" . $new_reference . '.contigsInbin.fas';
			#			my @bigger_b = `
			!system("$blast_opt $reference -i $contigs_inBin -o blast.out") or die "ERROR: Could not find 'blastall' , please install blast in your working path (else case)\n\n$blast_opt $reference -i $contigs_inBin -o blast.out\n";
		  }
		if ($pick_primer == 1)
		  {
			print " DESIGNING PRIMERS FOR GAP CLOSURE...\n";
			
			my $qq = "$new_query_file\_$new_reference.fasta";
                        #print $path_toPass, "\t", $chk_uniq, "....\n"; exit;
			
                        if ($debug) {print "Passing inputs ...\n";}
                      
			primerP::pickPrimers($qq, $ref_loc, $flank, $path_toPass, $chk_uniq,$qualAvailable,@Quality_Array );	
		}
	}
}
################################################################################


