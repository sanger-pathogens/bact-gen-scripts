#!/usr/local/bin/perl -w
# Copyright (C) 2009 Genome Research Limited. All Rights Reserved.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#ABACAS.1.1.1
#--------------------------------
# UPDATES: This version inserts 100 N's between overlapping contigs and creats a feature file for gaps

#Please report bugs to:
#sa4@sanger.ac.uk & tdo@sanger.ac.uk
 

use strict;
use DaemonsLog::Log;
use warnings;
use POSIX qw(ceil floor);
use Getopt::Std;

#-------------------------------------------------------------------------------

if (@ARGV < 1) { usage();}

my ($help, $reference, $query_file, $choice, $sen, $seed, $mlfas, $fasta_bin, $avoid_Ns,
   $tbx, $min_id, $min_cov, $diff_cov, $min_len,$add_bin_2ps, $pick_primer,
   $flank, $chk_uniq,$redoMummer, $is_circular,$escapeToPrimers, $debug, $gaps_2file, $prefix)
   =checkUserInput( @ARGV );

my $ref_inline; 
pickPrimers ($reference, $query_file, $flank, $chk_uniq);
exit;

#----------------------------------------------- PRIMER DESIGN ---------------------------------------------------
sub pickPrimers
{
          #$ps = pseudo molecule,$rf = reference, $flan = flanking region size
          my ($ps,$rf, $flan, $passed_path, $chk_uniq,$qualAvailable, @Quality_Array);
          if (@_==4){
                  ($rf,$ps, $flan, $chk_uniq) = @_;
                   print "Primers without ordering..\n";
                   print $rf;
                   $passed_path = "nucmer";
                   $qualAvailable =0;
                   @Quality_Array = [];
          }
          else #(@_ == 7)
          {
                ($ps,$rf, $flan, $passed_path, $chk_uniq, $qualAvailable, @Quality_Array) = @_;
          }
          
          my $dna='';
          my @gappedSeq;
          my $records='';
          my @sequence;
          my $input='';
	  #tdo
	  my @gappedQual;
	  #my $quality='';
          my $path_toPass = $passed_path;
          my @fasta;
          my $query='';
          my @exc_regions;
          my $ch_name;
          #my $flank = $flan; 
          open (FH, $rf) or die "Could not open reference file\n";
          open (FH2, $ps) or die "Could not open query/pseudomolecule file\n";
          my $ref; #print ".... ", $rf; exit; 
          my @r = <FH>;
          my @qry = <FH2>;
          my $dn = join ("", @qry);
          $ref = join ("", @r);
          $dna = Fasta2ordered ($dn);
          #check if primer3 is installed
          my $pr3 = "primer3_core";
          my ($pr3_path, $pr3_Path, $pr3_path_dir) = checkProg ($pr3);
          #my @check_prog = `which primer3_core`;     
          open (PRI, '>primer3.summary.out') or die "Could not open file for write\n";
          
         # 
          
          #print $ref; exit;
          #PARSING FOR PRIMER3 INPUT
          my ($opt,$min,$max,$optTemp,$minTemp,$maxTemp,$flank,$lowRange,$maxRange,$gcMin,$gcOpt,$gcMax,$gclamp,$exclude,$quality) = getPoptions($qualAvailable); 
          my ($gap_position,@positions, %seq_hash);
          
          my $exc1 = $flank -$exclude;  #start of left exclude
          print "Please wait... extracting target regions ...\n";
          #regular expression extracts dna sequence before and after gaps in sequence (defined by N)
          while($dna=~ /([atgc]{$flank,$flank}N+[atgc]{$flank,$flank})/gi)
          {               
	        $records= $1;
			push (@gappedSeq, $records);
			$gap_position = index($dna, $records);
			push @positions, $gap_position;
			$seq_hash{$gap_position}=$records;
			#dna
			if ($qualAvailable) {
			  my $res;
			  for (my $nn=($gap_position-1); $nn <= ($gap_position-1+length($records)-1); $nn++) {
				$res.="$Quality_Array[$nn] ";
			  }
			  push @gappedQual, $res;			  
			}
          }
        #loop prints out primer targets into a file format accepted by primer3
        my $count=1;
        my $identify='';
        my $seq_num = scalar @gappedSeq;
        my $name= " ";
        
        my ($totalp, @left_name, @right_name, @left_seq, @right_seq);

        my ($leftP_names, $rightP_names, $leftP_seqs, $rightP_seqs, $leftP_start, $leftP_lens, $rightP_ends, $rightP_lens,$left_Exclude,$right_Exclude, $primers_toExc, $prod_size)=
        ("","","","","","","","","","", "", "");    
        
        print $seq_num, " gaps found in target sequence\n"; 
        print "Please wait...\nLooking for primers...\n";
        print "Running Primer3 and checking uniquness of primers...\nCounting left and right primer hits from a nucmer mapping (-c 15 -l 15)\n";
        
        for (my $i=0; $i<$seq_num; $i+=1)
        {            
                    $identify = $count++;
                    if (defined $ch_name)
                    {
                            $name = $ch_name;
                    }
                    my $len = length($gappedSeq[$i]);
                    my $exc2 = $len - $flank;
                    open(FILE, '>data') or die "Could not open file\n";
                    #tdo
                    my $qual='';
                    if ($qualAvailable) {
                      $qual="PRIMER_SEQUENCE_QUALITY=$gappedQual[$i]\nPRIMER_MIN_QUALITY=$quality\n";
                    }
				
#WARNING: indenting the following lines may cause problems in primer3 
print FILE "PRIMER_SEQUENCE_ID=Starting_Pos $positions[$i] 
SEQUENCE=$gappedSeq[$i]
PRIMER_OPT_SIZE=$opt
PRIMER_MIN_SIZE=$min
PRIMER_MAX_SIZE=$max
PRIMER_OPT_TM=$optTemp
PRIMER_MIN_TM=$minTemp
PRIMER_MAX_TM=$maxTemp
PRIMER_NUM_NS_ACCEPTED=1
PRIMER_PRODUCT_SIZE_RANGE=$lowRange-$maxRange
PRIMER_MIN_GC=$gcMin
PRIMER_GC_CLAMP =$gclamp
PRIMER_OPT_GC_PERCENT=$gcOpt
PRIMER_MAX_GC=$gcMax
PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION=$exc1,$exclude $exc2,$exclude
".$qual."Number To Return=1
=\n";
close FILE;

        #runs primer3 from commandline 
        
        ################# NOTE: PRIMER3 SHOULD BE IN YOUR WORKING PATH #########
      
                    my  @Pr3_output = `$pr3_path -format_output <data`;
                    #print $positions[$i], "\t", $i, " ", $path_toPass, "  ", $rf, $exc1, " ",$exc2, "\n";
                    my $fil = join (":%:", @Pr3_output);
                   my ($uniq_primer, $string,$left_nm,$right_nm,$left_sq, $right_sq,$left_strt,$left_ln, $right_End,$right_ln,$primers_toExclude, $product_size)
                    = check_Primers ($fil, $positions[$i], $i,$path_toPass, $rf, $exc1, $exc2);
                    
                    print PRI $string;
                    if ($uniq_primer ==1)
                    {
                              $leftP_names.=$left_nm."\n";
                              $rightP_names.=$right_nm."\n";
                              $leftP_seqs.=$left_sq."\n";
                              $rightP_seqs.=$right_sq."\n";
                              $leftP_start.=$left_strt."\n";
                              $leftP_lens.=$left_ln."\n";
                              $rightP_ends.=$right_End."\n";
                              $rightP_lens.=$right_ln."\n";
                              $left_Exclude.=$exc1."\n";
                              $right_Exclude.=$exc2."\n";
                              $prod_size.=$product_size."\n";
                    }
                    if ($primers_toExclude ne "")
                    {
                              $primers_toExc.= $primers_toExclude; #."\n";
                    }
          
          }
        write_Primers ($leftP_names, $rightP_names, $leftP_seqs, $rightP_seqs, $leftP_start, $leftP_lens, $rightP_ends, $rightP_lens,$primers_toExc,$left_Exclude,$right_Exclude, $prod_size);     
        #write_Primers (@left_name, @right_name, @left_seq, @right_seq,@left_start, @left_len, @right_end, @right_len, @left_exclude, @right_exclude, $primers_toExclude);
        
}

#checks the uniqueness of primers
#input an array with promer3 output for each gap
sub check_Primers
{
          
          my ($fil, $position, $index,$path_toPass, $rf, $exc1, $exc2) = @_;
          my  @Pr3_output = split /:%:/, $fil;           
          my ($left_name, $right_name, $left_seq, $right_seq, $left_start,$left_len,$right_end,$right_len,$left_exclude,$right_exclude) = ("", "", "", "", "", "", "", "", "", "");
          my $primers_toExclude ="";
          my $product_size ="";
          my $string ="";
          my $uniq_primer = 0;
          $string.="=========================================================================================\n";
          $string.="Primer set for region starting at ".$position."\n";
          
          if (defined $Pr3_output[5] && defined $Pr3_output[6])
          {
                    if ($Pr3_output[5]=~ /LEFT PRIMER/)
                    {
                           # print $Pr3_output[5];
                            #check uniquness of primer against the genome
                            my @splits_1 = split (/\s+/, $Pr3_output[5]);
                            my $left_primer = $splits_1[8];
                            my $left_st = $splits_1[2];
                            my $left_length = $splits_1[3];
                            
                            my @splits_2 = split (/\s+/, $Pr3_output[6]);
                            my $right_primer = $splits_2[8];
                            my $right_st = $splits_2[2];
                            my $right_length = $splits_2[3];
                            
                            open (QRY_1, '>./left_query'); # open a file for left primers
                            print QRY_1 ">left_01\n";    #
                            print QRY_1 $left_primer,"\n";
                            open (QRY_2, '>./right_query');
                            print QRY_2 ">right_01\n";
                            print QRY_2 $right_primer,"\n";
                            
                            my ($left_count, $right_count);
                            #if ($chk_uniq eq "nucmer")
                            #{
                                my $options = "-c 15 --coords  -l 15 ";
                                my $rq = "right_query";
                                my $lq = "left_query";
                                my (@right_ps, @left_ps);
                               # print $path_toPass, "\t", $options, "\n";
                                
                        
                                my @Rrun = `$path_toPass $options -p R $rf  $rq &> /dev/null`;
                                print ".";
                                my $f1 = "R.coords";
                                open (RP, $f1) or die "Could not open file $f1 while checking uniqueness of right primer\n";                           
                                  while (<RP>)
                                {
                                    chomp;
                                    if ($_ =~ /right_01$/)
                                    {
                                        push @right_ps, $_;
                                    }
                                }
                                close (RP);
                                my @Lrun = `$path_toPass $options -p L $rf  $lq &> /dev/null`;                               
                                print ".";
                                my $f2 = "L.coords";
                                open (LQ, $f2) or die "Could not open file $f2\n";
                                while (<LQ>)
                                {
                                    chomp;
                                    if ($_ =~ /left_01$/)
                                    {
                                        push @left_ps, $_;
                                    }
                                }
                                close (LQ);
                                $right_count = scalar (@right_ps);  
                                $left_count = scalar(@left_ps); 
                              #check if a primer is not in the excluded region::
                              my $primer_NearEnd =0;
                              if ($left_st > $exc1 || $right_st < $exc2)
                              {
                                        $primer_NearEnd = 1;
                              }
                              
                            if ($left_count < 2 && $right_count<2 && $primer_NearEnd ==0)
                            {
                                    $string.=$left_count."\t".$Pr3_output[5]."\n";
                                    $string.=$right_count."\t".$Pr3_output[6]."\n";
                                    $string.="***************************** PRIMER3 OUTPUT **************************\n";
                                    foreach (@Pr3_output) {$string.=$_;}
                                    
                                    my @prod_size_split = split /\s+/, $Pr3_output[10];
                                    
                                    $product_size = substr($prod_size_split[2], 0, -1);
                                    $left_name = $position;
                                    $right_name = $position;
                                    my $lp_uc = uc ($left_primer);
                                    my $rp_uc = uc($right_primer);
                                        #print $left_count, "..", $right_count, "\t";
                                    $left_seq = $lp_uc;
                                    $right_seq= $rp_uc;
                                    
                                    $left_start= $left_st;
                                    $left_len = $left_length;
                                    
                                    $right_end = $right_st;
                                    $right_len = $right_length;
                                    
                                    $left_exclude = $exc1;
                                    $right_exclude =$exc2;
                                    $uniq_primer =1;
                            }
                            else
                            {
                                        if ($primer_NearEnd ==1)
                                        {
                                             $string.="One of the oligos is near the end of a contig\n";     
                                        }
                                        else
                                        {
                                                  $string.="Primer set not unique\n";
                                        }
                                    $primers_toExclude.=">L.".$position."\n".$left_primer."\n";
                                    $primers_toExclude.=">R.".$position."\n".$right_primer."\n";
                            }
                            
                    }
                    else
                    {
                            $string.="No Primers found\n";
                    }
          }
                    
          return ($uniq_primer, $string,$left_name,$right_name,$left_seq, $right_seq,$left_start,$left_len, $right_end,$right_len,$primers_toExclude, $product_size);
          
          
}

###------------------------------------
# Writes primers and their regions to file
sub  write_Primers {
          my ($leftP_names, $rightP_names, $leftP_seqs, $rightP_seqs, $leftP_start, $leftP_lens, $rightP_ends, $rightP_lens,$primers_toExclude,$left_Exclude,$right_Exclude, $product_sizes) = @_;
          my (@left_name, @right_name, @left_seq, @right_seq, @left_start, @left_len, @right_end, @right_len, @left_exclude, @right_exclude, @product_size);
         
          #open files to read
          @left_name = split /\n/, $leftP_names;
          @right_name= split /\n/, $rightP_names;
          @left_seq = split /\n/, $leftP_seqs;
          @right_seq = split /\n/, $rightP_seqs;
          
          @left_start = split /\n/, $leftP_start;
          @left_len = split /\n/, $leftP_lens;
          @right_end = split/\n/, $rightP_ends;
          @right_len = split /\n/, $rightP_lens;
          @left_exclude = split /\n/, $left_Exclude;
          @right_exclude = split /\n/,$right_Exclude;
          @product_size = split /\n/, $product_sizes;
            
          my $primers_withSize ="";       
          open (SEN, '>sense_primers.out') or die "Could not open file for write\n";
          open (ASEN, '>antiSense_primers.out') or die "Could not open file for write\n";
          open (REG_1, '>sense_regions.out') or die "Could not open file for write\n";
          open (REG_2, '>antiSense_regions.out') or die "Could not open file for write\n";
          
          if ($primers_toExclude ne "")
          {
                    open (PEX, '>primers_toExclude.out') or die "Could not open file for write\n";
                    print PEX $primers_toExclude;
          }
          
          
          my $totalp = scalar (@left_name);
          
          #print $totalp, "\n"; exit;
          
          my $well_pos;
          my $max_plates = ceil($totalp/96);
          #print "MAX Ps ", $max_plates, "\n";
          my $plate=1;
          my $sen ="";
          my $asen ="";
          my $plate_counter =0;
          my $wells = 96;
          for (my $index =0; $index < $totalp; $index += $wells)
          {
                   my $do = $index;
                   my $upper_bound= $index + $wells;
                   if ($upper_bound > $totalp)
                   {
                              $upper_bound = $totalp;
                   }
                   
                   for (my $j=$index; $j <= ($upper_bound-1); $j+=1)
                    {
                       my $i = $j;
                       if ($j < 96)
                       {
                              $well_pos = get_WellPosition ($j);
                       }
                       else
                       {
                            $well_pos = get_WellPosition ($j - $wells)  
                       }
                          
                           #$primers_withSize.=$product_size[$i]."\t"."Plate_".$plate. "\t\tS.".$i."\tS.".$left_name[$i]."\t".
                           print SEN "Plate_".$plate, "\t\t","S.", $i, "\tS.", $left_name[$i], "\t", $left_seq[$i], "\t\t+", "\t", $well_pos, "\n"; 
                           print ASEN "Plate_".$plate, "\t\t","AS.", $i, "\tAS.", $right_name[$i], "\t", $right_seq[$i], "\t\t-","\t", $well_pos,"\n";
                           print REG_1 "Plate_".$plate, "\t\t","S.", $i, "\t", $left_start[$i], "\t", $left_len[$i], "\n";
                           print REG_2 "Plate_".$plate, "\t\t","AS.", $i, "\t", $right_end[$i], "\t",$right_len[$i], "\n";
                       
                    }
                    $plate +=1;
          }
             
        #delete tmp. files
        #my $rm = "rm -f";
        system ("rm -f data left_query right_query R.delta R.cluster R.coords L.delta L.cluster L.coords");
        print "\nPRIMER DESIGN DONE\n\n";
	# end of primer design program
}#//
#####
# returns a well position for oligos 
sub get_WellPosition{
          
          my $j = shift;
          my $well_pos;
          if ($j < 12)
          {
                    $well_pos = "a".($j+1);
          }
          elsif ($j>11 && $j<24) {
                    $well_pos = "b". (($j+1) -12);
          }
          elsif ($j>23 && $j<36) {
                    $well_pos = "c". (($j+1) -24);
          }
          elsif ($j>35 && $j<48) {
                    $well_pos = "d". (($j+1) - 36);
          }
          elsif($j>47 && $j<60) {
                    $well_pos = "e". (($j+1) -48);
          }
          elsif ($j>59 && $j<72)
          {
                    $well_pos = "f". (($j+1) - 60);
          }
          elsif ($j>71 && $j< 84)
          {
                    $well_pos = "g". (($j+1) - 72);
          }
          elsif ($j>83 && $j<96)
          {
                    $well_pos = "h". (($j+1) - 84);
          }
          return $well_pos;
}


####################################################################
#get options for primer design
#----------------------
sub getPoptions{
          
          my $qualAvailable = shift;
          #### USER INPUTS ##########
        #ask for optimum primer size
        print "\nEnter Optimum Primer size (default 20 bases):";
        my $opt=<STDIN>;
        chomp $opt;
        if($opt eq '')
        {
                $opt = 20;
        }
        #ask for minimum primer size
        print "\nEnter Minimum Primer size (default 18 bases):";
        my $min=<STDIN>;
        chomp $min;
        if($min eq '')
        {
                $min = 18;
        }
        #ask for maximum primer size
        print "\nEnter Maximum Primer size (default 27 bases):";
        my $max= <STDIN>;
        chomp $max;
        if($max eq '')
        {
                $max= 27;
        }
        #ask for optimum primer temperature
        print "\nEnter Optimum melting temperature (Celcius) for a primer oligo (default 60.0C):";
        my $optTemp=<STDIN>;
        chomp $optTemp;
        if($optTemp eq '')
        {
                $optTemp = 60.0;
        }
        #ask for minimum primer temperature
        print "\nEnter Minimum melting temperature (Celcius) for a primer oligo (default 57.0C):";
        my $minTemp=<STDIN>;
        chomp $minTemp;
        if($minTemp eq '')
        {
                $minTemp = 57.0;
        }
        #ask for maximum primer temperature
        print "\nEnter Maximum melting temperature (Celcius) for a primer oligo (default 63.0C):";
        my $maxTemp=<STDIN>;
        chomp $maxTemp;
        if($maxTemp eq '')
        {
                $maxTemp = 63.0;
        }
        print "\nEnter flanking region size (default 1000 bases): ";
        my $flank=<STDIN>;
        chomp $flank;
        if ($flank eq '')
        {
          $flank = 1000;
        }
        #ask for primer product range
        print "\nEnter minimum product size produced by primers (default =flanking size):";
        my $lowRange=<STDIN>;
        chomp $lowRange;
        if($lowRange eq '')
        {
                $lowRange = $flank;
        }
        print "\nEnter maxmimum product size produced by primers (default 7000):";
        my $maxRange=<STDIN>;
        chomp $maxRange;
        if($maxRange eq '')
        {
                $maxRange = 7000;
        }
        #ask for minimum GC content in primers
        print "\nEnter minimum GC content in primers (default 20%):";
        my $gcMin=<STDIN>;
        chomp $gcMin;
        if($gcMin eq '')
        {
                $gcMin = 20.0;
        }
        #ask for optimum GC content in primers
        print "\nEnter optimum GC content in primers (default 50%):";
        my $gcOpt=<STDIN>;
        chomp $gcOpt;
        if($gcOpt eq '')
        {
                $gcOpt = 50.0;
        }
        #ask for maximum GC content in primers
        print "\nEnter maximum GC content in primers (default 80%):";
        my $gcMax=<STDIN>;
        chomp $gcMax;
        if($gcMax eq '')
        {
                $gcMax = 80.0;
        }
        print "\nEnter GC clamp  (default 1):";
        my $gclamp=<STDIN>;
        chomp $gclamp;
        if($gclamp eq '')
        {
                $gclamp = 1;
        }
	print "\nEnter size of region to exclude at the end of contigs (default 100 bases):";
	my $exclude=<STDIN>;
	chomp $exclude;
	if ($exclude eq '')
	{
		$exclude = 100;
	}


	  #tdo
          my $quality='';
	  if ($qualAvailable)
          {
		
                    print "\nEnter minimum quality for primer pick (default 40):";
                    $quality=<STDIN>;
                    chomp $quality;
                    if($quality eq '')
		  {
			$quality = 40;
		  }
	  }

          
return ($opt,$min,$max,$optTemp,$minTemp,$maxTemp,$flank,$lowRange,$maxRange,$gcMin,$gcOpt,$gcMax,$gclamp,$exclude, $quality);
          
}
###############
#-----------------------------------------------------END of PRIMER DESIGN ----------------------------------------------------------------
#-----------------------------------------------------END OF ABACAS -----------------------------------------------------------------------


