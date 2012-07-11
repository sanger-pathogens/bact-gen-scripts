#! /usr/bin/perl -w
#
# File: cigar00_2plot.pl
# Time-stamp: <29-Apr-2009 16:20:41 tdo>
# $Id: $
#
# Copyright (C) 2009 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto
#
# Description:
#


if (scalar(@ARGV)<3) {
  die "usage: <cigar file> <resultname> <readlength> optonal:<bases to cut on both ends>";
  
}
my $cigar=shift;
my $resultName=shift;
my $readLength =shift;


my $BORDER=shift;
if (!defined($BORDER)) {
  $BORDER=0;
}

my %PLOT;

if( $cigar =~ /\.gz$/ )
  {
	open(F, "gunzip -c $cigar |" ) or die "Cant open cigar file $cigar: $!\n";
  }
else {
  open (F,$cigar) or die "cigar not in\n";
}


#@F=<F>;

#foreach (@F) {
while (<F>) {
  my @ar=split(/\s+/);
  
  if (defined($ar[5])
	  # the read must map at least with readlength - 5 bps
	  &&	  ( (abs($ar[3]-$ar[2])+1)>= ($readLength - 5))
	 ){
	foreach (($ar[6]+$BORDER)..($ar[7]-$BORDER)) {
	  $PLOT{$ar[5]}[($_-1)]++;
	}
  }
}

my $res;
foreach my $chr (keys %PLOT) {
  $res='';
  
  open (F,"> $resultName") or die "please give resultnam $resultName \n";
#  print "working on $chr\n";
  
  foreach my $i (@{ $PLOT{$chr} }) {
	if (defined($i)) {
	  $res.="$i\n";
	}
	else {
	  $res.="0\n";
	}
  }
  print F $res;
  close(F);
  
}



