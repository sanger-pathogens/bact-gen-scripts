#!/usr/bin/perl

$filenameA = $ARGV[0];
$filenameOut = $ARGV[1];

open $FILEA, "< $filenameA";

open $OUTFILE, "> $filenameOut";

while(<$FILEA>) {
	$_ =~ s/\@/\>/g;
	print $OUTFILE $_;
	$_ = <$FILEA>;
	print $OUTFILE $_; 
	$_ = <$FILEA>;
	$_ = <$FILEA>;

}
