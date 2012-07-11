#!/usr/bin/perl -w


$filenameA = $ARGV[0];
$filenameB = $ARGV[1];
$filenameOut = $ARGV[2];

if ($filenameA =~ /\.gz$/){
	open $FILEA, "gunzip -c $filenameA |" or die $!;
}
else {
	open $FILEA, "< $filenameA";
}
if ($filenameA =~ /\.gz$/){
	open $FILEB, "gunzip -c $filenameB|" or die $!;
}
else {
	open $FILEB, "< $filenameB";
}

open $OUTFILE, "> $filenameOut";

while(<$FILEA>) {
	print $OUTFILE $_;
	$_ = <$FILEA>;
	print $OUTFILE $_; 
	$_ = <$FILEA>;
        print $OUTFILE $_;
	$_ = <$FILEA>;
        print $OUTFILE $_;	

	$_ = <$FILEB>;
	print $OUTFILE $_; 
	$_ = <$FILEB>;
	print $OUTFILE $_;
	$_ = <$FILEB>;
        print $OUTFILE $_;
        $_ = <$FILEB>;
        print $OUTFILE $_;
}
