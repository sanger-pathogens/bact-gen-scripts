#!/usr/bin/perl -w
use strict;


if (@ARGV != 1) {
    print "$0 fastq.gz|fastq\n" ; 
	exit ;
}



if ( $ARGV[0] =~ /.gz/ ) {
    open (IN, "zcat $ARGV[0] | ") or die "can not open $ARGV[0]\n" ;
}
else {
    open (IN , "$ARGV[0]") or die "can not open $ARGV[0]\n" ;
}

open my $out, ">" , "$ARGV[0].fasta" or die "oops\n" ;
open my $out2, ">" , "$ARGV[0].fasta.qual" or die "oops\n" ;

while (<IN>)
{
  	if (/^@(\S+)/) {
		print $out ">$1\n" ;
		print $out2 ">$1\n" ;

		my $read = <IN> ;
		print $out "$read" ;
		$read = <IN> ;
		#print "+\n" ;
		$read = <IN> ;
		#print "$read" ;

		chomp($read) ;
		my @read_quals = split "" , $read ;
		my @read_quals_converted = map( ord($_) - 33 , @read_quals) ;
		print $out2 "@read_quals_converted\n" ;

	}

}


