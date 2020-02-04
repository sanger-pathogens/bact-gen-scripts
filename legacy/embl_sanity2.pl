#!/usr/bin/perl

use lib "/software/pathogen/external/lib/perl/lib/site_perl/5.8.8";
#use lib "/Users/tc7/Software/perl/lib/perl5";
use Bio::SeqIO;
use Bio::Seq;


	$efile = $ARGV[0];
	$efile2=$ARGV[1];
	print $efile;	
	

my $embl = Bio::SeqIO->new(-file => $efile, -format => "EMBL");
my $seq_out = Bio::SeqIO->new('-file' => ">$efile2",
                                       '-format' => "EMBL");

         # write each entry in the input file to the output file
         while (my $inseq = $embl->next_seq) {
            $seq_out->write_seq($inseq);
         }


