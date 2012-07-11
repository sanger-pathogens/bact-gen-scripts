#!/usr/local/bin/perl -w

use strict;
use Bio::PSU::SeqFactory;

my $file = shift;
my $seqi = Bio::PSU::SeqFactory->make(-file   => $file,  
				      -format => 'embl');
my $obj=Bio::PSU::IO::FTHandler->new;

my(@gen1_gen2_ortho, @gen1_only_ortho, @gen2_only_ortho, @no_ortho, $feature );




# genome1 first comparator genome (gen1), second comparator genome (gen2) are those you are comparing to your genome $gen1 and $gen2 relate to the identifiers you have used to store the orthologue information. You must also fill in your generic gename pattern YE**** and the qualifier under which it is held eg locus_tag or systematic_id.

#first ortho qualifier
my $gen1= "E2348" . "_orthologue";
#second ortho qualifier
my $gen2= "K12" . "_orthologue";
#sytematic ID
my $genid= "ROD";
#my $genqualifier="locus_tag";
my $genqualifier="systematic_id";
  
                        while (my $seq = $seqi->next_seq)
               
                               {
				   if (my $fnum = $seq->has_features)
				   {
				       foreach $feature ($seq->features)
					   
				       { 
					   if ($feature->key() eq "CDS")
					   {
					       
					       
					       
					       
					       if ($feature->qexists("$gen1") and $feature->qexists("$gen2")) {
						   
						 #  push (@BPP_BB_ortho, $feature);
						    push (@gen1_gen2_ortho, $feature);
						   
						   
					      
					       
					       } elsif  ($feature->qexists("$gen1") xor  $feature->qexists("$gen2")) 
						   
					       {
						   
						   if ($feature->qexists("$gen1")) {
						       
						       
						       push (@gen1_only_ortho, $feature);
						       
						       
						   } else { 
						       
						       push (@gen2_only_ortho, $feature);
						       
						       
						   }
					       
					   
						   
						       
					       } else { 
							 
						       push (@no_ortho, $feature);
						   
						   }
					       
					   }
				       }
				   }
			       }
my $n=0;	
my $m=0;
my $l=0;
my $y=0;

print "*****************Your genome ( $genid ) unique genes ******************************:\n";

foreach my $feature1 (@no_ortho) {

    my @genelist = $feature1->qvalues("$genqualifier");	
    (my $generic_genename) = grep /$genid\d+/, @genelist;	
    print "$generic_genename\n";
    $n++;

}


print "****************genes shared with genome1 ( $gen1 ) and genome2 ( $gen2 ) *****************:\n";

			       
foreach my $feature2 (@gen1_gen2_ortho) {
    
    my @genelist = $feature2->qvalues("$genqualifier");
    my @hit1list = $feature2->qvalues("$gen1");
    my @hit2list = $feature2->qvalues("$gen2");
    (my $generic_genename) = grep /$genid\d+/, @genelist;
    (my $gen1_name) = $hit1list[0];
    (my $gen2_name) = $hit2list[0];	
    print "$generic_genename\t$gen1_name\t$gen2_name\n";
    $m++;
}


print "*****************genes shared with genome 1 ( $gen1 ) only******************************:\n";
       
foreach  my $feature3 (@gen1_only_ortho) {
    
    my @genelist = $feature3->qvalues("$genqualifier");
    my @hitlist = $feature3->qvalues("$gen1");	
    (my $generic_genename) = grep /$genid\d+/, @genelist;
    (my $gen1_name) = $hitlist[0];
    print "$generic_genename\t$gen1_name\n";
    $l++;

}



print "****************genes shared with genome 2 ( $gen2 )only *******************************:\n";

			       
foreach my $feature4 (@gen2_only_ortho) {

    my @genelist = $feature4->qvalues("$genqualifier");
    my @hitlist = $feature4->qvalues("$gen2");
    (my $generic_genename) = grep /$genid\d+/, @genelist;
    (my $gen2_name) = $hitlist[0];
    print "$generic_genename\t$gen2_name\n";
    $y++;

}

print "********************************************************************************\n";

print "there are $n unique ( $genid ) genes \n";
print "there are $m genes shared with genome1( $gen1 ) and genome2 ( $gen2 )\n";
print "there are $l genes shared with genome1 ( $gen1 ) only \n";
print "there are $y genes shared with genome2 ( $gen2 ) only \n"
