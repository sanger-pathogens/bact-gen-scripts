#!/usr/bin/env	perl

use lib "/software/pathogen/external/lib/perl/lib/site_perl/5.8.8";
use Bio::SeqIO;
use Getopt::Long;

my $nmode = 0;
my $ev = 0.1;
my $pid = 50;
my $graph = 0;
my $bm = 0;
my $protein = 0;

Getopt::Long::Configure ('bundling');
GetOptions ('i=s' => \$infile,
            'f=s' => \$format,
            'b=s' => \$blastdir,
            'o=s' => \$outputfile1,
            'm' => \$mode,
            'n' => \$nmode,
            'g' => \$graph,
            'P' => \$protein,
            'e=s' => \$ev,
            'p=s' => \$pid,
            'B=s' => \$bm,
            't' => \$test,
            'h' => \$help) or help_msg(); 

if($help == 1)
{
	help_msg();
	
	exit 0;
}

# script to take two files, create a blast database and then blast the embl against the fasta file
my $out;
if($format eq "embl")
{
 cds_tofasta($infile);
 $out = $infile.".concat.fas";	
}
else
{
 $out = $infile; 	
}



$seqio_obj = Bio::SeqIO->new(-file => $out, -format => $format );
$out2 = $out."translated.fas";
$seqio_obj2 = Bio::SeqIO->new(-file => ">$out2", -format => $format ); 
				if($nmode == 1)
				{
					$out2 = $out;
				}
				else
				{
				while ( my $seq = $seqio_obj->next_seq() ) {
			if ( $protein == 0){
    			$seq = $seq->translate;
			}
    			$seqio_obj2->write_seq($seq);
    			}  	
				}

$out = $out2;


# so we now have a reference file, is time to do the blast searches

# create the blast databases
# and then run the blast searches

$jobno = int(rand(64000));
$scriptarray1 = int(rand(64000)).".sh";
$scriptarray2 = int(rand(64000)).".sh";
$scriptarray1v = int(rand(64000))."sh";
$scriptarray2v = int(rand(64000))."sh";


 			open(FILE, ">>$scriptarray2");
     		print FILE "\#! /bin/bash\n\#\$LSB_JOBINDEX\n\# create an array which has the commands to run\n declare -a runme\n";
     		close (FILE);     		
     		
$startbit = "bsub -J ".$jobno;

 @files = <$blastdir/*.nsq>;
 $filecount = 1;
my $keyfile = $outputfile1."-taxakey.csv";
my $taxanm = "";
 foreach $file (@files) {
  	$file =~ s/.nsq//;
  	print $file;
	$taxanm = $file;
	$taxanm =~ s/$blastdir//g;
	$taxanm =~ s/\///g;     
  	#formatdb
  	system("ln -s ".$file.".nsq tempblast".$taxanm.".nsq"); 
  	system("ln -s ".$file.".nin tempblast".$taxanm.".nin");
  	system("ln -s ".$file.".nhr tempblast".$taxanm.".nhr");
  	system("ln -s ".$file.".dna tempblast".$taxanm.".dna");

     open(FILE, ">>$scriptarray2");
     if($nmode == 1)
   	{
   		print FILE "runme[$filecount]=\"-p blastn -d tempblast".$taxanm." -i ".$out." -m0 -e ".$ev." -F F -o tempblast".$taxanm.$jobno.".out\"\n";

   	}
   	elsif($graph == 1)
   	{
   		print FILE "runme[$filecount]=\"-p tblastn -d tempblast".$taxanm." -i ".$out." -m0 -e ".$ev." -F F -o tempblast".$taxanm.$jobno.".out\"\n";
   	}
   	else
   	{
   		print FILE "runme[$filecount]=\"-p tblastn -d tempblast".$taxanm." -i ".$out." -m0 -e ".$ev." -F F -o tempblast".$taxanm.$jobno.".out\"\n";

   	}
   	
          close (FILE);	
     $filecount++;
  } 
     open(FILE, ">>$scriptarray2");
     print FILE 'echo blastall ${runme[$LSB_JOBINDEX]}'."\n".'blastall ${runme[$LSB_JOBINDEX]}';
     close (FILE);
     
  system("chmod 755 ".$scriptarray2);
  
  if($bm != 0)
  {
  system("bsub -R\"select[mem>".$bm."000] rusage[mem=".$bm."000]\" -M ".$bm."000 -J \"".$scriptarray2v."[1-$filecount]\" -o blasting.out%J-%I ./".$scriptarray2);
  }
  else
  {
  system("bsub -R\"select[mem>2000] rusage[mem=2000]\" -M 2000 -J \"".$scriptarray2v."[1-$filecount]\" -o blasting.out%J-%I ./".$scriptarray2);  	
  }
  $scripttidy = "scripttidy.sh";
  unlink($scripttidy);
  print "bsub -R\"select[mem>2000] rusage[mem=2000]\" -M 2000 -J \"".$scriptarray2v."[1-$filecount]\" -o blasting.out%J-%I ./".$scriptarray2."\n";
  print "bsub -R\"select[mem>".$bm."000] rusage[mem=".$bm."000]\" -M ".$bm."000 -J \"".$scriptarray2v."[1-$filecount]\" -o blasting.out%J-%I ./".$scriptarray2."\n";
   open(FILE, ">>$scripttidy");
   print FILE "\#!/bin/bash\nrm -rf ".$scriptarray1."\nrm -rf ".$scriptarray2."\nrm -f *.nsq\nrm -f *.nin\nrm -f *.nhr\nrm -f *.dna\n exit \$?";  	
   close (FILE);
  # system("chmod 775 ".$scripttidy);
  print "bsub -w \"ended(".$scriptarray2v.")\" ./".$scripttidy."\n";
   system("bsub -w \"ended(".$scriptarray2v.")\" ./".$scripttidy);
   if($mode == 1)
   {
   		system("bsub -q long -w \"ended(".$scriptarray2v.")\" -R\"select[mem>10000] rusage[mem=10000]\" -M 10000 perl /nfs/users/nfs_t/tc7/scripts/pub/tblastparse3_stats.pl ".$outputfile1." ".$blastdir);
   		print "bsub -q long -w \"ended(".$scriptarray2v.")\" -R\"select[mem>10000] rusage[mem=10000]\" -M 10000 perl /nfs/users/nfs_t/tc7/scripts/pub/tblastparse3_stats.pl ".$outputfile1." ".$blastdir."\n";
   }
   else
   {
   	if($nmode == 1)
   	{
   		system("bsub -q long -w \"ended(".$scriptarray2v.")\" -R\"select[mem>10000] rusage[mem=10000]\" -M 10000 -o parseout /nfs/users/nfs_t/tc7/scripts/pub/nblastparse3.pl ".$outputfile1." ".$blastdir." ".$infile." ".$pid);
   	}
   	elsif($graph == 1)
   	{
   		system("bsub -q long -w \"ended(".$scriptarray2v.")\" -R \"select[mem>10000] rusage[mem=10000]\" -M 10000 -q long -o parseout /nfs/users/nfs_t/tc7/scripts/pub/gblastparse3.pl ".$outputfile1." ".$blastdir." ".$infile." ".$pid);
   	}
   	else
   	{
   		system("bsub -q long -w \"ended(".$scriptarray2v.")\" -R \"select[mem>10000] rusage[mem=10000]\" -M 10000 -o parseout /nfs/users/nfs_t/tc7/scripts/pub/tblastparse3.pl ".$outputfile1." ".$blastdir." ".$infile." ".$pid);
   	}
   }
   sub cds_tofasta {
   	my $efile=$_[0];
my $embl = Bio::SeqIO->new(-file => $efile, -format => "EMBL");
my $annotation = $embl->next_seq;
# system("formatdb -i ".$file." -p F"); 
# create fasta query files
my $maxsqes = 0;
my $out = $efile.".concat.fas";
foreach $gene ($annotation->get_SeqFeatures) 
{
	if ($gene->primary_tag eq "CDS") 
	{
		$seq = $gene->spliced_seq->seq;
		if ($gene->has_tag('locus_tag')) 
		{
         for my $val ($gene->get_tag_values('locus_tag'))
         {
            $value = $val;
            # e.g. 'NDP', from a line like '/gene="NDP"'
         }
		}
		elsif ($gene->has_tag('systematic_id')) 
		{
         	for my $val ($gene->get_tag_values('systematic_id'))
         	{
            	$value = $val;
            	# e.g. 'NDP', from a line like '/gene="NDP"'
         		}
			}
			else
			{
				$value = $maxsqes;
			}
				$maxsqes++;			
 				open(FILE, ">>$out");
     			print FILE ">".$value."\n".$seq."\n";
     			close (FILE);
 				
		}
	}
   	
   }

sub help_msg
{

	print "mlsttblastn.pl Options;
			-i input file with genes/features to search,
            -f format of input file (fasta/embl),
            -b directory of all your blast databases,
            -o output filename,
            -m add this to just get the identity info from the samples
            -n run with blastn
            -g Do a comparison to produce a graph of hits from a big database,
	    -P input files are proteins
            -e minimum e value required,
            -p minimum identity % length of the sequence required to hit to be included,            
            -h show this message\n";
           
} 
