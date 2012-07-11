#!/usr/bin/perl

use lib "/software/pathogen/external/lib/perl/lib/site_perl/5.8.8";
#use lib "/Users/tc7/Software/perl/lib/perl5";
use Bio::SeqIO;
use Getopt::Long;

# script for parsing a nexus file, tree file and producing an xml file and starting a BEAST analysis

# read in the files

# get the options

# produce the xml file


#read in the csv with the dates and write these out into a taxa list

#beast output filename
my $out;
#datefile
my $datefile;
#sequence file
my $sequencefile;
#sequence format
my $format;
#how frequently stuff is logged
my $logtime; 
#the filename for the outputs
my $filename; 
#how long to run for
my $chainlength; 
my $startingtree = "";
my $relaxed;
my $help;
my $patterns;
my $overwrite;
print "run_BEAST.pl.  Script to produce BEAST XML files, written by T R Connor, WTSI, 2011\n\n";
Getopt::Long::Configure ('bundling');
GetOptions ('b=s' => \$out,
            'd=s' => \$datefile,
            'i=s' => \$sequencefile,
            'f=s' => \$format,
            'l=i' => \$logtime,
            'o=s' => \$filename,
            'c=i' => \$chainlength,
            't=s' => \$startingtree,
            'p=s' => \$patterns,
            's' => \$relaxed,
            'r' => \$overwrite,
            'h' => \$help) or die print "run_BEAST.pl Options;\n -b beastfile\n -d datefile\n -s sequencefile\n -f format of sequence file\n -l log interval\n -o beast output filename\n -c chain length\n -t starting tree file\n -r (if you want to run with a relaxed clock)\n";

if($help == 1)
{
	help_msg();
	
	exit 0;
}

# check arguments
if(length($datefile) == 0)
{
print "Date File Missing; this argument is required in order to run beast\n";
help_msg();
exit 0;
}
else
{
	unless (-e $datefile) {
 	print $datefile." Does Not Exist!\n";
 	exit 1;
 	} 
}

if(length($sequencefile) == 0)
{
print "No Sequence File Specified; this argument is required in order to run beast\n";
help_msg();
exit 0;
}
else
{
	unless (-e $sequencefile) {
 	print $sequencefile." Does Not Exist!\n";
 	exit 1;
 	} 
}

if(length($filename) == 0)
{
print "No prefix for the BEAST logs/trees specified; this argument is required in order to run beast\n";
help_msg();
exit 0;
}



if(length($out) == 0)
{
print "No BEAST xml file specified; this argument is required in order to run beast\n";
help_msg();
exit 0;
}

if($chainlength == 0)
{
print "No chain length specified, using a chain length of 100,000,000 states\n";
$chainlength = 100000000;
}
if($logtime == 0)
{
print "No logtime specified, using log to file/screen once every 1000 states\n";
$logtime = 1000;
}
if(length($format) == 0)
{
print "No format specified, assuming fasta\n";
$format = "fasta";
}

if($overwrite > 0)
{
	print "Overwriting ".$out."\n";
	unlink($out);
}



open(FILE, ">>$out");
print FILE '<?xml version="1.0" standalone="yes"?>'."\n";
print FILE '<beast>'."\n";

print "Reading dates\n";
#create the list of dates
	print FILE '<taxa id="taxa">'."\n";
	open (INFO, $datefile);
	my @lines = <INFO>;
	close (INFO);
	
	for(my $i = 0 ; $i <= $#lines; $i++)
	{
		$lines[$i] =~ s/\r/\n/;
		$lines[$i] =~ s/\n\n/\n/;
		if($lines[$i] =~ m/(.*),(.*),(.*)\n/)
		{
			print FILE '<taxon id="'.$1.'">'."\n".'<date value="'.$2.'" direction="forwards" units="'.$3.'"/>'."\n".'</taxon>'."\n";	
		}	
		elsif($lines[$i] =~ m/(.*),(.*),(.*)/)
		{
			print FILE '<taxon id="'.$1.'">'."\n".'<date value="'.$2.'" direction="forwards" units="'.$3.'"/>'."\n".'</taxon>'."\n";	
		}
		
	}		
	print FILE '</taxa>'."\n";
print "Dates written to file\n";
#read the sequences and write them out as xml
my $seqs;
my $id;

print FILE '<alignment id="alignment" dataType="nucleotide">'."\n";
print "Opening sequence file\n";
my $aln = Bio::SeqIO->new(-file => $sequencefile, -format => $format);
while ( my $gene = $aln->next_seq ) {
	$seqs = $gene->seq;
	$id = $gene->primary_id;
	#$seqs =~ s/N/?/gi;
	print FILE '<sequence>'."\n".'<taxon idref="'.$id.'"/>'."\n\t\t".$seqs."\n".'</sequence>'."\n";
}
print FILE 	'	</alignment>'."\n";

print "Sequence file read and outputted, writing rest of file\n";
	
#specify the pop size bit

#patterns bit
if(length($patterns) > 0)
{
	$patterns =~ s/,/ /g;
		
	print FILE '<mergePatterns id="patterns">'."\n".'<patterns from="1" every="1">'."\n".'<alignment idref="alignment"/>'."\n".'</patterns>'."\n".'<constantPatterns>'."\n".'<alignment idref="alignment"/>'."\n".'<counts>'."\n".'<parameter value=" '.$patterns.' "/>'."\n".'</counts>'."\n".'</constantPatterns>'."\n".'</mergePatterns>';
	
}
else
{
print FILE '		<!-- The unique patterns from 1 to end                                       -->'."\n".'	<patterns id="patterns" from="1">'."\n".'		<alignment idref="alignment"/>'."\n".'	</patterns>'."\n".''."\n";	
}

print FILE '	<!-- A prior assumption that the population size has remained constant       -->'."\n".'	<!-- throughout the time spanned by the genealogy.                           -->'."\n".'	<constantSize id="constant" units="years">'."\n".'		<populationSize>'."\n".'			<parameter id="constant.popSize" value="300.0" lower="0.0" upper="Infinity"/>'."\n".'		</populationSize>'."\n".'	</constantSize>'."\n";

if(length($startingtree) > 0)
{
	
}
else
{ 
#use a UPGMA starting tree
print FILE '<!-- Construct a rough-and-ready UPGMA tree as an starting tree              -->'."\n".'	<upgmaTree id="startingTree">'."\n".'		<distanceMatrix correction="JC">'."\n".'			<patterns>'."\n".''."\n".'				<!-- To generate UPGMA starting tree, only use the 1st aligment, which may be 1 of many aligments using this tree.-->'."\n".'				<alignment idref="alignment"/>'."\n".'			</patterns>'."\n".'		</distanceMatrix>'."\n".'	</upgmaTree>'."\n";
}
#finally add the rest of the special sauce

# for the relaxed;
if($relaxed == 0)
{
print FILE '	<!-- Generate a tree model                                                   -->'."\n".'	<treeModel id="treeModel">'."\n".'		<upgmaTree idref="startingTree"/>'."\n".'		<rootHeight>'."\n".'			<parameter id="treeModel.rootHeight"/>'."\n".'		</rootHeight>'."\n".'		<nodeHeights internalNodes="true">'."\n".'			<parameter id="treeModel.internalNodeHeights"/>'."\n".'		</nodeHeights>'."\n".'		<nodeHeights internalNodes="true" rootNode="true">'."\n".'			<parameter id="treeModel.allInternalNodeHeights"/>'."\n".'		</nodeHeights>'."\n".'	</treeModel>'."\n".''."\n".'	<!-- Generate a coalescent likelihood                                        -->'."\n".'	<coalescentLikelihood id="coalescent">'."\n".'		<model>'."\n".'			<constantSize idref="constant"/>'."\n".'		</model>'."\n".'		<populationTree>'."\n".'			<treeModel idref="treeModel"/>'."\n".'		</populationTree>'."\n".'	</coalescentLikelihood>'."\n".''."\n".'	<!-- The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut, 2006) -->'."\n".'	<discretizedBranchRates id="branchRates">'."\n".'		<treeModel idref="treeModel"/>'."\n".'		<distribution>'."\n".'			<logNormalDistributionModel meanInRealSpace="true">'."\n".'				<mean>'."\n".'					<parameter id="ucld.mean" value="7.5E-4" lower="0.0" upper="1.0"/>'."\n".'				</mean>'."\n".'				<stdev>'."\n".'					<parameter id="ucld.stdev" value="0.3333333333333333" lower="0.0" upper="Infinity"/>'."\n".'				</stdev>'."\n".'			</logNormalDistributionModel>'."\n".'		</distribution>'."\n".'		<rateCategories>'."\n".'			<parameter id="branchRates.categories" dimension="112"/>'."\n".'		</rateCategories>'."\n".'	</discretizedBranchRates>'."\n".'	<rateStatistic id="meanRate" name="meanRate" mode="mean" internal="true" external="true">'."\n".'		<treeModel idref="treeModel"/>'."\n".'		<discretizedBranchRates idref="branchRates"/>'."\n".'	</rateStatistic>'."\n".'	<rateStatistic id="coefficientOfVariation" name="coefficientOfVariation" mode="coefficientOfVariation" internal="true" external="true">'."\n".'		<treeModel idref="treeModel"/>'."\n".'		<discretizedBranchRates idref="branchRates"/>'."\n".'	</rateStatistic>'."\n".'	<rateCovarianceStatistic id="covariance" name="covariance">'."\n".'		<treeModel idref="treeModel"/>'."\n".'		<discretizedBranchRates idref="branchRates"/>'."\n".'	</rateCovarianceStatistic>'."\n".''."\n".'	<!-- The general time reversible (GTR) substitution model                    -->'."\n".'	<gtrModel id="gtr">'."\n".'		<frequencies>'."\n".'			<frequencyModel dataType="nucleotide">'."\n".'				<frequencies>'."\n".'					<parameter id="frequencies" value="0.25 0.25 0.25 0.25"/>'."\n".'				</frequencies>'."\n".'			</frequencyModel>'."\n".'		</frequencies>'."\n".'		<rateAC>'."\n".'			<parameter id="ac" value="1.0" lower="0.0" upper="Infinity"/>'."\n".'		</rateAC>'."\n".'		<rateAG>'."\n".'			<parameter id="ag" value="1.0" lower="0.0" upper="Infinity"/>'."\n".'		</rateAG>'."\n".'		<rateAT>'."\n".'			<parameter id="at" value="1.0" lower="0.0" upper="Infinity"/>'."\n".'		</rateAT>'."\n".'		<rateCG>'."\n".'			<parameter id="cg" value="1.0" lower="0.0" upper="Infinity"/>'."\n".'		</rateCG>'."\n".'		<rateGT>'."\n".'			<parameter id="gt" value="1.0" lower="0.0" upper="Infinity"/>'."\n".'		</rateGT>'."\n".'	</gtrModel>'."\n".''."\n".'	<!-- site model                                                              -->'."\n".'	<siteModel id="siteModel">'."\n".'		<substitutionModel>'."\n".'			<gtrModel idref="gtr"/>'."\n".'		</substitutionModel>'."\n".'		<gammaShape gammaCategories="4">'."\n".'			<parameter id="alpha" value="0.5" lower="0.0" upper="1000.0"/>'."\n".'		</gammaShape>'."\n".'	</siteModel>'."\n".'	<treeLikelihood id="treeLikelihood" useAmbiguities="false">'."\n".'		<patterns idref="patterns"/>'."\n".'		<treeModel idref="treeModel"/>'."\n".'		<siteModel idref="siteModel"/>'."\n".'		<discretizedBranchRates idref="branchRates"/>'."\n".'	</treeLikelihood>'."\n".''."\n".'	<!-- Define operators                                                        -->'."\n".'	<operators id="operators">'."\n".'		<scaleOperator scaleFactor="0.75" weight="0.1">'."\n".'			<parameter idref="ac"/>'."\n".'		</scaleOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="0.1">'."\n".'			<parameter idref="ag"/>'."\n".'		</scaleOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="0.1">'."\n".'			<parameter idref="at"/>'."\n".'		</scaleOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="0.1">'."\n".'			<parameter idref="cg"/>'."\n".'		</scaleOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="0.1">'."\n".'			<parameter idref="gt"/>'."\n".'		</scaleOperator>'."\n".'		<deltaExchange delta="0.01" weight="0.1">'."\n".'			<parameter idref="frequencies"/>'."\n".'		</deltaExchange>'."\n".'		<scaleOperator scaleFactor="0.75" weight="0.1">'."\n".'			<parameter idref="alpha"/>'."\n".'		</scaleOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="3">'."\n".'			<parameter idref="ucld.mean"/>'."\n".'		</scaleOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="3">'."\n".'			<parameter idref="ucld.stdev"/>'."\n".'		</scaleOperator>'."\n".'		<subtreeSlide size="30.0" gaussian="true" weight="15">'."\n".'			<treeModel idref="treeModel"/>'."\n".'		</subtreeSlide>'."\n".'		<narrowExchange weight="15">'."\n".'			<treeModel idref="treeModel"/>'."\n".'		</narrowExchange>'."\n".'		<wideExchange weight="3">'."\n".'			<treeModel idref="treeModel"/>'."\n".'		</wideExchange>'."\n".'		<wilsonBalding weight="3">'."\n".'			<treeModel idref="treeModel"/>'."\n".'		</wilsonBalding>'."\n".'		<scaleOperator scaleFactor="0.75" weight="3">'."\n".'			<parameter idref="treeModel.rootHeight"/>'."\n".'		</scaleOperator>'."\n".'		<uniformOperator weight="30">'."\n".'			<parameter idref="treeModel.internalNodeHeights"/>'."\n".'		</uniformOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="3">'."\n".'			<parameter idref="constant.popSize"/>'."\n".'		</scaleOperator>'."\n".'		<upDownOperator scaleFactor="0.75" weight="3">'."\n".'			<up>'."\n".'				<parameter idref="ucld.mean"/>'."\n".'			</up>'."\n".'			<down>'."\n".'				<parameter idref="treeModel.allInternalNodeHeights"/>'."\n".'			</down>'."\n".'		</upDownOperator>'."\n".'		<swapOperator size="1" weight="10" autoOptimize="false">'."\n".'			<parameter idref="branchRates.categories"/>'."\n".'		</swapOperator>'."\n".'		<randomWalkIntegerOperator windowSize="1" weight="10">'."\n".'			<parameter idref="branchRates.categories"/>'."\n".'		</randomWalkIntegerOperator>'."\n".'		<uniformIntegerOperator weight="10">'."\n".'			<parameter idref="branchRates.categories"/>'."\n".'		</uniformIntegerOperator>'."\n".'	</operators>'."\n".''."\n".'	<!-- Define MCMC                                                             -->'."\n".'	<mcmc id="mcmc" chainLength="'.$chainlength.'" autoOptimize="true">'."\n".'		<posterior id="posterior">'."\n".'			<prior id="prior">'."\n".'				<gammaPrior shape="0.05" scale="10.0" offset="0.0">'."\n".'					<parameter idref="ac"/>'."\n".'				</gammaPrior>'."\n".'				<gammaPrior shape="0.05" scale="20.0" offset="0.0">'."\n".'					<parameter idref="ag"/>'."\n".'				</gammaPrior>'."\n".'				<gammaPrior shape="0.05" scale="10.0" offset="0.0">'."\n".'					<parameter idref="at"/>'."\n".'				</gammaPrior>'."\n".'				<gammaPrior shape="0.05" scale="10.0" offset="0.0">'."\n".'					<parameter idref="cg"/>'."\n".'				</gammaPrior>'."\n".'				<gammaPrior shape="0.05" scale="10.0" offset="0.0">'."\n".'					<parameter idref="gt"/>'."\n".'				</gammaPrior>'."\n".'				<exponentialPrior mean="0.3333333333333333" offset="0.0">'."\n".'					<parameter idref="ucld.stdev"/>'."\n".'				</exponentialPrior>'."\n".'				<oneOnXPrior>'."\n".'					<parameter idref="constant.popSize"/>'."\n".'				</oneOnXPrior>'."\n".'				<coalescentLikelihood idref="coalescent"/>'."\n".'			</prior>'."\n".'			<likelihood id="likelihood">'."\n".'				<treeLikelihood idref="treeLikelihood"/>'."\n".'			</likelihood>'."\n".'		</posterior>'."\n".'		<operators idref="operators"/>'."\n".''."\n".'		<!-- write log to screen                                                     -->'."\n".'		<log id="screenLog" logEvery="1000">'."\n".'			<column label="Posterior" dp="4" width="12">'."\n".'				<posterior idref="posterior"/>'."\n".'			</column>'."\n".'			<column label="Prior" dp="4" width="12">'."\n".'				<prior idref="prior"/>'."\n".'			</column>'."\n".'			<column label="Likelihood" dp="4" width="12">'."\n".'				<likelihood idref="likelihood"/>'."\n".'			</column>'."\n".'			<column label="rootHeight" sf="6" width="12">'."\n".'				<parameter idref="treeModel.rootHeight"/>'."\n".'			</column>'."\n".'			<column label="ucld.mean" sf="6" width="12">'."\n".'				<parameter idref="ucld.mean"/>'."\n".'			</column>'."\n".'		</log>'."\n".''."\n".'		<!-- write log to file                                                       -->'."\n".'		<log id="fileLog" logEvery="'.$logtime.'" fileName="'.$filename.'.log" overwrite="false">'."\n".'			<posterior idref="posterior"/>'."\n".'			<prior idref="prior"/>'."\n".'			<likelihood idref="likelihood"/>'."\n".'			<parameter idref="treeModel.rootHeight"/>'."\n".'			<parameter idref="constant.popSize"/>'."\n".'			<parameter idref="ac"/>'."\n".'			<parameter idref="ag"/>'."\n".'			<parameter idref="at"/>'."\n".'			<parameter idref="cg"/>'."\n".'			<parameter idref="gt"/>'."\n".'			<parameter idref="frequencies"/>'."\n".'			<parameter idref="alpha"/>'."\n".'			<parameter idref="ucld.mean"/>'."\n".'			<parameter idref="ucld.stdev"/>'."\n".'			<rateStatistic idref="meanRate"/>'."\n".'			<rateStatistic idref="coefficientOfVariation"/>'."\n".'			<rateCovarianceStatistic idref="covariance"/>'."\n".'			<treeLikelihood idref="treeLikelihood"/>'."\n".'			<coalescentLikelihood idref="coalescent"/>'."\n".'		</log>'."\n".''."\n".'		<!-- write tree log to file                                                  -->'."\n".'		<logTree id="treeFileLog" logEvery="'.$logtime.'" nexusFormat="true" fileName="'.$filename.'.trees" sortTranslationTable="true">'."\n".'			<treeModel idref="treeModel"/>'."\n".'			<discretizedBranchRates idref="branchRates"/>'."\n".'			<posterior idref="posterior"/>'."\n".'		</logTree>'."\n".'	</mcmc>'."\n".'	<report>'."\n".'		<property name="timer">'."\n".'			<mcmc idref="mcmc"/>'."\n".'		</property>'."\n".'	</report>'."\n".'</beast>';
}
else
{
		print FILE '<!-- Generate a tree model                                                   -->'."\n".'	<treeModel id="treeModel">'."\n".'		<tree idref="startingTree"/>'."\n".'		<rootHeight>'."\n".'			<parameter id="treeModel.rootHeight"/>'."\n".'		</rootHeight>'."\n".'		<nodeHeights internalNodes="true">'."\n".'			<parameter id="treeModel.internalNodeHeights"/>'."\n".'		</nodeHeights>'."\n".'		<nodeHeights internalNodes="true" rootNode="true">'."\n".'			<parameter id="treeModel.allInternalNodeHeights"/>'."\n".'		</nodeHeights>'."\n".'	</treeModel>'."\n".''."\n".'	<!-- Generate a coalescent likelihood                                        -->'."\n".'	<coalescentLikelihood id="coalescent">'."\n".'		<model>'."\n".'			<constantSize idref="constant"/>'."\n".'		</model>'."\n".'		<populationTree>'."\n".'			<treeModel idref="treeModel"/>'."\n".'		</populationTree>'."\n".'	</coalescentLikelihood>'."\n".''."\n".'	<!-- The strict clock (Uniform rates across branches)                        -->'."\n".'	<strictClockBranchRates id="branchRates">'."\n".'		<rate>'."\n".'			<parameter id="clock.rate" value="1.0E-4" lower="0.0" upper="1.0"/>'."\n".'		</rate>'."\n".'	</strictClockBranchRates>'."\n".''."\n".'	<!-- The general time reversible (GTR) substitution model                    -->'."\n".'	<gtrModel id="gtr">'."\n".'		<frequencies>'."\n".'			<frequencyModel dataType="nucleotide">'."\n".'				<frequencies>'."\n".'					<parameter id="frequencies" value="0.25 0.25 0.25 0.25"/>'."\n".'				</frequencies>'."\n".'			</frequencyModel>'."\n".'		</frequencies>'."\n".'		<rateAC>'."\n".'			<parameter id="ac" value="1.0" lower="0.0" upper="Infinity"/>'."\n".'		</rateAC>'."\n".'		<rateAG>'."\n".'			<parameter id="ag" value="1.0" lower="0.0" upper="Infinity"/>'."\n".'		</rateAG>'."\n".'		<rateAT>'."\n".'			<parameter id="at" value="1.0" lower="0.0" upper="Infinity"/>'."\n".'		</rateAT>'."\n".'		<rateCG>'."\n".'			<parameter id="cg" value="1.0" lower="0.0" upper="Infinity"/>'."\n".'		</rateCG>'."\n".'		<rateGT>'."\n".'			<parameter id="gt" value="1.0" lower="0.0" upper="Infinity"/>'."\n".'		</rateGT>'."\n".'	</gtrModel>'."\n".''."\n".'	<!-- site model                                                              -->'."\n".'	<siteModel id="siteModel">'."\n".'		<substitutionModel>'."\n".'			<gtrModel idref="gtr"/>'."\n".'		</substitutionModel>'."\n".'		<gammaShape gammaCategories="4">'."\n".'			<parameter id="alpha" value="0.5" lower="0.0" upper="1000.0"/>'."\n".'		</gammaShape>'."\n".'	</siteModel>'."\n".'	<treeLikelihood id="treeLikelihood" useAmbiguities="false">'."\n".'		<patterns idref="patterns"/>'."\n".'		<treeModel idref="treeModel"/>'."\n".'		<siteModel idref="siteModel"/>'."\n".'		<strictClockBranchRates idref="branchRates"/>'."\n".'	</treeLikelihood>'."\n".''."\n".'	<!-- Define operators                                                        -->'."\n".'	<operators id="operators">'."\n".'		<scaleOperator scaleFactor="0.75" weight="0.1">'."\n".'			<parameter idref="ac"/>'."\n".'		</scaleOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="0.1">'."\n".'			<parameter idref="ag"/>'."\n".'		</scaleOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="0.1">'."\n".'			<parameter idref="at"/>'."\n".'		</scaleOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="0.1">'."\n".'			<parameter idref="cg"/>'."\n".'		</scaleOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="0.1">'."\n".'			<parameter idref="gt"/>'."\n".'		</scaleOperator>'."\n".'		<deltaExchange delta="0.01" weight="0.1">'."\n".'			<parameter idref="frequencies"/>'."\n".'		</deltaExchange>'."\n".'		<scaleOperator scaleFactor="0.75" weight="0.1">'."\n".'			<parameter idref="alpha"/>'."\n".'		</scaleOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="3">'."\n".'			<parameter idref="clock.rate"/>'."\n".'		</scaleOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="3">'."\n".'			<parameter idref="treeModel.rootHeight"/>'."\n".'		</scaleOperator>'."\n".'		<uniformOperator weight="30">'."\n".'			<parameter idref="treeModel.internalNodeHeights"/>'."\n".'		</uniformOperator>'."\n".'		<scaleOperator scaleFactor="0.75" weight="3">'."\n".'			<parameter idref="constant.popSize"/>'."\n".'		</scaleOperator>'."\n".'		<upDownOperator scaleFactor="0.75" weight="3">'."\n".'			<up>'."\n".'				<parameter idref="clock.rate"/>'."\n".'			</up>'."\n".'			<down>'."\n".'				<parameter idref="treeModel.allInternalNodeHeights"/>'."\n".'			</down>'."\n".'		</upDownOperator>'."\n".'	</operators>'."\n".''."\n".'	<!-- Define MCMC                                                             -->'."\n".'	<mcmc id="mcmc" chainLength="'.$chainlength.'" autoOptimize="true">'."\n".'		<posterior id="posterior">'."\n".'			<prior id="prior">'."\n".'				<gammaPrior shape="0.05" scale="10.0" offset="0.0">'."\n".'					<parameter idref="ac"/>'."\n".'				</gammaPrior>'."\n".'				<gammaPrior shape="0.05" scale="20.0" offset="0.0">'."\n".'					<parameter idref="ag"/>'."\n".'				</gammaPrior>'."\n".'				<gammaPrior shape="0.05" scale="10.0" offset="0.0">'."\n".'					<parameter idref="at"/>'."\n".'				</gammaPrior>'."\n".'				<gammaPrior shape="0.05" scale="10.0" offset="0.0">'."\n".'					<parameter idref="cg"/>'."\n".'				</gammaPrior>'."\n".'				<gammaPrior shape="0.05" scale="10.0" offset="0.0">'."\n".'					<parameter idref="gt"/>'."\n".'				</gammaPrior>'."\n".'				<oneOnXPrior>'."\n".'					<parameter idref="constant.popSize"/>'."\n".'				</oneOnXPrior>'."\n".'				<coalescentLikelihood idref="coalescent"/>'."\n".'			</prior>'."\n".'			<likelihood id="likelihood">'."\n".'				<treeLikelihood idref="treeLikelihood"/>'."\n".'			</likelihood>'."\n".'		</posterior>'."\n".'		<operators idref="operators"/>'."\n".''."\n".'		<!-- write log to screen                                                     -->'."\n".'		<log id="screenLog" logEvery="'.$logtime.'">'."\n".'			<column label="Posterior" dp="4" width="12">'."\n".'				<posterior idref="posterior"/>'."\n".'			</column>'."\n".'			<column label="Prior" dp="4" width="12">'."\n".'				<prior idref="prior"/>'."\n".'			</column>'."\n".'			<column label="Likelihood" dp="4" width="12">'."\n".'				<likelihood idref="likelihood"/>'."\n".'			</column>'."\n".'			<column label="rootHeight" sf="6" width="12">'."\n".'				<parameter idref="treeModel.rootHeight"/>'."\n".'			</column>'."\n".'			<column label="clock.rate" sf="6" width="12">'."\n".'				<parameter idref="clock.rate"/>'."\n".'			</column>'."\n".'		</log>'."\n".''."\n".'		<!-- write log to file                                                       -->'."\n".'		<log id="fileLog" logEvery="'.$logtime.'" fileName="'.$filename.'.log" overwrite="false">'."\n".'			<posterior idref="posterior"/>'."\n".'			<prior idref="prior"/>'."\n".'			<likelihood idref="likelihood"/>'."\n".'			<parameter idref="treeModel.rootHeight"/>'."\n".'			<parameter idref="constant.popSize"/>'."\n".'			<parameter idref="ac"/>'."\n".'			<parameter idref="ag"/>'."\n".'			<parameter idref="at"/>'."\n".'			<parameter idref="cg"/>'."\n".'			<parameter idref="gt"/>'."\n".'			<parameter idref="frequencies"/>'."\n".'			<parameter idref="alpha"/>'."\n".'			<parameter idref="clock.rate"/>'."\n".'			<treeLikelihood idref="treeLikelihood"/>'."\n".'			<coalescentLikelihood idref="coalescent"/>'."\n".'		</log>'."\n".''."\n".'		<!-- write tree log to file                                                  -->'."\n".'		<logTree id="treeFileLog" logEvery="'.$logtime.'" nexusFormat="true" fileName="'.$filename.'.trees" sortTranslationTable="true">'."\n".'			<treeModel idref="treeModel"/>'."\n".'			<strictClockBranchRates idref="branchRates"/>'."\n".'			<posterior idref="posterior"/>'."\n".'		</logTree>'."\n".'	</mcmc>'."\n".'	<report>'."\n".'		<property name="timer">'."\n".'			<mcmc idref="mcmc"/>'."\n".'		</property>'."\n".'	</report>'."\n".'</beast>';
}
close (FILE);
print $out." scuccessfuly written\n";
exit 0;


sub help_msg
{

	print "run_BEAST.pl Options; 
	-b beast xml file 
	-d csv file of dates, 3 columns first is taxa id, second is year/day since some point in the past, third is unit of measurement (year/day) 
	-i sequencefile in fasta format 
	-f format of sequence file (default fasta)
	-l interval at which BEAST should log results to file (default, every 1000) 
	-o the filename prefix for BEAST log and tree files
	-c the chain length (default 100,000,000 states)
	-t the starting tree file 
	-k the clock model to use. Choose from strict, relaxedln (lognormal), relaxedexp (exponential) or random
	-p option to specify comma seperated site patterns
	-r delete old BEAST xml files with the same name and replace with a new file\n\n";
	
}
