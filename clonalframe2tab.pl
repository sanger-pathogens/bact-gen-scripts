#! /usr/bin/perl -w
use strict;
use Bio::TreeIO;
use Getopt::Long;

# usage

my $usage = "\nUsage:\nclonalframe2tab.pl -c [clonal frame file] -t [threshold for recombination prediction - default 0.95] -o [output file prefix]\n\n";

sub croak {
	print STDERR "$usage";
	exit(0);
}

# options

my $outfile = "clonalframe";
my $threshold = 0.95;
my $cf_file;
my $help = 0;

GetOptions (
        "c=s" => \$cf_file,
        "t=s" => \$threshold,
        "o=s" => \$outfile,
        "h" => \$help
);

if ($help == 1 || !(defined($cf_file))) {
	croak();
}

# parse CF file

my $in_tree = 0;
my $in_taxa = 0;
my $in_rec = 0;
my $in_ref = 0;

my $tree;
my %taxa;
my %downstream_taxa;
my $n = 1;
my %ref;
my %recstart;
my %recend;
my %refstart;
my %refend;
my %recnode;
my $within_rec = 0;
my $confirmed_rec = 0;
my $reccount = 0;
my $position = 1;
my $firstref = 0;
my $lastref = 0;

open CF, "$cf_file" or croak();

while (<CF>) {
	chomp;
	if (/^\#constree/) {
		$in_tree = 1;
	} elsif (/^\#names/) {
		$in_tree = 0;
		$in_taxa = 1;
	} elsif (/^\#consevents/) {
		$n = 1;
		$in_taxa = 0;
		$in_rec = 1;
	} elsif (/^#consinfo/) {
		if ($within_rec == 1 && $confirmed_rec == 1) {
			$refend{$reccount} = $position-1;
			$recnode{$reccount} = $n-1;
			$within_rec = 0;
			$confirmed_rec = 0;
		} elsif ($within_rec == 1) {
			$reccount--;
		} elsif ($within_rec == 0) {
			$reccount--;
		}
		$in_rec = 0;
	} elsif (/^\#poly/) {
		$in_ref = 1;
		$position = 1;
	} elsif (/^\#ll/) {
		$in_ref = 0;
		#close CF;
	} elsif ($in_tree == 1) {
		$tree = $_;
	} elsif ($in_taxa == 1) {
		$taxa{$n} = $_;
		$n++;
	} elsif ($in_rec == 1) {
		if ($within_rec == 1 && $confirmed_rec == 1) {
			$refend{$reccount} = $position-1;
			$recnode{$reccount} = $n-1;
			$reccount++;
		}
		$within_rec = 0;
		$confirmed_rec = 0;
		my @sites = split(/\s+/,$_);
		my $r = 1;
		$position = 1;
		foreach my $site (@sites) {
			if ($r == 0) {
				$r = 1;
			} else {
				if ($site > $threshold) {
					$confirmed_rec = 1;
				}
				if ($within_rec == 1 && $site <= 0.5) {
					if ($confirmed_rec == 1) {
						$refend{$reccount} = $position-1;
						$recnode{$reccount} = $n;
						$reccount++;
					} else {
						delete($refstart{$reccount});
					}
					$within_rec = 0;
					$confirmed_rec = 0;
				} elsif ($within_rec == 0 && $site > 0.5) {
					$refstart{$reccount} = $position;
					$within_rec = 1;
				}
				$r = 0;
				$position++;
			}
		}
		$n++;
	} elsif ($in_ref == 1) {
		if ($firstref == 0 && $_ > 0) {
			$firstref = $position;
		}
		if ($_ > 0) {
			$lastref = $position;
		}
		$ref{$position} = $_;
		$position++;
	}
}

close CF;

# identify real recombination boundaries

foreach my $reccount (keys %refend) {
	if (defined($ref{$refstart{$reccount}}) && $ref{$refstart{$reccount}} > 0) {
		$recstart{$reccount} = $ref{$refstart{$reccount}};
	} else {
		if ($refstart{$reccount} < $firstref) {
			$refstart{$reccount} = $firstref;
		} else {
			while ($ref{$refstart{$reccount}} == 0 && $refstart{$reccount} > $firstref) {
				$refstart{$reccount}--;
			}
		}
		$recstart{$reccount} = $ref{$refstart{$reccount}};
	}
	if (defined($ref{$refend{$reccount}}) && $ref{$refend{$reccount}} > 0) {
		$recend{$reccount} = $ref{$refend{$reccount}};
	} else {
		if ($refend{$reccount} > $lastref) {
			$refend{$reccount} = $lastref;
		} else {
			while ($ref{$refend{$reccount}} == 0 && $refend{$reccount} < $lastref) {
				$refend{$reccount}++;
			}
		}
		$recend{$reccount} = $ref{$refend{$reccount}};
	}
	
}

undef(%ref);

# process tree

open TMP, "> $outfile.tmptree" or croak();

print TMP "($tree);\n";

close TMP;

my $treeIO = Bio::TreeIO->new(-file => "$outfile.tmptree", -format => "newick");
my $treeObject = $treeIO->next_tree;

# first change taxon names

for my $node ($treeObject->get_nodes) {
	if ($node->is_Leaf) {
		push(@{$downstream_taxa{$node->id}},$taxa{$node->id});
		$node->id($taxa{$node->id});
	}
}

# then find downstream taxa of all internal nodes

for my $node ($treeObject->get_nodes) {
	if ($node->is_Leaf) {
		@{$downstream_taxa{$node->id}} = ($node->id); #edit
	} else {
		unless ($node eq $treeObject->get_root_node) {
			my $nodenum = $node->id;
			foreach my $subnode ($node->get_all_Descendents) {
				if ($subnode->is_Leaf) {
					my $subnode_id = $subnode->id;
					push(@{$downstream_taxa{$nodenum}},$subnode_id);
				}
			}
		}
	}
}

# print out final tree

my $treeoutfile = Bio::TreeIO->new(-file => "> $outfile.tree", format => "Newick");
$treeoutfile->write_tree($treeObject);

system "rm $outfile.tmptree";

# print out tabfile of recombinations

open TAB, "> $outfile.tab" or croak();

foreach $reccount (keys %recstart) {
	print TAB "FT   misc_feature    $recstart{$reccount}..$recend{$reccount}\n";
	if ($#{$downstream_taxa{$recnode{$reccount}}} == -1) {
		print TAB "FT                   /colour=4\n";
	} else {
		print TAB "FT                   /colour=2\n";
	}
	my $nodenum = $recnode{$reccount};
	print TAB "FT                   /node=$nodenum\nFT                   /taxa=\"@{$downstream_taxa{$nodenum}}\"\n";
}

close TAB;
