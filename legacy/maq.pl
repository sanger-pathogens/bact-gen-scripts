#!/usr/bin/perl -w

# Author: lh3

use strict;
use warnings;
use Getopt::Std;
use File::Copy;
use File::Spec;
use Cwd qw/getcwd abs_path/;

my $version = '0.4.7';
&usage if (@ARGV < 1);
# global variables
my @nucl_type = (0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
my @conv = split('', 'XACMGRSVTWYHKDBN');

my $command = shift(@ARGV);
my %func = (easyrun=>\&easyrun, chrpt2snp=>\&chrpt2snp, fastq2bfq=>\&fastq2bfq,
			cat2pair=>\&cat2pair, SNPfilter=>\&SNPfilter, ucsc2snp=>\&ucsc2snp,
			statmap=>\&statmap, sv=>\&sv, demo=>\&demo, splitread=>\&splitread,
			indelpe=>\&indelpe);

die("Unknown command \"$command\".\n") if (!defined($func{$command}));
&{$func{$command}}();
exit(0);

#
# easyrun command
#

sub easyrun
{
	my %opts = (d=>'easyrun', 1=>0, 2=>0, A=>'', n=>2000000, a=>250, e=>3, q=>40, m=>2,
				D=>256, E=>20, N=>2, w=>5, b=>60, B=>2);
	getopts('d:1:2:A:n:a:pe:q:E:D:m:N:b:B:Su:', \%opts);
	die qq(
Usage:   maq.pl easyrun [options] <ref.fa|bfa> <reads1.fq|bfq> [reads2.fastq]

Options: -d DIR    output directory [$opts{d}]
         -n INT    number of reads in a split file [$opts{n}]
         -1 INT    length of read1 [$opts{1}]

         -A FILE   file that contains the 3'adapter [null]
         -m INT    number of mismatches in the seed [$opts{m}]
         -N INT    number of haplotypes [$opts{N}]
         -b INT    maximum sum of qualities of mismatches [$opts{b}]
         -u FILE   create fasta file of unmapped reads [null]

         -e INT    minimum read depth to call a SNP [$opts{e}]
         -D INT    maximum read depth to call a SNP [$opts{D}]
         -E INT    minimum neighbour consensus quality to call a SNP [$opts{E}]
         -q INT    quality threshold for the final SNP calls [$opts{q}]
         -w INT    SNP within +/- INT distant to an indel to be filtered [$opts{w}]
         -B INT    maximum number of SNPs in a window [$opts{B}]
         -S        apply indel detection with the splitread (VERY SLOW)

         -p        mate-pair alignment (PE mode)
         -2 INT    length of read2 (PE only) [$opts{2}]
         -a INT    max insert size (PE only) [$opts{a}]

Notes: For pooled data, please specify the number of haplotypes in the
       pool with `-N' and set `-E 0' at the same time. Do not try
       `-N 1', which only causes segfault.

) if (@ARGV < 2);
	if (defined($opts{p}) && @ARGV != 3) {
	  die("** for PE mode, two reads files should be specified.\n");
	}
	my $d = $opts{d};
	my $cwd = getcwd;
	my $mapopt = " -n $opts{m} -e ".($opts{b}+10);
	map { $_ = File::Spec->rel2abs($_) } @ARGV; # change to absolute path
    $mapopt .= $opts{1}? " -1 $opts{1}" : '';
    $mapopt .= $opts{2}? " -2 $opts{2}" : '';
    $mapopt .= $opts{u}? " -u $opts{u} " : '';
    $mapopt .= $opts{A}? " -d " . File::Spec->rel2abs($opts{A}) : '';
	my $filtopt = "-d $opts{e} -D $opts{D} -n $opts{E}";
	my $exe = gwhich('maq') || die("** Cannot find 'maq' executable.");
	my $pl = gwhich($0) || die("** Cannot find 'maq.pl' script. Should be a bug actually.");
	mkdir($d) unless (-d $d);
	# run fasta2bfa
	my ($is_gzip, $first) = &test_file($ARGV[0]);
	if (ord($first) == ord('>')) { # fasta format
	  if ($is_gzip) { # gzipped
		&run_cmd("gzip -dc $ARGV[0] | $exe fasta2bfa - $d/ref.bfa 2> /dev/null");
	  } else { # plain
		&run_cmd("$exe fasta2bfa $ARGV[0] $d/ref.bfa 2> /dev/null");
	  }
	} else { # create a symbolic link
	  &run_cmd("ln -sf $ARGV[0] $d/ref.bfa");
	}
	# run fastq2bfq && match
	shift(@ARGV);
	my $i = 0;
	my (@map_files, $map_file);
	foreach my $file (@ARGV) {
	  ++$i;
	  ($is_gzip, $first) = &test_file($file);
	  if (ord($first) == ord('@')) { # fastq format
		if ($is_gzip) { # gzipped
		  &run_cmd("gzip -dc $file | $exe fastq2bfq -n $opts{n} - $d/read$i");
		} else { # plain fastq
		  &run_cmd("$exe fastq2bfq -n $opts{n} $file $d/read$i");
		}
	  } elsif ($is_gzip) { # probably bfq
		&run_cmd("ln -sf $file $d/read$i.bfq");
	  } else {
		die("** Cannot guess the format of file '$file'.");
	  }
	}
	# run alignment
	if (defined $opts{p}) { # paired reads
	  my @bfq_files = `(cd $d; find . -name "read1*.bfq")`;
	  foreach (@bfq_files) {
		chomp;
		/^\.\/read1(.*)\.bfq$/;
		my $tag = $1;
		my $bfqs = (-f "$d/read2$tag.bfq")? "read1$tag.bfq read2$tag.bfq" : "read1$tag.bfq";
		&run_cmd("(cd $d; $exe map $mapopt aln$tag.map ref.bfa $bfqs 2> aln$tag.map.log)");
		push(@map_files, "aln$tag.map");
	  }
	} else { # unpaired reads
	  my @bfq_files = `(cd $d; find . -name "read*.bfq")`;
	  foreach (@bfq_files) {
		chomp;
		/^\.\/read(.*)\.bfq$/;
		my $tag = $1;
		&run_cmd("(cd $d; $exe map $mapopt -u unmap$tag.txt aln$tag.map ref.bfa read$tag.bfq 2> aln$tag.map.log)");
		push(@map_files, "aln$tag.map");
	  }
	}
	if (@map_files > 1) { # we should merge them
	  &run_cmd("(cd $d; $exe mapmerge all.map @map_files)");
	} else {
	  &run_cmd("(cd $d; mv $map_files[0] all.map)");
	}
	$map_file = "all.map";
	&run_cmd("(cd $d; $exe mapcheck ref.bfa $map_file > mapcheck.txt)");
	&run_cmd("(cd $d; $exe assemble -N $opts{N} -Q $opts{b} consensus.cns ref.bfa $map_file 2> assemble.log)");
	&run_cmd("$exe cns2fq $d/consensus.cns > $d/cns.fq");
	&run_cmd("$exe cns2snp $d/consensus.cns > $d/cns.snp");
	&run_cmd("$exe cns2win $d/consensus.cns > $d/cns.win");
	&run_cmd("$exe indelsoa $d/ref.bfa $d/$map_file > $d/cns.indelse");
	if (defined $opts{p}) {
	  &run_cmd("$exe indelpe $d/ref.bfa $d/$map_file > $d/cns.indelpe");
	  &run_cmd("$pl SNPfilter -q $opts{q} -w $opts{w} -N $opts{B} -f $d/cns.indelse -F $d/cns.indelpe $filtopt -Q60 $d/cns.snp > $d/cns.final.snp");
	} else {
	  &run_cmd("(cd $d; touch unmap.indel)");
	  if (!defined($opts{S})) {
		&run_cmd("$pl SNPfilter -q $opts{q} -w $opts{w} -N $opts{B} -f $d/cns.indelse $filtopt $d/cns.snp > $d/cns.final.snp");
	  } else { # call splitread indels
		&run_cmd("(cd $d; cat unmap*.txt | $pl splitread -p unmap - ref.bfa 2> unmap.log)");
		&run_cmd("$pl SNPfilter -q $opts{q} -w $opts{w} -N $opts{B} -f $d/cns.indelse -S $d/unmap.indel $filtopt $d/cns.snp > $d/cns.final.snp");
	  }
	}
	&run_cmd("(cd $d; ln -s cns.final.snp cns.filter.snp)");
#	my $maq_plot = gwhich("maq_plot.pl");
#	my $gnuplot = gwhich("gnuplot");
#	if ($maq_plot && $gnuplot) {
#		&run_cmd("$maq_plot depth -x 0.2 $d/depth $d/cns.win");
#	}
	&run_cmd("$pl statmap $d/*.map.log");
}

sub statmap {
  my ($is_paired, $n_reads1, $n_reads2, $n_moved_high, $n_moved_low, $n_mapped1, $n_mapped2, $n_paired, $n_added);
  $is_paired = $n_reads1 = $n_reads2 = $n_moved_high = $n_moved_low = $n_mapped1 = $n_mapped2 = $n_paired = $n_added = 0;

  while (<>) {
	if (/^--\s*\(total,\s*isPE,.*\) .*=\s*\((\d+),\s*(\d+),\s*(\d+),\s*(\d+)\)/) {
	  $is_paired = $2;
	  if ($is_paired) { $n_reads2 += $1; $n_mapped2 += $3; }
	  else { $n_reads1 += $1; $n_mapped1 += $3; }
	  $n_paired += $4;
	} elsif (/match_data2mapping.* (\d+) pairs are added/) {
	  $n_added += $1;
	} elsif (/match_data2mapping.* first.*\((\d+), (\d+)\).*second.*\((\d+), (\d+)\)/) {
	  $n_moved_high += $1+$3;
	  $n_moved_low += $2+$4;
	}
  }
  my $tot = $n_reads1 + $n_reads2;
  my $map1_ratio = ($n_reads1 == 0)? 'NA' : int(10000 * $n_mapped1 / $n_reads1) / 100;
  my $map2_ratio = ($n_reads2 == 0)? 'NA' : int(10000 * $n_mapped2 / $n_reads2) / 100;
  my $pair_ratio = ($n_mapped2 == 0)? 'NA' : int(10000 * $n_paired / $n_mapped2) / 100;
  my $add_ratio = ($n_paired == 0)? 'NA' : int(10000 * $n_added / $n_paired) / 100;
  my $high_ratio = ($n_paired == 0)? 'NA' : int(10000 * $n_moved_high / $n_paired) / 100;
  my $low_ratio = ($n_paired == 0)? 'NA' : int(10000 * $n_moved_low / $n_paired) / 100;
  print qq(
-- == statmap report ==\n
-- # single end (SE) reads: $n_reads1
-- # mapped SE reads: $n_mapped1 (/ $n_reads1 = $map1_ratio%)
-- # paired end (PE) reads: $n_reads2
-- # mapped PE reads: $n_mapped2 (/ $n_reads2 = $map2_ratio%)
-- # reads that are mapped in pairs: $n_paired (/ $n_mapped2 = $pair_ratio%)
-- # Q>=30 reads that are moved to meet mate-pair requirement: $n_moved_high (/ $n_paired = $high_ratio%)
-- # Q<30 reads that are moved to meet mate-pair requirement: $n_moved_low ($low_ratio%)
\n);
}

#
# SNPfilter command
#

sub SNPfilter {
  my %opts = (f=>'', S=>'', s=>3, m=>1, Q=>40, d=>3, F=>'', w=>5, D=>256, N=>2, W=>10, n=>20, c=>-1.0, q=>20);
  getopts('af:s:m:Q:d:D:w:F:W:N:c:n:S:q:', \%opts);
  die(qq{
Usage:   maq.pl SNPfilter [options] <cns2snp.snp>

Options: -d INT        minimum depth to call a SNP [$opts{d}]
         -D INT        maximum depth (<=254), otherwise ignored [$opts{D}]
         -n INT        minimum neighbouring quality [$opts{n}]
         -Q INT        required max mapping quality of the reads covering the SNP [$opts{Q}]
         -q INT        minimum consensus quality [$opts{q}]
         -w INT        size of the window in which SNPs should be filtered out [$opts{w}]
         -S FILE       splitread output [null]
         -F FILE       indelpe output [null]
         -f FILE       indelsoa output [null]
         -s INT        indelsoa score (= left_clip + right_clip - across) [$opts{s}]
         -m INT        indelsoa: max number of reads mapped across the indel [$opts{m}]
         -W INT        window size for filtering dense SNPs [$opts{W}]
         -N INT        maximum number of SNPs in a window [$opts{N}]
         -a            alternative filter for single end reads
\n}) unless (@ARGV);
  my (%hash, $fh);
  my $skip = $opts{w};
  if ($opts{f}) { # for indelsoa
	my $n = 0;
	open($fh, $opts{f}) || die;
	while (<$fh>) {
	  my @t = split;
	  next unless ($t[4]+$t[5]-$t[3] >= $opts{s} && $t[3] <= $opts{m}); # a simple filter
	  ++$n;
	  if ($t[2] < 0) { # potential deletion
		for (my $x = $t[1] + $t[2] - $skip; $x <= $t[1] + $skip; ++$x) {
		  $hash{$t[0],$x} = 1;
		}
	  } else { # potential insertion
		for (my $x = $t[1] - $skip; $x <= $t[1] + $t[2] + $skip; ++$x) {
		  $hash{$t[0],$x} = 1;
		}
	  }
	}
	close($fh);
	warn("-- $n potential soa-indels pass the filter.\n");
  }
  if ($opts{F}) { # for indelpe
	my $n = 0;
	open($fh, $opts{F}) || die;
	while (<$fh>) {
	  my @t = split;
	  next unless ($t[2] eq '*' || $t[2] eq '+');
	  ++$n;
	  for (my $x = $t[1] - 1; $x < $t[1] + 2*$skip; ++$x) {
		$hash{$t[0],$x} = 1;
	  }
	}
	close($fh);
	warn("-- $n potential pe-indels pass the filter.\n");
  }
  if ($opts{S}) { # for splitread
	my $n = 0;
	open($fh, $opts{S}) || die;
	while (<$fh>) {
	  next unless (/^VR/);
	  my @t = split;
	  next if ($t[5] < 2 || ($t[5] < 4 && $t[4] == 0));
	  ++$n;
	  $hash{$t[1],$_} = 1 for ($t[2] .. $t[3]);
	}
	close($fh);
	warn("-- $n potential splitread-indels pass the filter.\n");
  }
  my $is_alter = defined($opts{a});
  my (@last, $last_chr);
  $last_chr = '';
  while (<>) {
	my @t = split;
	next if ($hash{$t[0],$t[1]});
	my $is_good;
	if (!$is_alter) { # the default filter
	  $is_good = ($t[5] >= $opts{d} && $t[5] <= $opts{D} && $t[6] > $opts{c} && $t[7] >= $opts{Q} && $t[8] >= $opts{n})? 1 : 0;
	} else { # the alternative filter for SE reads
	  $is_good = ($t[5] >= $opts{d} && $t[5] <= $opts{D} && $t[6] > $opts{c} && $t[6] <= 4.0 && $t[8] >= $opts{n})? 1 : 0;
	}
	$is_good = 0 unless ($t[4] >= $opts{q} || ($t[2] ne $t[9] && $t[4]+$t[10] >= $opts{q})); # consensus quality filter
	next unless ($is_good); # drop
	if ($t[0] ne $last_chr) { # a different chr, print
	  map { print $_->{L} if ($_->{F}) } @last;
	  @last = ();
	  $last_chr = $t[0];
	}
	# The following block provided by Nathans Weeks. Thanks, Nathans.
	push(@last, {L => $_, X => $t[1], F => 1}); # Enqueue current SNP
	if ($#last == $opts{N}) {                   # number of SNPs in queue is N+1
	  if ($last[$#last]{X} - $last[0]{X} < $opts{W}) { # if all within window W
		map {$_->{F} = 0} @last; # all SNPs in the window of size W are "bad"
	  }
	  print STDOUT $last[0]{L} if ($last[0]{F}); # print first SNP if good
	  shift @last # dequeue first SNP
	}
  }
  # print the last few lines if applicable
  map { print $_->{L} if ($_->{F}) } @last;
}

#
# splitread
#

sub splitread {
  my %opts = (1=>12, 2=>12, p=>'splitread', a=>20);
  getopts('1:2:p:a:', \%opts);
  die(qq/
Usage:   maq.pl splitread [options] <input.unmap> <in.bfa>

Options: -1 INT     length of read1 [$opts{1}]
         -2 INT     length of read2 [$opts{2}]
         -a INT     maximum indel size [$opts{a}]
         -p STR     prefix [$opts{p}]

/) if (@ARGV < 2);
  my $maq = gwhich('maq') || die("[splitread] fail to locate 'maq'. Abort!\n");
  my ($l1, $l2) = ($opts{1}, $opts{2});
  my ($n, $fh, $fh1, $fh2);
  # generate fastq
  warn("-- generate paired end reads\n");
  open($fh, ($ARGV[0]=~/\.gz$/)? "gzip -dc $ARGV[0]" : $ARGV[0]) || die;
  open($fh1, ">$opts{p}1.fq") || die;
  open($fh2, ">$opts{p}2.fq") || die;
  while (<$fh>) {
	my @t = split;
	next if ($t[1] != 99);
	my $l = length($t[2]);
	die("[splitread] reads are too short!\n") if ($l < $l1 + $l2);
	# sequence
	my ($s1, $s2) = (substr($t[2], 0, $l1), substr($t[2], -$l2));
	# discard poor reads
	next if ($s1 =~ /N/i || $s2 =~ /N/i);
	$n = ($s1 =~ tr/A/A/);
	next if ($n == $l1);
	$n = ($s2 =~ tr/A/A/);
	next if ($n == $l2);
	# reverse $s2
	$s2 = reverse($s2);
	$s2 =~ tr/ACGTN/TGCAN/;
	# quality
	my ($q1, $q2) = (substr($t[3], 0, $l1), substr($t[3], -$l2));
	$q2 = reverse($q2);
	# output
	print $fh1 "\@$t[0]_$l/1\n$s1\n+\n$q1\n";
	print $fh2 "\@$t[0]_$l/2\n$s2\n+\n$q2\n";
  }
  close($fh); close($fh1); close($fh2);
  # paired end alignment
  my $ins_size = $l1 + $l2 + $opts{a};
  &run_cmd("$maq fastq2bfq $opts{p}1.fq $opts{p}1.bfq");
  &run_cmd("$maq fastq2bfq $opts{p}2.fq $opts{p}2.bfq");
  &run_cmd("$maq map -Pn1 -a $ins_size $opts{p}.map $ARGV[1] $opts{p}1.bfq $opts{p}2.bfq");
  # get indel sites
  my @last = ('', 0, 0, 0, 0); # name, begin, end, 
  open($fh, "$maq mapview $opts{p}.map |") || die;
  open($fh1, ">$opts{p}.indel") || die;
  print $fh1 qq{CC\tCC comments
CC\tRD read_id chr begin end ins_size
CC\tVR chr begin end ins_size n_reads
CC
CC\tIn VR lines, zero ins_size indicates that there are indels of different
CC\tlengths. In addition, indels supported by one read is not reliable.
CC
};
  while (<$fh>) {
	my @t = split;
	next if ($t[5] != 18 || $t[8] == 0 || $t[4] < 0);
	$t[0] =~ /^(\S+)_(\d+)\/[12]$/;
	my ($name, $l) = ($1, $2);
	next if ($t[4] == $l); # no indel
	if ($last[0] ne $t[1]) {
	  if ($last[0]) {
		print $fh1 join("\t", "VR", @last), "\n";
	  }
	  @last = ($t[1], 0, 0, 0, 0);
	}
	my ($b, $e, $s) = ($t[2], $t[2]+$t[4]-1, $l - $t[4]);
	if ($last[1] == 0) { # the first
	  @last = ($t[1], $b, $e, $s, 1);
	} elsif ($b <= $last[2]) { # overlap
	  if ($s != $last[3]) { # conflict
		$last[3] = 0;
	  }
	  ++$last[4];
	  $last[1] = $b;
	} else { # non-overlap
	  print $fh1 join("\t", "VR", @last), "\n";
	  @last = ($t[1], $b, $e, $s, 1);
	}
	print $fh1 join("\t", 'RD', $name, $t[1], $b, $e, $s), "\n";
  }
  close($fh); close($fh1);
}

#
# for Sanger's PE read format only
#

sub cat2pair
{
  my %opts = (1=>0, 2=>0);
  getopts('1:2:', \%opts);
  die qq(Usage: maq.pl cat2pair [-1 $opts{1}] [-2 $opts{2}] <read1_len> <input.fastq> [output.fastq]\n) if (@ARGV < 2);
  my ($tl1, $tl2) = ($opts{1}, $opts{2});
  my $size1 = $ARGV[0];
  my $fn = $ARGV[1];
  my $fn_out = (@ARGV >= 3)? $ARGV[2] : $ARGV[1];
  $fn_out =~ s/^.*\/([^\/\s]+)$/$1/ if ($fn_out =~ /\//);
  $fn_out =~ s/\.gz$// if ($fn_out =~ /\.gz$/);
  mkdir("read1"); mkdir("read2");
  my ($fh1, $fh2, $fhin);
  $fn = "gzip -dc $fn |" if ($fn =~ /\.gz$/);
  open($fhin, $fn) || die;
  open($fh1, ">read1/$fn_out") || die;
  open($fh2, ">read2/$fn_out") || die;
  while (<$fhin>) {
	if (/^@/) {
	  chomp; print $fh1 "$_/1\n"; print $fh2 "$_/2\n";
	  $_ = <$fhin>; chomp; # sequence
	  print $fh1 substr($_, $tl1, $size1-$tl1), "\n+\n";
	  print $fh2 substr($_, $size1+$tl2), "\n+\n";
	  <$fhin>; $_ = <$fhin>; chomp; # qualities
	  print $fh1 substr($_, $tl1, $size1-$tl1), "\n";
	  print $fh2 substr($_, $size1+$tl2), "\n";
	}
  }
  close($fhin); close($fh1); close($fh2);
}

#
# for Sanger's farm only
#

sub fastq2bfq {
  my %opts = (s=>1000000, e=>'fastq', d=>'', r=>'', a=>250);
  getopts('s:e:Ed:r:a:', \%opts);
  die qq(
Usage:   maq.pl fastq2bfq [-s nreads] [-e ext] <src_dir> <dst_dir>\n
Options: -s INT       number of reads per file [1000000]
         -e STR       extension of the read files [fastq]
         -d FILE      adapter sequence file [null]
         -r FILE      reference genome [hg18_male.bfa]
         -a INT       maximum insert size [$opts{a}]
         -E           the input in is Solexa's "export" format
\n) if (@ARGV < 2);
  my ($size, $ext) = ($opts{s}, $opts{e});
  my $is_export = (defined $opts{E})? 1 : 0;
  my ($src_dir, $dst_dir) = ($ARGV[0], $ARGV[1]);
  $dst_dir =~ s/\/$//;
  my ($fh_lst, $fh_pl);
  open($fh_lst, ">$dst_dir.lst") || die;
  open($fh_pl, ">$dst_dir.pl") || die;
  my $tmp_dir = "$dst_dir/tmp";
  my $cwd = getcwd;
  die("FATAL ERROR: source directory '$src_dir' does not exist!\n") unless (-d $src_dir);
  mkdir($dst_dir) unless (-d $dst_dir);
  mkdir($tmp_dir) unless (-d $tmp_dir);
  chdir($src_dir);
  my @list = `(find . -name "*.$ext" -follow; find . -name "*.$ext.gz" -follow)`; # get the list of files
  chdir($cwd);
  $size = $size >> 1 << 1;
  # run 'maq fastq2bfq'
  my $maq = gwhich("maq") || die("ERROR: Cannot find maq\n");
  my $fq_all2std = gwhich("fq_all2std.pl") || die("ERROR: Cannot find fq_all2std.pl\n");
  foreach (@list) {
	chomp; s/^\.\///;
	my $ori = $_;
	my $prog = (/\.gz$/i? 'zcat' : 'cat') . " $src_dir/$ori";
	if ($is_export) {
	  next if ($ori !~ /export/i);
	  $prog = "$prog | $fq_all2std export2sol | $maq sol2sanger - -";
	}
	s/\//-/g; s/\.gz$//i; s/\.$ext$//;
	my ($cur_size, $name);
	if ($is_export) {
	  # IMPORTANT: I do not know what single end reads look like.
	  #   possibly this part does not work well for SE reads.
	  if (/s_\d+_([12])_/) { # paired
		$cur_size = int($size/2);
		$name = "reads$1";
		s/(s_\d+_)[12]_/$1/;
	  }
	} else {
	  if (/read(s?)[12]/) { # paired
		$cur_size = int($size/2);
		$name = (/read(s?)1/i)? "reads1" : "reads2";
		s/read(s?)[12]//i;
	  } else {
		$cur_size = $size;
		$name = "reads.bfq";
	  }
	}
	s/--/-/g; s/^-//;
	$name = (($name =~ /[12]/)? "$_-PE" : "$_-SE") . ":$name";
	&run_cmd("$prog | $maq fastq2bfq -n $cur_size - $tmp_dir/$name");
  }
  # move files to separate directories
  my $n_jobs = 0;
  @list = `(cd $tmp_dir; ls)`;
  foreach (@list) {
	chomp;
	if (/(.*):(reads[12]?)\@(\d+)\.bfq$/) {
	  my $d = "$dst_dir/$1\@$3";
	  my $t = "$d/$2.bfq";
	  unless (/reads2\@\d+\.bfq/) {
		++$n_jobs;
		print $fh_lst "$d\n";
		mkdir($d);
	  }
	  move("$tmp_dir/$_", $t);
	}
  }
  system("rm -fr $tmp_dir");
  close($fh_lst);
  # generate configuration file for farm-run.pl
  $opts{d} = "-d ".abs_path($opts{d}) if ($opts{d});
  if ($opts{r}) {
	$opts{r} = abs_path($opts{r});
  } else {
	$opts{r} = '../../human_male.bfa';
  }
  print $fh_pl qq(\%fr_config =
(
  # The list of directories that should be processed. Do NOT change this.
  run_list=>'$dst_dir.lst',
  # number of LSF jobs at a time
  n_jobs=>$n_jobs,
  LSF_queue=>'long',
  LSF_resource=>'select[mem>800 && type==X86_64] rusage[mem=800]',
  # dir of the executables
  binary_path=>"$ENV{HOME}/lsf-prog",
  # dir of the scripts, if it is different from binary_path
  script_path=>"$ENV{HOME}/lsf-prog/scripts",
  # in each working directory, this function will be called.
  action=>\\&func
);

sub func
{
  my \$dir = shift; # the working directory
  # The current working dir is the one where *.bfq are staying.
  # usually the reference sequence is put two-level higher than
  # the current working dir.
  my \$maq_comm = "maq map -a $opts{a} $opts{d} $dst_dir.map $opts{r}";
  if (-f "reads.bfq") { # single end
    system("\$maq_comm reads.bfq");
  } elsif (-f "reads1.bfq" && ! -f "reads2.bfq") { # single end
    system("\$maq_comm reads1.bfq");
  } elsif (-f "reads1.bfq" && -f "reads2.bfq") { # paired end
    system("\$maq_comm reads1.bfq reads2.bfq");
  }
}
);
  close($fh_pl);
}

sub sv {
  my %opts = (i=>150, l=>35, q=>35, s=>7);
  getopts('s:i:l:q:', \%opts);
  die("
Usage:   maq.pl sv <in.mapview>\n
Options: -i INT    maximum insert size [$opts{i}]
         -l INT    average read length [$opts{l}]
         -q INT    minimum alternative mapping quality [$opts{q}]
         -s INT    minimum length of a region [$opts{s}]\n
") unless (@ARGV);
  my $d = $opts{i} - $opts{l};
  my ($begins, $beginc, $lasts, $lastc) = ('', -1, '', -1);
  my (@regs, %read, @reg_name);
  my @reg_seq;

  warn("-- read the mapview output\n");
  while (<>) {
	my @t = split;
	next if ($t[8] <= $opts{q});
	next if ($t[5] == 18 || $t[5] == 64 || $t[5] == 130);
	my $do_break = ($t[1] ne $lasts || $t[2] - $lastc > $d)? 1 : 0;
	if ($do_break) {
	  if ($lastc - $beginc > $opts{s}) { # skip short/unreliable regions
		my $k = @regs;
		my $flag = ($lastc - $beginc < $opts{i})? '*' : '.';
		push(@reg_name, "$begins\t$beginc\t$lastc\t$flag");
		my $p = \@{$regs[$k]};
		foreach (@reg_seq) {
		  push(@$p, $_);
		  my @s = split;
		  push(@{$read{$s[0]}}, $k);
		}
	  }
	  ($begins, $beginc) = @t[1..2];
	  @reg_seq = ();
	}
	$t[0] =~ s/\/[12]$//;
	push(@reg_seq, join(" ", @t[0..3,5,8,13]));
	($lasts, $lastc) = @t[1..2];
  }

  # build connections
  warn("-- link regions\n");
  my %link;
  foreach my $x (keys %read) {
	my $p = $read{$x};
	next if (@$p != 2);
	my $key = sprintf("%.10d %.10d", $p->[0], $p->[1]); # $p->[0] < $p->[1] always stands
	if (defined $link{$key}) {
	  ++$link{$key};
	} else {
	  $link{$key} = 1;
	}
  }

  # print connected regions
  warn("-- print result\n");
  foreach (reverse sort{$link{$a}<=>$link{$b}} keys %link) {
	my @s = split;
	my $x = $_;
	my @count;
	$count[$_] = 0 for (0..3);
	foreach my $y (@{$regs[$s[0]]}) {
	  my @t = split(" ", $y);
	  ++$count[($t[3] eq '+')? 0 : 1];
	}
	foreach my $y (@{$regs[$s[1]]}) {
	  my @t = split(" ", $y);
	  ++$count[($t[3] eq '+')? 2 : 3];
	}
	# infer flag
	my $flag = ($s[0] == $s[1])? 'LOP' : '';
	if (!$flag) {
	  my @t1 = split(/\s+/, $reg_name[$s[0]]);
	  my @t2 = split(/\s+/, $reg_name[$s[1]]);
	  if ($t1[0] ne $t2[0]) {
		$flag = 'DIF';
	  } elsif (($count[1] == 0 && 2*$count[2] < $count[3])
			   || ($count[2] == 0 && 2*$count[1] < $count[0])) {
		$flag = 'DEL';
	  } else {
		$flag = 'AMB';
	  }
	}
	print "$flag\t$link{$x}\t$reg_name[$s[0]]\t$count[0]\t$count[1]\t$reg_name[$s[1]]\t$count[2]\t$count[3]\n";
  }
}

#
# chrpt2snp command
#

sub chrpt2snp
{
	# calculate @mm
	my @mm;
	die("Usage: maq.pl chrpt2snp <chr_rpt.dbSNP>\n") unless (@ARGV);
	for (my $i = 0; $i != 0x10000; ++$i) {
		my ($x, $k) = (1, 0);
		for (my $l = 0; $l != 16; ++$l, $x <<= 1) {
			++$k if ($i & $x);
		}
		$mm[$i] = $k;
	}
	my $fh;
	open($fh, "| sort -k1,1 -k2,2n");
	while (<>) {
		next unless (/^\d+/ && /reference$/);
		my @t = split("\t", $_);
		next unless ($t[11] =~ /^\d+/); # no position
		my $vs = $mm[$t[16]&0xffff] + $mm[$t[16]>>16&0xffff];
		printf $fh ("chr$t[6]\t$t[11]\tN\tN\t%d\t%.3f\n", $vs*10, ($t[13] ne ' ')? $t[13] : 0.0);
	}
	close($fh);
}

#
# ucsc2snp command
#

sub ucsc2snp {
  warn("-- only one-basepair substitutions and indels will be retained.\n");
  my %conv = (AC=>'M', CA=>'M', AG=>'R', GA=>'R',
			  AT=>'W', TA=>'W', CG=>'S', GC=>'S',
			  CT=>'Y', TC=>'Y', GT=>'K', TG=>'K');
  my ($n_not2, $n_not1) = (0, 0);
  while (<>) {
	my @t = split;
	if ($t[9] !~ /-/) {
	  if (length($t[9]) != 3 || length($t[7]) != 1) {
		++$n_not1;
		next;
	  }
	}
	my ($a1, $a2) = split("/", $t[9]);
	my $SNP;
	if ($t[11] eq 'single') {
	  $SNP = $conv{"$a1$a2"};
	} else {
	  if ($t[7] eq '-') {
		$SNP = ($a1 ne '-')? $a1 : $a2;
	  } else {
		$SNP = '-';
	  }
	}
	my $qual = 0;
	if ($t[12] ne 'unknown') {
	  my @s = split(",", $t[12]);
	  $qual = @s * 10;
	}
	my $ref = (length($t[7]) == 1)? $t[7] : 'N';
	print "$t[1]\t$t[3]\t$ref\t$SNP\t$qual\t$t[13]\t$t[15]\n";
  }
  warn("-- discarded: ($n_not1, $n_not2)\n");
}

sub demo {
  my %opts = (N=>1000000, d=>'maqdemo');
  getopts('N:d:hs', \%opts);
  die("
Usage:   maq.pl demo [-N npairs] [-d outdir] [-h] <in.fa> <in.simudat>\n
Options: -N INT    number of read pairs [$opts{N}]
         -d DIR    output directory [$opts{d}]
         -s        single-end mode in alignment
         -h        haploid mode in simulation\n
") if (@ARGV < 2);
  my $peopt = (defined $opts{s})? '' : '-p';
  my $simuopt = "-N $opts{N}";
  $simuopt .= " -h" if (defined $opts{h});
  my $maq = gwhich("maq");
  my $maq_pl = gwhich("maq.pl");
  my $eval_pl = gwhich("maq_eval.pl");
  die("** 'maq', 'maq.pl' and 'maq_eval.pl' MUST be on the \$PATH\n") unless ($maq && $maq_pl && $eval_pl);
  &run_cmd("mkdir -p $opts{d}");
  &run_cmd("$maq simulate $simuopt $opts{d}/r1.fq $opts{d}/r2.fq $ARGV[0] $ARGV[1] > $opts{d}/true.snp");
  &run_cmd("$maq fasta2bfa $ARGV[0] $opts{d}/ref.bfa");
  &run_cmd("(cd $opts{d}; $maq_pl easyrun $peopt -d easyrun ref.bfa r1.fq r2.fq)");
  &run_cmd("(cd $opts{d}; $maq simustat easyrun/all.map > eval.simustat)");
  &run_cmd("(cd $opts{d}; $eval_pl sub -p eval.sub true.snp true.snp easyrun/cns.filter.snp)");
  &run_cmd("(cd $opts{d}; $eval_pl indelpe true.snp easyrun/cns.indelpe > eval.indelpe)") if ($peopt);
  &run_cmd("(cd $opts{d}; $eval_pl indelsoa true.snp easyrun/cns.indelse > eval.indelse)");
  warn("++ $opts{d}/easyrun/cns.filter.snp gives the SNPs that passes most of filters.\n");
  warn("++ $opts{d}/easyrun/cns.final.snp is the high-quality subset of the previous SNPs.\n");
  warn("++ $opts{d}/easyrun/mapcheck.txt gives some statistics about qualities.\n");
  warn("++ $opts{d}/true.snp contains all variants used in simulation.\n");
  warn("++ $opts{d}/eval.* give various benchmarks. Their formats will be documented later.\n");
}

#
#
#

sub indelpe {
  my %opts = (m=>6, n=>4);
  getopts('m:n:', \%opts);
  die("Usage: maq.pl indelpe [-m $opts{m}] [-n $opts{n}] <in.indelpe>\n") if (@ARGV == 0 && -t STDIN);
  my ($om, $on) = ($opts{m}, $opts{n});
  while (<>) {
	my @t = split;
	# calculate length of homopolymer tract
	$_ = $t[8];
	my $c = substr($_, 0, 1);
	/^($c+)/;
	my $n_hpt = length($1);
	$n_hpt = length($1) if (/^(${c}{3,}[^$c\s]{1,2}${c}{3,})/);
	$n_hpt = length($1) if (/^(${c}{4,}[^$c\s]{1,4}${c}{4,})/);
	# get rid of reads afftected by homopolymer tract
	my (@b, @e, @mm, @nb, @ne, @nmm);
	$_ = $t[10]; s/^@//; @b = split(',', $_);
	$_ = $t[11]; s/^@//; @e = split(',', $_);
	$_ = $t[12]; s/^@//; @mm = split(',', $_);
	for my $i (0 .. $#b) {
	  if ($b[$i] > $om && $e[$i] > $om && $e[$i] > $n_hpt + $on && $mm[$i] <= 2) {
		push(@nb, $b[$i]);
		push(@ne, $e[$i]);
		push(@nmm, $mm[$i]);
	  }
	}
	$t[3] = $t[5] + $t[6] + @nb;
	$t[9] = @nb;
	$t[10] = "@".join(",", @nb).",";
	$t[11] = "@".join(",", @ne).",";
	$t[12] = "@".join(",", @nmm).",";
	print join("\t", @t), "\n";
  }
}

#
# Usage
#

sub usage
{
	die qq(
Program: maq.pl (helper script for maq)
Version: $version
Contact: Heng Li <lh3\@sanger.ac.uk>

Usage:   maq.pl <command> [options] <input> [...]

Command: easyrun      simple pipeline for small dataset
         demo         demonstration of maq functionalities (for SE only)
         SNPfilter    filter SNPs
         statmap      extract statistics from the error output of 'maq map'
         chrpt2snp    convert dbSNP's chr_rpt file to .snp file
         ucsc2snp     convert UCSC's SNP dump to .snp file
         fastq2bfq    convert fastq in batch
         cat2pair     convert paired end reads format
         sv           call structural variations
         indelpe      post-processing indelpe output
\n);
}

#
# Other utilities
#

sub test_file
{
	my ($file) = @_;
	my $gzip = gwhich('gzip') || die;
	my ($is_gzip, $fh, $first);
	if (system("$gzip -l $file >/dev/null 2>&1")) {
		$is_gzip = 0;
		open($fh, $file);
	} else {
		$is_gzip = 1;
		open($fh, "$gzip -dc $file |");
	}
	read($fh, $first, 1);
	close($fh);
	return ($is_gzip, $first);
}

sub run_cmd
{
	my ($cmd) = @_;
	warn("-- CMD: $cmd\n");
	system("$cmd") && die("** fail to run command '$cmd'");
}

# the following codes are copied from treefam::generic

sub dirname
{
	my $prog = shift;
	my $cwd = getcwd;
	return $cwd if ($prog !~ /\//);
	$prog =~ s/\/[^\s\/]+$//g;
	return $prog;
}
sub which
{
	my $file = shift;
	my $path = (@_)? shift : $ENV{PATH};
	return if (!defined($path));
	foreach my $x (split(":", $path)) {
		$x =~ s/\/$//;
		return "$x/$file" if (-x "$x/$file" && -f "$x/$file");
	}
	return;
}
sub gwhich
{
	my $progname = shift;
	my $addtional_path = shift if (@_);
	my $dirname = &dirname($0);
	my $tmp;

	chomp($dirname);
	if (-x $progname && -f $progname) {
		return abs_path($progname);
	} elsif (defined($addtional_path) && ($tmp = &which($progname, $addtional_path))) {
		return $tmp; # lh3: Does it work? I will come back to this later
	} elsif (defined($dirname) && (-x "$dirname/$progname" && -f "$dirname/$progname")) {
		return abs_path("$dirname/$progname");
	} elsif (($tmp = &which($progname))) { # on the $PATH
		return $tmp;
	} else {
		warn("[gwhich] fail to find executable $progname anywhere.");
		return;
	}
}
