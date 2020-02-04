#!/usr/local/bin/perl

# Read a fastq file, extract the numeric tag from the read name, write
# the fastq reads to separate files by tag.

# NOTE: This script was cloned from split_tags.pl. split_tags.pl
# expects the tag to be the last bases of the first read. It processes
# fastq files in pairs, because the tag is found in only the first
# read fastq. This script, OTOH, expects a numeric tag to be part of
# the read name.  It can therefore process a single fastq file. All of
# those differences argued for creating a new script, rather than
# further encumbering split_tags.pl with more processing options.

# See &usage for a description of the output file naming convention.

use strict;
use warnings;
use English qw(-no_match_vars);
use Carp;
use Getopt::Long;

my $VERSION = '20090824.01';

sub usage;
sub initialise;
sub process;
sub read_fastq($);

my $opts = initialise;
my $verbose = $opts->{'verbose'} || 0;

process;

exit;

# ----------------------------------------------------------------------
sub usage {

  ## no critic

  my ($short_name) = $PROGRAM_NAME =~ /([^\/]+)$/xm;

  print STDERR "\n";
  print STDERR "$short_name version $VERSION\n";
  print STDERR "\n";
  print STDERR "  options:\n";
  print STDERR "\n";
  print STDERR "    --verbose      include step-by-step stats\n";
  print STDERR "    --stem         output filename prefix (see below)\n";
  print STDERR "    --tags         comma-separated list of tags to output (def all)\n";
  print STDERR "    --no-zero      do not output tag zero (shorter form of '--tags 1,2,3,4,5,6,7,8,...')\n";
  print STDERR "    --help         print this message and quit\n";
  print STDERR "\n";
  print STDERR "  Multiple input fastq files can be specified on the command line.\n";
  print STDERR "  The output file names will by default be the first input file name\n";
  print STDERR "  with the tag number inserted before the last field in the name\n";
  print STDERR "  -- e.g., 1234_4_1.fastq becomes 1234_4_1_9.fastq for tag 9.\n";
  print STDERR "  That behaviour can be overridden using the --stem command line\n";
  print STDERR "  parameter. If '--stem xyz' is specified, the output file name for\n";
  print STDERR "  tag 9 will be xyz.9.fastq\n";
  print STDERR "\n";
  print STDERR "  In any case, output files will be writen to the current directory.\n";
  print STDERR "\n";

}

# ----------------------------------------------------------------------
sub initialise {

  my %opts;
  my $rc = GetOptions(\%opts, 'help', 'verbose', 'stem:s', 'tags:s', 'no-zero');
  if ( ! $rc) {
    print {*STDERR} "\nerror in command line parameters\n" or croak 'print failed';
    usage;
    exit 1;
  }

  if (exists $opts{help}) {
    usage;
    exit;
  }

  return \%opts;

}

# ----------------------------------------------------------------------
sub process {

  my $stem;                            # see &usage
  if (exists $opts->{'stem'}) {
    $stem = $opts->{'stem'} . '.fastq';
  } elsif (@ARGV > 0) {
    $stem = $ARGV[0];
  } else {
    usage;
    exit 1;
  }

  my $fh = *ARGV;
  my @out_handles;
  my $nozero = exists $opts->{'no-zero'};

  my %want_tags = ();
  if (exists $opts->{'tags'}) {
    %want_tags = map {$_ => undef} split /,/xm, $opts->{'tags'};
  }

  while (my $lines = read_fastq($fh)) { # sets $lines

    if ($lines->[0] =~ /\#(\d+)/xm) {

      my $tag = $1;
      next if ($nozero and $tag == 0);
      next if (%want_tags > 0 and ! exists $want_tags{$tag});

      if ( ! $out_handles[$tag]) {      # is there an oputput file open for this tag?

        my ($front, $back) = $stem =~ /([^\/]+)\.([^\.]+)$/xm;
        my $outfile = "${front}_${tag}.${back}";
        open $out_handles[$tag], '>', $outfile or croak "could not open $outfile";
        $verbose && print STDERR "opened $outfile\n";

      }

      print {$out_handles[$tag]} join ("\n", @$lines), "\n";

    } else {
      croak 'no numeric index tag in read name: ' . $lines->[0];
    }

  }

  foreach my $handle (@out_handles) {
    if ($handle) {
      close $handle or croak "close failed";
    }
  }

}

# ----------------------------------------------------------------------
# Read one 4-line read from a fastq file whose filehandle is passed as
# a parameter. Check for @ and + in the expected places. Return a ref
# to a 4-element array, or undef at eof.

sub read_fastq ($) {

  my ($fh) = @_;
  my @lines;
  my $which = 0;
  
  while (my $line = <$fh>) { # sets $line

    chomp $line;
    push @lines, $line;

    if ($which == 0 and $line !~ /^\@/) {
      croak "not a read name at $.: $line\n";
    }
    if ($which == 2 and $line !~ /^\+/) {
      croak "not a + line at $.: $line\n";
    }

    last if ++$which >= 4;

  }

  return undef if $which == 0;
  croak "unexpected eof on fastq file" unless $which == 4;
  return \@lines;

}

