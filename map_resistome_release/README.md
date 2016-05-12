# map_resistome

PLEASE NOTE: THIS SCRIPT IS DEPRECATED. WE RECOMMEND THE USE OF ARIBA INSTEAD.

A python script for identifying matches to resistance alleles in genome assemblies. Unlike many other methods map_resistome will highlight partial matches at contig boundaries that may be the result of assembly issues for repetitive resistance genes.

Requirements:

Requires local intallations of SMALT and samtools. If these are not in your default PATH, please edit the SMALT_DIR and SAMTOOLS_DIR locations in the two python files.

Usage: 
map_resistome.py [options]

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -c FILE, --contigs=FILE
                        multifasta containing contigs to search in
  -f FILE, --forward=FILE
                        forward fastq file (may be zipped, but must end .gz)
  -r FILE, --reverse=FILE
                        reverse fastq file (may be zipped, but must end .gz)
  -g FILE, --genes=FILE
                        multifasta containing genes to search for
  -o FILE, --output=FILE
                        output prefix
  -i float, --id=float  minimum id to report match (excluding clipping due to
                        contig breaks). Must be between 0 and 1. [default =
                        0.9]
