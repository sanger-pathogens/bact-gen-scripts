# TreeGubbins

A python script for identifying matches to resistance alleles in genome assemblies. Unlike many other methods map_resistome will highlight partial matches at contig boundaries that may be the result of assembly issues for repetitive resistance genes.

Requirements:

Requires python 2.x and the dendropy python library

Usage: TreeGubbins.py [options]

Options:
  -h, --help            show this help message and exit
  -t FILE, --tree=FILE  Input tree file (NOTE: tree must be fully bifurcating)
  -o OUTPUT, --output_prefix=OUTPUT
                        Output prefix
  -s FLOAT, --significance=FLOAT
                        Significance cutoff level [default= 0.05]
  -p INT, --permutations=INT
                        Number of permutations to run to test significance
                        [default= 100]
  -m, --midpoint        Midpoint root output tree pdf [default= False]
  -l LADDERISE, --ladderise=LADDERISE
                        ladderise tree (choose from right or left) [default=
                        none]
  -i, --iterative       Run in iterative mode, which allows nested clusters.
                        Default = do not run iteratively (much faster)

