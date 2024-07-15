# 2FP magus workflow for systematical genomic analysis of the coral holobiont

## Strategy

The goal here is to provide already quality-controlled, deep sequencing data from a coral and get out 1) tables of taxonomic abundances from short read alignment 2) Putative coral/bacterial/viral/algal genomes and summary statistics (e.g., likely taxonomy) 3) maybe some functional stuff too

That works out to fourish workflows, all ideally written in python and callable at the command line:

magus taxonomy

magus assemble

magus bin

magus genelevelanalysis

We might also want to have a "compare" function to look at known vs. unknown holobiont genomes. 

Each of those should have detailed documentation and be able to take multiple input types (e.g., short reads, long reads, both, etc).

The first step, though, is to build some python infrastructure to accomodate each of these modules. I would recommend we model this after gtdb-tk (https://github.com/Ecogenomics/GTDBTk/tree/master/gtdbtk). Pretty straightforward, each subcommand is a different python script. I can organize this repository accordingly pretty easily. 
