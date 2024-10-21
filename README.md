








# 2FP magus workflow for systematic genomic analysis of the coral holobiont (and other similar ecosystems)

## Strategy

The goal here is to provide already quality-controlled, deep sequencing data from a coral and get out 1) tables of taxonomic abundances from short read alignment 2) Putative coral/bacterial/viral/algal genomes and summary statistics (e.g., likely taxonomy) 3) maybe some functional stuff too.

The core challenge is that sequencing data from corals contains genomes from across the tree of life: corals, bacteria, viruses, algal symbionts, and other eukaryotes. We need to first quantify the abundance of each high-level (eg domain) and low-level (eg species) taxonomy by short read alignment (likely with kraken2), then assemble the reads, then bin the genomes. We will also functionally annotate the genomes.

STRUCTURE

I'd like to model the codebase structure on GTDBTK (https://github.com/Ecogenomics/GTDBTk/tree/master/gtdbtk/external) -- so we'll need to have classes for each of the external tools, like they do in the link. A good place to start is to write these scripts.

I think we can generate four workflows, all ideally written in python and callable at the command line:

magus taxonomy
  - This gets taxonomy on short read data and computes the domain and species level composition 

  - I think I want to run kraken2 and mimic the approach used by phanta (https://github.com/bhattlab/phanta), which implements some nice coverage cutoffs etc. We'll have a "general" database that uses all genomes out there, and we'll also build a coral-specific one. The option to align to different databases (path specific by user to start) will be critical.
  
magus assemble
  - This assembles everything with megahit or metaspades. The choice is up to the user at the CLI level. Specific parameters will be in a config file.

magus resolve genomes
 - This runs binning to separate out the components of the holobiont into different genomes (so it runs metabat2, checkv, a eukaryotic binner of some kind?)

magus genelevelanalysis
 - Functional/gene analysis of the bins

We might also want to have a "compare" function to look at known vs. unknown holobiont genomes. 

Each of those should have detailed documentation and be able to take multiple input types (e.g., short reads, long reads, both, etc).

The first step, though, is to build some python infrastructure to accomodate each of these modules. I would recommend we model this after gtdb-tk (https://github.com/Ecogenomics/GTDBTk/tree/master/gtdbtk). Pretty straightforward, each subcommand is a different python script. I can organize this repository accordingly pretty easily. 
