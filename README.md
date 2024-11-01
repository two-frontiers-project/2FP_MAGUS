# MAGUS: Pan-domain, holobiont characterization via co-assembly binning emphasizing low abundance organisms

note from BT -- this probably won't install correctly yet, I haven't updated yml 

## Background

The term "holobiont" refers to the assemblage of all organisms that make up a single meta-organism. This could be humans and their resident microbes, corals and their native viruses and algae, or any other number of metagenomic ecosystems. Given the complex interplay between microbes and their hosts, studying the macroscale holobiont instead of individual systems in isolation is critical for understanding ecosystem dynamics on a biologically meaningful scale. 

However, in DNA sequencing-based metagenomic analysis, we tend to bias our efforts to studying only high abundance organisms within a specific branch of the tree of life. For example, metagenomic binning in pursuit of resolving species genomes tends to emphasize bacteria, despite the fact that bacteria rarely exist in nature only with other bacteria -- this usually only happens in lab settings that are designed by us humans. Further, our metagenomic sequencing, when considered on a sample-by-sample basis, is highly biased to high-abundance microbes, ones that dominate the captured sequences from a sample. As a result, binning and other methods miss some ecologically critical, low-abundance organisms.

Here, we provide MAGUS -- a pipeline designed for pan-domain analysis of the holobiont, capturing low abundance organisms via a co-assembly-based method. Taking sequencing datasets as input, we return 1) bacterial/archaeal, 2) viral, and 3) putative eukaryotic MAGs. Our initial focus is on systematic characterization of coral holobionts, but in principle this toolkit can be used for any ecosystem where low abundance organisms are of interest.

## Approach

![Alt text](images/magus_workflow.png)

MAGUS takes a multi-pass, single then co-assembly approach to identify putative metagenomic bins. For assembly, we use a modified megahit implementation that [GABE DESCRIBE]. Following single sample assembly, we run MetaBAT2 in order to identify putative bins. Samples are then selected for coassembly based on jaccard distance between assembled contigs, with the hypothesis being that samples with a certain degree of similarity will be more likely to assemble low abundance bins that were missed in single assembly. Following coassembly, binning is attempted both with aggregated coverage (i.e.,, without alignment to compute individual coverages) as well as distributed coverage (using alignment). We use CheckM2 to identify putative bacteria/archaea, CheckV to get viruses, and two Eukaryotic binners (EukCC and EukRep) to identify putative eukaryotic bins. Viral genomes are dereplicated at the 90% identity level. Identical bacterial/archael genomes are consolidated between coassembled and single assembled samples. Eukaryotic genomes are not dereplicated. 

## Installation

MAGUS was built and runs quite happily on [Fedora Linux 40 (Workstation Edition)]. 

```bash
git clone https://github.com/two-frontiers-project/2FP_MAGUS.git   
cd 2FP_MAGUS
conda create -n magus -f magus.yml
conda activate magus
pip install .
```

You'll also need to install some databases. Use this function:

```
magus install_db --path /path/to/dbfiles
```

## Database and config setup

You'll notice MAGUS is designed to not be run in a single click (we have no end-to-end implmentation) -- this is intentional, as not all users will need to run it fully, and the co-assembly steps are extraordinarily memory intensive. 

Additionally, we parameterize the different functions based on config files (located by default in the config directory). These provide paths to the sequencing files you want to process, as well as the raw database locations (this is to avoid muddying up your paths and prevent having to manually specify database locations in each steps). 

So, before running MAGUS **be sure that you update the raw_config and db_locs config files with the appopropriate paths to raw data and databases on your system.**

## Input and output

In its maximal form, when you run MAGUS you'll end up with a single directory that contains:

1. Bacterial MAGS
2. Viral Genomes
3. Putative Eukaryotic MAGS
4. Taxonomic information on every MAG
5. Summarized quality statistics for each MAG

You'll also end up, of course, with all your assemblies, and in the near future you'll have gene catalogs, functional annotations, and phylogenies.

## A note on runtimes

This is a memory intensive piece of software. Single assemblies are easy enough, but co-assemblies can easily require 3+ terabytes of RAM. It can take weeks to work through only a few hundred samples, even if they're sub 100M PE reads. We're working on some methods for clever downsampling to speed things up in the coassembly step without losing critical information, but in the meantime we recommend leveraging HPC systems, cloud credits, and leveraging spot instances. If you have specific challenges, please feel free to reach out to ```info at two frontiers dot org.```

## Commands and arguments

The sequence of commands to run the full pipeline is as follows:

| Command       | Description                                                                                                  |
|----------------|--------------------------------------------------------------------------------------------------------------|
| magus qc |   Read quality control and compression    |
| magus single-assembly |    Assemble samples one at a time        |
| magus single-binning |  Bin genomes with MetaBAT2  and run CheckM2  |
| magus cluster-contigs |   Identify samples for potential co-assembly  |
| magus coassembly |  Co-assemble samples|
| magus coassembly-binning |    Run MetaBAT2 and CheckM2 on coassembled bins    |
| magus finalize-bacterial-mags |  Filter redundant bacterial/archaeal MAGS identified in both single and coassembly binning.  |
| magus find-viruses | Identify viral contigs with CheckV |
| magus find-euks | Identify putative eukaryotic bins with EukRep and EukCC   |

In the near future we'll release our gene catalog modules that enable functional comparison of various genomes and construction of phylogenetic trees.

## Command arguments

| Command                  | Argument                   | Description                                                    |
|--------------------------|----------------------------|----------------------------------------------------------------|
| **qc**      | `--config`                | Location of the configuration file containing the raw reads                            |
|                          | `--max_workers`           | Number of parallel jobs to run simultaneously                  |
|                          | `--threads`               | Number of threads assigned per job                             |
| **single-assembly**      | `--config`                | Location of the configuration file containing the qc'd reads                            |
|                          | `--max_workers`           | Number of parallel jobs to run simultaneously                  |
|                          | `--threads`               | Number of threads assigned per job                             |
| **single-binning**       | `--config`                | Location of the configuration file containing the qc'd reads                             |
|                          | `--threads`               | Number of threads assigned per job                             |
|                          | `--asmdir`                | Directory for the assembly output                              |
|                          | `--max_workers`           | Number of parallel jobs to run simultaneously                  |
|                          | `--tmp_dir`               | Temporary directory for intermediate files                     |
|                          | `--test_mode`             | Enables test mode for debugging                                |
| **cluster-contigs**      | `--config`                | Location of the configuration file containing the qc'd reads                              |
|                          | `--threads`               | Number of threads assigned per job                             |
|                          | `--contig_dir`            | Directory containing contigs for clustering                    |
|                          | `--combined_output`       | File to store combined contig output                           |
|                          | `--tmp_dir`               | Temporary directory for intermediate files                     |
| **coassembly**           | `--config`                | Location of the configuration file containing the qc'd reads                              |
|                          | `--coasm_todo`            | To-do list or specification file for co-assembly tasks         |
|                          | `--outdir`                | Output directory for co-assembly results                       |
|                          | `--tmp_dir`               | Temporary directory for intermediate files                     |
|                          | `--threads`               | Number of threads assigned per job                             |
|                          | `--test_mode`             | Enables test mode for debugging                                |
| **coassembly-binning**   | `--config`                | Location of the configuration file containing the qc'd reads                              |
|                          | `--coasm_outdir`          | Directory to store co-assembly output                          |
|                          | `--tmp_dir`               | Temporary directory for intermediate files                     |
|                          | `--threads`               | Number of threads assigned per job                             |
|                          | `--checkm_db`             | Path to CheckM database for quality assessment                 |
|                          | `--max_workers`           | Number of parallel jobs to run simultaneously                  |
|                          | `--test_mode`             | Enables test mode for debugging                                |
| **find-viruses**         | `--asm_dir`               | Directory for single-assembly files                            |
|                          | `--coasm_dir`             | Directory for co-assembly files                                |
|                          | `--combined_contig_file`  | File containing combined contigs for analysis                  |
|                          | `--filtered_contig_file`  | File containing filtered contigs based on length and quality   |
|                          | `--min_length`            | Minimum length cutoff for viral contigs                        |
|                          | `--max_length`            | Maximum length cutoff for viral contigs                        |
|                          | `--threads`               | Number of threads assigned per job                             |
|                          | `--quality`               | Quality threshold for viral identification                     |
|                          | `--tmp_dir`               | Temporary directory for intermediate files                     |
|                          | `--dblocs`                | Database locations for viral identification                    |
| **find-euks**            | `--coasm_dir`             | Directory for co-assembly files                                |
|                          | `--asm_dir`               | Directory for single-assembly files                            |
|                          | `--size_threshold`        | Size threshold for eukaryotic genome detection                 |
|                          | `--euk_binning_outputdir` | Output directory for eukaryotic binning results                |
|                          | `--dblocs`                | Database locations for eukaryotic genome identification        |
|                          | `--max_workers`           | Number of parallel jobs to run simultaneously                  |
|                          | `--threads`               | Number of threads assigned per job                             |
|                          | `--skip_eukrep`           | Option to skip EukRep-based filtering                          |
|                          | `--eukrep_env`            | Conda environment for running EukRep                           |
|                          | `--skip_eukcc`            | Option to skip EukCC-based filtering                           |
| **finalize-bacterial-mags** | `--singleassembly_mag_dir` | Directory for single-assembly MAGs                          |
|                          | `--coasm_mag_dir`         | Directory for co-assembly MAGs                                 |
|                          | `--outdir`                | Output directory for finalized bacterial MAGs                  |
|                          | `--threads`               | Number of threads assigned per job                             |
|                          | `--tmp_dir`               | Temporary directory for intermediate files                     |


## Additional software requirements

The external software that we use (e.g., our version of MegaHIT, etc) is all found in the bin/ directory. This should be added to the path on installation. A brief description of each tool is here:

| Software       | Description                                                                                                  |
|----------------|--------------------------------------------------------------------------------------------------------------|
| **shi7_trimmer** | Trims adapter sequences from raw sequencing reads.                                                         |
| **minigzip**   | Compresses files using a faster gzip algorithm.                              |
| **checkm2**    | Assesses the quality and completeness of metagenome-assembled genomes (MAGs).                                |
| **megahit-g**    | Custom megahit implementation that XXX. Performs metagenomic assembly, constructing longer sequences (contigs) from short sequencing reads.          |
| **sorenson-g** | Estimates sequencing coverage of contigs using read alignments.                                              |
| **metabat2**   | Bins assembled contigs into putative MAGs based on coverage and sequence composition.                        |
| **fac**        | Filters contigs based on length and coverage.                                                                |
| **lingenome**  | Generates concatenated genomes from individual FASTA files.                                                  |
| **akmer100b**  | Calculates k-mer frequencies and distances for genome comparisons.                                           |
| **bestmag2** | Selects the "best" MAGs based on quality metrics and coverage information.                          |
| **spamw2**     | Clusters genomes based on pairwise jaccard distances.                                                        |





