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

MAGUS is designed to be run on [XXX].

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

## Deployment

You'll notice MAGUS is designed to not be run in a single click (we have no end-to-end implmentation) -- this is intentional, as not all users will need to run it fully, and the co-assembly steps are extraordinarily memory intensive. 

Additionally, we parameterize the different functions based on config files (located by default in the config directory). These provide paths to the sequencing files you want to process, as well as the raw database locations (this is to avoid muddying up your paths and prevent having to manually specify database locations in each steps). 

So, before running MAGUS **be sure that you update the raw_config and db_locs config files with the appopropriate paths to raw data and databases on your system.**








