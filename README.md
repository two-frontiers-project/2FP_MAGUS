# MAGUS: Pan-domain, holobiont characterization via co-assembly binning emphasizing low abundance organisms

## Background

The term "holobiont" refers to the assemblage of all organisms that make up a single meta-organism. This could be humans and their resident microbes, corals and their native viruses and algae, or any other number of metagenomic ecosystems. Given the complex interplay between microbes and their hosts, studying the macroscale holobiont instead of individual systems in isolation is critical for understanding ecosystem dynamics on a biologically meaningful scale. 

However, in DNA sequencing-based metagenomic analysis, we tend to bias our efforts to studying only high abundance organisms within a specific branch of the tree of life. For example, metagenomic binning in pursuit of resolving species genomes tends to emphasize bacteria, despite the fact that bacteria rarely exist in nature only with other bacteria -- this usually only happens in lab settings that are designed by us humans. Further, our metagenomic sequencing, when considered on a sample-by-sample basis, is highly biased to high-abundance microbes, ones that dominate the captured sequences from a sample. As a result, binning and other methods miss some ecologically critical, low-abundance organisms.

Here, we provide MAGUS -- a pipeline designed for pan-domain analysis of the holobiont, capturing low abundance organisms via a co-assembly-based method. Taking sequencing datasets as input, we return 1) bacterial/archaeal, 2) viral, and 3) putative eukaryotic MAGs. Our initial focus is on systematic characterization of coral holobionts, but in principle this toolkit can be used for any ecosystem where low abundance organisms are of interest.

## Approach

![Alt text](images/magus_workflow.png)

MAGUS takes a multi-pass, single then co-assembly approach to identify putative metagenomic bins. For assembly, we use a modified megahit implementation that [GABE DESCRIBE]. Following single sample assembly, we run MetaBAT2 in order to identify putative bins. Samples are then selected for coassembly based on jaccard distance between assembled contigs, with the hypothesis being that samples with a certain degree of similarity will be more likely to assemble low abundance bins that were missed in single assembly. Following coassembly, binning is attempted both with aggregated coverage (i.e.,, without alignment to compute individual coverages) as well as distributed coverage (using alignment). We use CheckM2 to identify putative bacteria/archaea, CheckV to get viruses, and two Eukaryotic binners (EukCC and EukRep) to identify putative eukaryotic bins. Viral genomes are dereplicated at the 90% identity level. Identical bacterial/archael genomes are consolidated between coassembled and single assembled samples. Eukaryotic genomes are not dereplicated. 

## Installation

MAGUS was built and runs quite happily on Fedora Linux 40 (Workstation Edition). 

When we preprint the paper, we'll have a docker container you can just pull.

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

Now, if you want to run EukRep, you will have to run it into its own conda environment. The models it uses for genome prediction were trained with an older version of ```scikit-learn``` than is required by CheckM2, which is a more critical tool for MAGUS overall. We have forked EukRep's repo and are playing around with updating the dependency versioning, but for the time being it's easier to install EukRep in a conda environment (eukrep-env is what MAGUS looks for by default, you can update this in the config files). In the find-euks step, MAGUS will activate then deactivate this environment as needed.

Here's the link to EukRep:

```https://github.com/patrickwest/EukRep```

To install:

```
conda create -y -n eukrep-env -c bioconda scikit-learn==0.19.2 eukrep
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

## A note on runtimes and memory usage

This is a **wildly** memory intensive piece of software. It is not meant to be run on personal machines. Single assemblies are easy enough, but co-assemblies can easily require 3+ terabytes of RAM. It can take weeks to work through only a few hundred samples, even if they're sub 100M PE reads. We're working on some methods for clever downsampling to speed things up in the coassembly step without losing critical information, but in the meantime we recommend leveraging HPC systems, cloud credits, and leveraging spot instances. If you have specific challenges, please feel free to reach out to ```info at two frontiers dot org.```

## Commands and arguments

MAGUS exposes the following sub-commands through `magus <command>`:

| Command | Description |
| --- | --- |
| `magus qc` | Run read QC and compression with `shi7_trimmer` and `minigzip`. |
| `magus assemble-hosts` | Generate host assemblies and clustering metadata used to mask dominant genomes prior to metagenomic processing. |
| `magus subsample-reads` | Down-sample reads with `seqtk` and emit an updated config file. |
| `magus filter-reads` | Remove reads matching host references using xtree `.perq` output. |
| `magus taxonomy` | Run xtree on raw reads to generate coverage, taxonomy, and abundance summaries. |
| `magus single-assembly` | Assemble each sample independently with MEGAHIT. |
| `magus binning` | Bin single-assembly contigs with MetaBAT2 and assess quality with CheckM2. |
| `magus cluster-contigs` | Collect filtered contigs, perform clustering, and build the co-assembly task list. |
| `magus coassembly` | Execute co-assemblies for grouped samples. |
| `magus coassembly-binning` | Bin co-assembly results and evaluate bins with CheckM2. |
| `magus finalize-bacterial-mags` | Merge redundant bacterial/archaeal MAGs across single and co-assemblies. |
| `magus find-viruses` | Merge contigs, run CheckV, and dereplicate viral genomes. |
| `magus find-euks` | Identify putative eukaryotic genomes using EukRep and EukCC. |
| `magus dereplicate` | Dereplicate MAGs with lingenome and canolax5. |
| `magus call-orfs` | Call and annotate open reading frames across MAG collections. |
| `magus build-gene-catalog` | Build pan-genome gene catalogs from ORF predictions. |
| `magus filter-mags` | Filter MAGs using xtree alignment statistics. |
| `magus build-tree` | Construct phylogenetic trees from concatenated single-copy genes. |

### Command arguments

Key options for each command are summarised below (see `magus <command> --help` for complete details):

- **qc** (`magus qc`)
  - `--config`: TSV describing raw reads.
  - `--slurm_config`: Optional Slurm settings for batch execution.
  - `--max_workers`: Parallel samples to process at once.
  - `--mode`: `local` or `slurm` execution mode.
- **assemble-hosts** (`magus assemble-hosts`)
  - `--config`: TSV of reads to subsample for host assembly.
  - `--threads` / `--max_workers`: Compute resources for preprocessing.
  - `--ksize`: Override the number of host clusters to recover.
  - `--output_config`: Destination for the host-filtering config.
  - `--tmpdir`: Working directory for intermediate host files.
- **subsample-reads** (`magus subsample-reads`)
  - `--config`: Input sequencing config.
  - `--outdir`: Directory for subsampled reads.
  - `--out_config`: Path for the updated config file.
  - `--depth`: Target read count per sample.
  - `--threads` / `--max_workers`: Parallelism controls.
- **filter-reads** (`magus filter-reads`)
  - `--config`: Sequencing config aligning filenames to read pairs.
  - `--perq_dir`: Directory containing xtree `.perq` files.
  - `--output_dir`: Where filtered reads are written.
  - `--min_kmers`: Minimum k-mer evidence retained per read.
  - `--threads` / `--max_workers`: Filtering parallelism.
- **taxonomy** (`magus taxonomy`)
  - `--config`: TSV mapping sample IDs to read locations.
  - `--db`: xtree database to query.
  - `--output`: Output directory for alignments and summaries.
  - `--threads` / `--max_workers`: Worker and thread counts.
  - `--coverage-cutoff` and `--skip-*` flags: Control downstream filtering and optional output files.
- **single-assembly** (`magus single-assembly`)
  - `--config`: QC'd read configuration.
  - `--threads`: Threads passed to MEGAHIT.
  - `--max_workers`: Number of samples assembled in parallel.
  - `--mode` / `--slurm_config`: Toggle and configure Slurm execution.
- **binning** (`magus binning`)
  - `--config`: Assembly config file (expects `filename`, `pe1`, `pe2`).
  - `--asmdir`: Location of assemblies to bin.
  - `--tmpdir`: Scratch directory for binning intermediates.
  - `--threads` / `--max_workers`: Resources for MetaBAT2 and CheckM2.
  - `--checkm_db`: Optional CheckM database override.
  - Quality thresholds (`--completeness`, `--contamination`, `--low-quality`, `--medium-quality`, `--high-quality`) and `--restart` flags to resume stages.
- **cluster-contigs** (`magus cluster-contigs`)
  - `--config`: Same config used for single assemblies.
  - `--asmdir` / `--magdir`: Source assemblies and MAGs.
  - `--contig_dir` / `--combined_output`: Output paths for filtered contigs.
  - `--threads`: CPU threads for clustering utilities.
  - `--tmpdir`: Workspace for clustering output.
- **coassembly** (`magus coassembly`)
  - `--config`: Sequencing config.
  - `--coasm_todo`: Task list generated by `cluster-contigs`.
  - `--outdir`: Destination for co-assembly results.
  - `--tmpdir`: Working directory for co-assembly intermediates.
  - `--threads`: Threads allocated per co-assembly.
  - `--test_mode`: Relax filters and generate smaller test runs.
- **coassembly-binning** (`magus coassembly-binning`)
  - `--config`: Sequencing config.
  - `--coasm_outdir`: Root of co-assembly outputs to process.
  - `--tmpdir`: Scratch directory for binning.
  - `--threads` / `--max_workers`: Resource controls.
  - `--checkm_db`: Optional CheckM database override.
  - `--test_mode` and `--restart`: Control relaxed thresholds and resume behaviour.
- **finalize-bacterial-mags** (`magus finalize-bacterial-mags`)
  - `--singleassembly_mag_dir` / `--coasm_mag_dir`: Input MAG locations.
  - `--outdir`: Final MAG export directory.
  - `--threads`: CPU threads for dereplication.
  - `--tmpdir`: Temporary directory for merging work.
- **find-viruses** (`magus find-viruses`)
  - `--asm_paths` *or* `--config`: Supply assemblies to scan.
  - `--checkv_db`: CheckV database location.
  - `--combined_contig_file` / `--filtered_contig_file`: Filenames for merged contigs.
  - `--quality`: CheckV quality tiers to retain.
  - `--tmpdir`: Scratch directory for CheckV runs and dereplication.
  - `--restart cleanup`: Resume downstream processing after CheckV completes.
- **find-euks** (`magus find-euks`)
  - `--bin_dirs`: Pipe-delimited directories to search for bins.
  - `--wildcards`: Patterns (also pipe-delimited) to match candidate bins.
  - `--size_threshold`: Minimum bin size included (bp).
  - `--euk_binning_outputdir`: Output directory for eukaryotic bins.
  - `--dblocs`: Mapping file for database paths (expects `eukccdb`).
  - `--max_workers` / `--threads`: Concurrency settings.
  - `--skip_eukrep`, `--skip_eukcc`, `--eukrep_env`, `--checkm2_file`: Control filtering stages.
- **dereplicate** (`magus dereplicate`)
  - `--mag_dir`: Glob pointing to MAGs to dereplicate.
  - `--tmpdir`: Working directory for lingenome and canolax intermediates.
  - `--threads`: Threads for canolax5.
  - `--extension`, `--wildcard`, `--output`, `--kmer_size`, `--max_genome_size`: Options controlling dereplication scope and behaviour.
- **call-orfs** (`magus call-orfs`)
  - Supports configs or globs via `--config` or `--mag_dir`/`--wildcard`.
  - Control domains, annotation targets, and runtime via `--domain`, `--output_directory`, `--max_workers`, `--threads`, `--extension`, and `--force`.
  - Annotation flags include `--hmmfile`, `--annotation-fullseq-evalue`, `--annotation-domain-evalue`, `--suffix`, `--eukdb`, `--cleanup`, and `--restart`.
- **build-gene-catalog** (`magus build-gene-catalog`)
  - `--summary-file`: ORF summary table.
  - `--faa-dir`: Directory of amino-acid FASTA files.
  - `--output-dir`: Catalog output root.
  - `--threads`: Threads for MMseqs2.
  - Filtering knobs: `--evalue-cutoff`, `--identity-threshold`, `--coverage-threshold`, `--identity-only`, `--multi-sample-single-copy`.
  - `--tmpdir`: Working directory handed to MMseqs2.
- **filter-mags** (`magus filter-mags`)
  - `--output-dir`: Filtered MAG destination.
  - `--perq-dir`: Directory housing xtree `.perq` files.
  - `--mag-dir`: MAG directory to filter.
  - `--kmer-threshold`: Evidence cutoff for retaining contigs.
- **build-tree** (`magus build-tree`)
  - Positional `fasta_dir`: Directory of per-genome alignments.
  - `--gene-list`: Single-copy gene definitions.
  - `--orf-data`: ORF summary file.
  - `--coverage-threshold`, `--evalue-cutoff`, `--trimal-cutoff`: Filtering parameters for alignments.
  - `--iqtree`: Toggle IQ-TREE instead of FastTree.
  - `--output-dir`: Output directory for trees.
  - `--genome-list`: Optional list of genomes to include.
  - `--threads`: Parallelism for tree building.

## Conda and python dependencies 

## Other software requirements

The external software that we use (e.g., tools not found in conda, like our version of MegaHIT) is all found in the bin/ directory. This should be added to the path on installation. A brief description of each tool is here:

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
| **akmer102**  | Calculates k-mer frequencies and distances for genome comparisons.                                           |
| **bestmag2** | Selects the "best" MAGs based on quality metrics and coverage information.                          |
| **spamw2**     | Clusters genomes based on pairwise jaccard distances.                                                        |

## License

MAGUS is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) license. This means it is free for academic and non-commercial use. For commercial use, please contact us to discuss licensing terms.

## Authors

Gabe Al-Ghalith ```(gabe at two frontiers dot org)```
Braden Tierney ```(braden at two frontiers dot org)```

## Contact

If you have questions, reach out to the authors and/or ```info at two frontiers dot org```, which will reach more of us at 2FP.

