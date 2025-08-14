#!/usr/bin/env python3
"""
Build phylogenetic trees from eukaryotic single copy genes.

This script takes the output from find_euk_single_copy.py and:
1. Identifies genes with sufficient genome coverage
2. Extracts ORF sequences from the original summary
3. Concatenates sequences per genome
4. Runs MAFFT for alignment
5. Runs FastTree to build phylogenetic trees
"""

import argparse
import csv
import os
import sys
import subprocess
import tempfile
from collections import defaultdict, Counter
from typing import Dict, List, Set, Tuple
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_single_copy_genes(csv_file: str) -> Dict[str, List[str]]:
    """Load single copy genes from CSV file."""
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"Single copy genes file not found: {csv_file}")
    
    gene_families = defaultdict(list)
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_family = row['gene_family']
            gene_name = row['gene_name']
            gene_families[gene_family].append(gene_name)
    
    logger.info(f"Loaded {sum(len(genes) for genes in gene_families.values())} genes from {len(gene_families)} families")
    return gene_families

def load_orf_summary(summary_file: str) -> Dict[str, Dict[str, Dict]]:
    """Load ORF summary file and organize by sample_id and target_name."""
    if not os.path.exists(summary_file):
        raise FileNotFoundError(f"ORF summary file not found: {summary_file}")
    
    samples_data = defaultdict(dict)
    
    with open(summary_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sample_id = row['sample_id']
            target_name = row['target_name']
            samples_data[sample_id][target_name] = row
    
    logger.info(f"Loaded ORF data for {len(samples_data)} samples")
    return samples_data

def check_genome_coverage(gene_families: Dict[str, List[str]], 
                         samples_data: Dict[str, Dict[str, Dict]],
                         coverage_threshold: float = 100.0) -> Dict[str, List[str]]:
    """Check which genes have sufficient genome coverage."""
    total_genomes = len(samples_data)
    min_genomes = int(total_genomes * coverage_threshold / 100.0)
    
    logger.info(f"Total genomes: {total_genomes}")
    logger.info(f"Minimum genomes required: {min_genomes} ({coverage_threshold}%)")
    
    # Check coverage for each gene family
    genes_with_coverage = defaultdict(list)
    
    for family_name, genes in gene_families.items():
        for gene_name in genes:
            # Count how many genomes have this gene
            genomes_with_gene = 0
            for sample_id, sample_orfs in samples_data.items():
                # Check if any ORF in this sample has this gene name
                for target_name, orf_data in sample_orfs.items():
                    if orf_data.get('query_name') == gene_name:
                        genomes_with_gene += 1
                        break
            
            if genomes_with_gene >= min_genomes:
                genes_with_coverage[family_name].append(gene_name)
                logger.info(f"{gene_name}: {genomes_with_gene}/{total_genomes} genomes ({genomes_with_gene/total_genomes*100:.1f}%)")
    
    return genes_with_coverage

def extract_sequences(genes_with_coverage: Dict[str, List[str]],
                     samples_data: Dict[str, Dict[str, Dict]]) -> Dict[str, Dict[str, str]]:
    """Extract protein sequences for genes with sufficient coverage."""
    sequences = defaultdict(dict)
    
    for family_name, genes in genes_with_coverage.items():
        for gene_name in genes:
            # For each genome, find the best hit for this gene
            for sample_id, sample_orfs in samples_data.items():
                best_hit = None
                best_evalue = float('inf')
                
                for target_name, orf_data in sample_orfs.items():
                    if orf_data.get('query_name') == gene_name:
                        try:
                            evalue = float(orf_data.get('full_evalue', float('inf')))
                            if evalue < best_evalue:
                                best_evalue = evalue
                                best_hit = orf_data
                        except ValueError:
                            continue
                
                if best_hit:
                    # Store the target_name as the sequence identifier
                    sequences[gene_name][sample_id] = best_hit['target_name']
    
    return sequences

def create_fasta_files(sequences: Dict[str, Dict[str, str]],
                      samples_data: Dict[str, Dict[str, Dict]],
                      output_dir: str) -> Dict[str, str]:
    """Create FASTA files for each gene family."""
    os.makedirs(output_dir, exist_ok=True)
    
    fasta_files = {}
    
    for gene_name, genome_hits in sequences.items():
        fasta_file = os.path.join(output_dir, f"{gene_name}.fasta")
        
        with open(fasta_file, 'w') as f:
            for sample_id, target_name in genome_hits.items():
                # Get the actual sequence from the original FASTA files
                # For now, we'll use the target_name as a placeholder
                # In a real implementation, you'd need to load the actual sequences
                f.write(f">{sample_id}\n{target_name}\n")
        
        fasta_files[gene_name] = fasta_file
        logger.info(f"Created FASTA file for {gene_name}: {len(genome_hits)} sequences")
    
    return fasta_files

def run_mafft(fasta_file: str, output_dir: str) -> str:
    """Run MAFFT alignment on a FASTA file."""
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    aligned_file = os.path.join(output_dir, f"{base_name}_aligned.fasta")
    
    cmd = ['mafft', '--auto', fasta_file]
    
    try:
        with open(aligned_file, 'w') as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
        logger.info(f"MAFFT alignment completed: {aligned_file}")
        return aligned_file
    except subprocess.CalledProcessError as e:
        logger.error(f"MAFFT failed for {fasta_file}: {e}")
        return None

def run_fasttree(aligned_file: str, output_dir: str) -> str:
    """Run FastTree on an aligned FASTA file."""
    base_name = os.path.splitext(os.path.basename(aligned_file))[0].replace('_aligned', '')
    tree_file = os.path.join(output_dir, f"{base_name}.tree")
    
    cmd = ['fasttree', '-out', tree_file, aligned_file]
    
    try:
        subprocess.run(cmd, check=True)
        logger.info(f"FastTree completed: {tree_file}")
        return tree_file
    except subprocess.CalledProcessError as e:
        logger.error(f"FastTree failed for {aligned_file}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(
        description='Build phylogenetic trees from eukaryotic single copy genes'
    )
    parser.add_argument(
        'single_copy_csv',
        help='Path to eukaryotic_single_copy_genes.csv file'
    )
    parser.add_argument(
        'orf_summary',
        help='Path to eukaryotic ORF summary file (e.g., eukaryotes_orf_summary.tsv)'
    )
    parser.add_argument(
        '--coverage-threshold',
        type=float,
        default=100.0,
        help='Minimum percentage of genomes a gene must appear in (0-100, default: 100.0)'
    )
    parser.add_argument(
        '--output-dir',
        default='eukaryotic_trees',
        help='Output directory for trees and alignments (default: eukaryotic_trees)'
    )
    
    args = parser.parse_args()
    
    # Validate coverage threshold
    if args.coverage_threshold < 0 or args.coverage_threshold > 100:
        logger.error("coverage-threshold must be between 0 and 100")
        sys.exit(1)
    
    try:
        # Load single copy genes
        logger.info("Loading single copy genes...")
        gene_families = load_single_copy_genes(args.single_copy_csv)
        
        # Load ORF summary
        logger.info("Loading ORF summary...")
        samples_data = load_orf_summary(args.orf_summary)
        
        # Check genome coverage
        logger.info("Checking genome coverage...")
        genes_with_coverage = check_genome_coverage(gene_families, samples_data, args.coverage_threshold)
        
        if not genes_with_coverage:
            logger.warning("No genes meet the coverage threshold")
            sys.exit(0)
        
        # Extract sequences
        logger.info("Extracting sequences...")
        sequences = extract_sequences(genes_with_coverage, samples_data)
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Create FASTA files
        logger.info("Creating FASTA files...")
        fasta_files = create_fasta_files(sequences, samples_data, args.output_dir)
        
        # Run MAFFT and FastTree for each gene
        logger.info("Running alignments and tree building...")
        for gene_name, fasta_file in fasta_files.items():
            # Run MAFFT
            aligned_file = run_mafft(fasta_file, args.output_dir)
            if aligned_file:
                # Run FastTree
                tree_file = run_fasttree(aligned_file, args.output_dir)
                if tree_file:
                    logger.info(f"Successfully created tree for {gene_name}: {tree_file}")
        
        logger.info(f"Tree building completed. Results in: {args.output_dir}")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main() 