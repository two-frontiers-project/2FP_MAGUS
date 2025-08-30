#!/usr/bin/env python3
"""
Build phylogenetic trees from eukaryotic single copy genes.

This script takes the output from find_euk_single_copy.py and:
1. Identifies genes with sufficient overall genome coverage
2. Extracts ORF sequences from the original summary
3. Concatenates ALL genes per genome into single sequences
4. Runs MAFFT for alignment on concatenated sequences
5. Trims the alignment with trimAl to remove gappy columns
6. Runs FastTree to build ONE phylogenetic tree
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
from Bio import SeqIO

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

def load_genome_list(genome_list_file: str = None) -> Set[str]:
    """Load list of genomes to include in analysis."""
    if not genome_list_file:
        return set()
    
    if not os.path.exists(genome_list_file):
        raise FileNotFoundError(f"Genome list file not found: {genome_list_file}")
    
    genomes = set()
    with open(genome_list_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                genomes.add(line)
    
    return genomes

def load_orf_summary(summary_file: str, genome_list: Set[str] = None) -> Dict[str, Dict[str, Dict]]:
    """Load ORF summary file and organize by sample_id and target_name, optionally filtering by genome list."""
    if not os.path.exists(summary_file):
        raise FileNotFoundError(f"ORF summary file not found: {summary_file}")
    
    samples_data = defaultdict(dict)
    
    with open(summary_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sample_id = row['sample_id']
            
            # Skip if genome list is provided and this sample is not in it
            if genome_list and sample_id not in genome_list:
                continue
                
            target_name = row['target_name']
            samples_data[sample_id][target_name] = row
    
    logger.info(f"Loaded ORF data for {len(samples_data)} samples")
    return samples_data

def load_fasta_sequences(fasta_dir: str, sample_id: str) -> Dict[str, str]:
    """Load protein sequences from a sample's .fas file using Biopython."""
    fas_file = os.path.join(fasta_dir, f"{sample_id}.fas")
    
    if not os.path.exists(fas_file):
        logger.warning(f"FASTA file not found: {fas_file}")
        return {}
    
    sequences = {}
    try:
        for record in SeqIO.parse(fas_file, "fasta"):
            sequences[record.id] = str(record.seq)
        logger.info(f"Loaded {len(sequences)} sequences from {fas_file}")
    except Exception as e:
        logger.error(f"Error loading FASTA file {fas_file}: {e}")
        return {}
    
    return sequences

def check_overall_genome_coverage(gene_families: Dict[str, List[str]], 
                                 samples_data: Dict[str, Dict[str, Dict]],
                                 coverage_threshold: float = 100.0) -> Tuple[List[str], Dict[str, Set[str]]]:
    """Check overall genome coverage across ALL genes combined."""
    total_genomes = len(samples_data)
    min_genomes = int(total_genomes * coverage_threshold / 100.0)
    
    logger.info(f"Total genomes: {total_genomes}")
    logger.info(f"Minimum genomes required: {min_genomes} ({coverage_threshold}%)")
    
    # Get all unique genes
    all_genes = []
    for genes in gene_families.values():
        all_genes.extend(genes)
    
    logger.info(f"Total genes to analyze: {len(all_genes)}")
    
    # For each genome, check which genes it has
    genome_gene_coverage = {}
    for sample_id, sample_orfs in samples_data.items():
        genes_in_genome = set()
        for target_name, orf_data in sample_orfs.items():
            query_name = orf_data.get('query_name', '')
            if query_name in all_genes:
                genes_in_genome.add(query_name)
        genome_gene_coverage[sample_id] = genes_in_genome
    
    # Count how many genomes have at least one gene
    genomes_with_any_gene = sum(1 for genes in genome_gene_coverage.values() if len(genes) > 0)
    
    logger.info(f"Genomes with at least one gene: {genomes_with_any_gene}/{total_genomes}")
    
    if genomes_with_any_gene < min_genomes:
        logger.warning(f"Insufficient genome coverage: {genomes_with_any_gene}/{total_genomes} < {min_genomes}")
        return [], {}
    
    return all_genes, genome_gene_coverage

def extract_and_concatenate_sequences(all_genes: List[str],
                                    genome_gene_coverage: Dict[str, Set[str]],
                                    samples_data: Dict[str, Dict[str, Dict]],
                                    fasta_dir: str,
                                    evalue_cutoff: float = 0.001) -> Dict[str, Dict]:
    """Extract and concatenate all gene sequences for each genome."""
    concatenated_sequences = {}
    
    # Sort genes to ensure consistent concatenation order
    sorted_genes = sorted(all_genes)
    
    for sample_id, genes_in_genome in genome_gene_coverage.items():
        if not genes_in_genome:  # Skip genomes with no genes
            continue
        
        # Load this sample's FASTA sequences
        sample_sequences = load_fasta_sequences(fasta_dir, sample_id)
        if not sample_sequences:
            continue
            
        # For each gene, find the best hit (lowest E-value) and extract sequence
        concatenated_seq = []
        gene_order = []
        
        for gene_name in sorted_genes:
            if gene_name in genes_in_genome:
                # Find best hit for this gene in this genome
                best_hit = None
                best_evalue = float('inf')
                best_target_name = None
                
                for target_name, orf_data in samples_data[sample_id].items():
                    if orf_data.get('query_name') == gene_name:
                        try:
                            evalue = float(orf_data.get('full_evalue', float('inf')))
                            if evalue < evalue_cutoff and evalue < best_evalue:
                                best_evalue = evalue
                                best_hit = orf_data
                                best_target_name = target_name
                        except ValueError:
                            continue
                
                if best_hit and best_target_name in sample_sequences:
                    # Extract the actual protein sequence
                    protein_seq = sample_sequences[best_target_name]
                    concatenated_seq.append(protein_seq)
                    gene_order.append(gene_name)
                    logger.debug(f"Added {gene_name} sequence (length: {len(protein_seq)}) for {sample_id}")
                else:
                    # No valid hit found
                    concatenated_seq.append('X' * 100)  # Placeholder for missing gene
                    gene_order.append(f"{gene_name}_MISSING")
            else:
                concatenated_seq.append('X' * 100)  # Placeholder for missing gene
                gene_order.append(f"{gene_name}_MISSING")
        
        # Join all sequences with a separator
        concatenated_sequences[sample_id] = {
            'sequence': '|'.join(concatenated_seq),
            'gene_order': gene_order,
            'genes_present': len([g for g in genes_in_genome if g in sorted_genes])
        }
    
    return concatenated_sequences

def create_concatenated_fasta(concatenated_sequences: Dict[str, Dict],
                             output_dir: str) -> str:
    """Create a single FASTA file with concatenated sequences."""
    os.makedirs(output_dir, exist_ok=True)
    
    fasta_file = os.path.join(output_dir, "concatenated_genes.fasta")
    
    with open(fasta_file, 'w') as f:
        for sample_id, data in concatenated_sequences.items():
            f.write(f">{sample_id}\n{data['sequence']}\n")
    
    # Also save gene order info
    info_file = os.path.join(output_dir, "gene_order.txt")
    with open(info_file, 'w') as f:
        f.write("Gene order in concatenated sequences:\n")
        for i, gene in enumerate(concatenated_sequences[list(concatenated_sequences.keys())[0]]['gene_order']):
            f.write(f"{i+1}: {gene}\n")
    
    logger.info(f"Created concatenated FASTA: {fasta_file}")
    logger.info(f"Gene order info: {info_file}")
    
    return fasta_file

def run_mafft(fasta_file: str, output_dir: str) -> str:
    """Run MAFFT alignment on concatenated FASTA file."""
    aligned_file = os.path.join(output_dir, "concatenated_genes_aligned.fasta")
    
    cmd = ['mafft', '--auto', fasta_file]
    
    try:
        with open(aligned_file, 'w') as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
        logger.info(f"MAFFT alignment completed: {aligned_file}")
        return aligned_file
    except subprocess.CalledProcessError as e:
        logger.error(f"MAFFT failed: {e}")
        return None

def run_trimal(aligned_file: str, output_dir: str) -> str:
    """Trim MAFFT alignment using trimAl with a 50% gap threshold."""
    trimmed_file = os.path.join(output_dir, "concatenated_genes_aligned_trimmed.fasta")
    cmd = ['trimal', '-in', aligned_file, '-out', trimmed_file, '-gt', '0.5']

    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logger.info(f"trimAl completed: {trimmed_file}")
        return trimmed_file
    except subprocess.CalledProcessError as e:
        logger.error(f"trimAl failed: {e}")
        return None

def run_fasttree(aligned_file: str, output_dir: str) -> str:
    """Run FastTree on aligned concatenated sequences."""
    tree_file = os.path.join(output_dir, "eukaryotic_phylogeny.tree")
    
    cmd = ['fasttree', '-out', tree_file, aligned_file]
    
    try:
        subprocess.run(cmd, check=True)
        logger.info(f"FastTree completed: {tree_file}")
        return tree_file
    except subprocess.CalledProcessError as e:
        logger.error(f"FastTree failed: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(
        description='Build phylogenetic trees from concatenated eukaryotic single copy genes'
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
        'fasta_dir',
        help='Directory containing the .fas files for each sample'
    )
    parser.add_argument(
        '--coverage-threshold',
        type=float,
        default=100.0,
        help='Minimum percentage of genomes that must have at least one gene (0-100, default: 100.0)'
    )
    parser.add_argument(
        '--evalue-cutoff',
        type=float,
        default=0.001,
        help='E-value cutoff for filtering gene hits (default: 0.001)'
    )
    parser.add_argument(
        '--output-dir',
        default='eukaryotic_trees',
        help='Output directory for trees and alignments (default: eukaryotic_trees)'
    )
    parser.add_argument(
        '--genome-list',
        help='File containing list of genome IDs to include (one per line, # for comments)'
    )
    
    args = parser.parse_args()
    
    # Validate coverage threshold
    if args.coverage_threshold < 0 or args.coverage_threshold > 100:
        logger.error("coverage-threshold must be between 0 and 100")
        sys.exit(1)
    
    try:
        # Load genome list if provided
        genome_list = load_genome_list(args.genome_list)
        if genome_list:
            logger.info(f"Filtering to {len(genome_list)} specified genomes")
        
        # Load single copy genes
        logger.info("Loading single copy genes...")
        gene_families = load_single_copy_genes(args.single_copy_csv)
        
        # Load ORF summary
        logger.info("Loading ORF summary...")
        samples_data = load_orf_summary(args.orf_summary, genome_list)
        
        # Check overall genome coverage across all genes
        logger.info("Checking overall genome coverage...")
        all_genes, genome_gene_coverage = check_overall_genome_coverage(gene_families, samples_data, args.coverage_threshold)
        
        if not all_genes:
            logger.warning("No genes meet the coverage threshold")
            sys.exit(0)
        
        # Extract and concatenate sequences
        logger.info("Extracting and concatenating sequences...")
        concatenated_sequences = extract_and_concatenate_sequences(
            all_genes, genome_gene_coverage, samples_data, args.fasta_dir, args.evalue_cutoff
        )
        
        # Create concatenated FASTA file
        logger.info("Creating concatenated FASTA file...")
        fasta_file = create_concatenated_fasta(concatenated_sequences, args.output_dir)
        
        # Run MAFFT on concatenated sequences
        logger.info("Running MAFFT alignment...")
        aligned_file = run_mafft(fasta_file, args.output_dir)

        if aligned_file:
            logger.info("Trimming alignment...")
            trimmed_file = run_trimal(aligned_file, args.output_dir)

            if trimmed_file:
                # Run FastTree on trimmed concatenated sequences
                logger.info("Running FastTree...")
                tree_file = run_fasttree(trimmed_file, args.output_dir)

                if tree_file:
                    logger.info(f"Successfully created phylogenetic tree: {tree_file}")
                    logger.info(f"Tree building completed. Results in: {args.output_dir}")
                else:
                    logger.error("FastTree failed")
            else:
                logger.error("trimAl failed")
        else:
            logger.error("MAFFT alignment failed")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
