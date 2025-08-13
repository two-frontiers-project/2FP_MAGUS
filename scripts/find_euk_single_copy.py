#!/usr/bin/env python3
"""
Find eukaryotic single copy genes from ORF summary output.

This script identifies candidate single copy genes based on specific gene families
and filters for full sequence hits with lowest E-values per gene.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict, Counter
from typing import Dict, List, Set, Tuple

def load_orf_summary(summary_file: str) -> Dict[str, List[Dict]]:
    """Load ORF summary file and group by sample_id."""
    if not os.path.exists(summary_file):
        raise FileNotFoundError(f"Summary file not found: {summary_file}")
    
    samples_data = defaultdict(list)
    
    with open(summary_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sample_id = row['sample_id']
            samples_data[sample_id].append(row)
    
    print(f"Loaded {len(samples_data)} samples from {summary_file}")
    return samples_data

def filter_hmm_hits(samples_data: Dict[str, List[Dict]], 
                    fullseq_evalue_cutoff: float = 0.001) -> Dict[str, List[Dict]]:
    """Filter for full sequence hits and keep lowest E-value hit per gene per sample."""
    filtered_data = {}
    
    for sample_id, orfs in samples_data.items():
        # Group by query_name (gene family) and keep lowest E-value hit
        gene_best_hits = {}
        
        for orf in orfs:
            # Skip if no HMM hit
            if not orf.get('query_name') or orf.get('query_name') == '':
                continue
            
            query_name = orf['query_name']
            full_evalue = orf.get('full_evalue', '')
            
            # Skip if no E-value or above cutoff
            if not full_evalue or full_evalue == '':
                continue
            
            try:
                evalue = float(full_evalue)
                if evalue > fullseq_evalue_cutoff:
                    continue
            except ValueError:
                continue
            
            # Keep best hit (lowest E-value) for this gene family in this sample
            if query_name not in gene_best_hits or evalue < gene_best_hits[query_name]['evalue']:
                gene_best_hits[query_name] = {
                    'query_name': query_name,
                    'target_name': orf['target_name'],
                    'full_evalue': evalue,
                    'sample_id': sample_id
                }
        
        filtered_data[sample_id] = list(gene_best_hits.values())
    
    return filtered_data

def get_single_copy_candidates(filtered_data: Dict[str, List[Dict]]) -> Dict[str, List[str]]:
    """Identify candidate single copy genes based on gene families."""
    
    # Define gene families of interest
    gene_families = {
        'RNA_pol_Rpb1': ['RNA_pol_Rpb1_1', 'RNA_pol_Rpb1_2', 'RNA_pol_Rpb1_3', 
                         'RNA_pol_Rpb1_4', 'RNA_pol_Rpb1_5', 'RNA_pol_Rpb1_6', 
                         'RNA_pol_Rpb1_7', 'Rpb1_R'],
        'RNA_pol_Rpb2': ['RNA_pol_Rpb2_1', 'RNA_pol_Rpb2_2', 'RNA_pol_Rpb2_3',
                         'RNA_pol_Rpb2_4', 'RNA_pol_Rpb2_5', 'RNA_pol_Rpb2_6',
                         'RNA_pol_Rpb2_7'],
        'RNA_pol_other': ['RNA_pol_Rpb4', 'RNA_pol_Rpb5_N', 'RNA_pol_Rpb5_C', 
                          'RNA_pol_Rpb6', 'RNA_pol_Rpb8', 'RNA_pol_RpbG'],
        'TBP': ['TBP', 'TBP-binding', 'TBP-TOTE'],
        'PCNA': ['PCNA'],
        'RFC': ['RFC1', 'RFC2', 'RFC3', 'RFC4', 'RFC5'],
        'CPSF': ['CPSF73-100_C', 'CPSF100_C', 'CPSF_A', 'RSLD_CPSF6'],
        'CSTF': ['CSTF1_dimer', 'CSTF2_hinge', 'CSTF_C', 'partial_CstF'],
        'RNase_PH': ['RNase_PH', 'RNase_PH_C'],
        'rRNA_tRNA': ['Fibrillarin', 'Nop1', 'GCD14']
    }
    
    candidates = defaultdict(list)
    
    # Count occurrences of each gene family across samples
    for sample_id, genes in filtered_data.items():
        for gene in genes:
            query_name = gene['query_name']
            
            # Find which family this gene belongs to
            for family_name, family_genes in gene_families.items():
                if any(family_gene in query_name for family_gene in family_genes):
                    candidates[family_name].append(query_name)
                    break
    
    return candidates

def analyze_single_copy_candidates(candidates: Dict[str, List[str]], 
                                 filtered_data: Dict[str, List[Dict]]) -> None:
    """Analyze and report single copy gene candidates."""
    
    print("\n=== SINGLE COPY GENE CANDIDATES ===\n")
    
    # Count total samples
    total_samples = len(filtered_data)
    print(f"Total samples analyzed: {total_samples}")
    
    for family_name, genes in candidates.items():
        print(f"\n{family_name}:")
        
        # Count occurrences of each gene in this family
        gene_counts = Counter(genes)
        
        for gene, count in sorted(gene_counts.items()):
            percentage = (count / total_samples) * 100
            print(f"  {gene}: {count}/{total_samples} samples ({percentage:.1f}%)")
            
            # Show which samples have this gene
            samples_with_gene = []
            for sample_id, sample_genes in filtered_data.items():
                if any(g['query_name'] == gene for g in sample_genes):
                    samples_with_gene.append(sample_id)
            
            if count <= 10:  # Only show sample names if count is small
                print(f"    Samples: {', '.join(samples_with_gene)}")
    
    # Summary statistics
    print(f"\n=== SUMMARY ===")
    print(f"Total samples: {total_samples}")
    
    for family_name, genes in candidates.items():
        unique_genes = set(genes)
        print(f"{family_name}: {len(unique_genes)} unique genes found")
        
        # Find genes that appear in most samples (potential single copy)
        gene_counts = Counter(genes)
        high_frequency_genes = [(gene, count) for gene, count in gene_counts.items() 
                               if count >= total_samples * 0.8]  # 80% or more samples
        
        if high_frequency_genes:
            print(f"  High frequency genes (≥80% samples):")
            for gene, count in sorted(high_frequency_genes, key=lambda x: x[1], reverse=True):
                print(f"    {gene}: {count}/{total_samples}")

def main():
    parser = argparse.ArgumentParser(
        description='Find eukaryotic single copy genes from ORF summary output'
    )
    parser.add_argument(
        'summary_file',
        help='Path to eukaryotic ORF summary file (e.g., eukaryotes_orf_summary.tsv)'
    )
    parser.add_argument(
        '--fullseq-evalue',
        type=float,
        default=0.001,
        help='Full sequence E-value cutoff (default: 0.001)'
    )
    
    args = parser.parse_args()
    
    try:
        # Load ORF summary data
        print(f"Loading ORF summary from: {args.summary_file}")
        samples_data = load_orf_summary(args.summary_file)
        
        # Filter for full sequence hits and keep best hit per gene per sample
        print(f"Filtering HMM hits (E-value cutoff: {args.fullseq_evalue})")
        filtered_data = filter_hmm_hits(samples_data, args.fullseq_evalue)
        
        # Find single copy gene candidates
        print("Identifying single copy gene candidates...")
        candidates = get_single_copy_candidates(filtered_data)
        
        # Analyze and report results
        analyze_single_copy_candidates(candidates, filtered_data)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main() 