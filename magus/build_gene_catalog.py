#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import subprocess
import tempfile
import shutil
from pathlib import Path
from collections import defaultdict
import pandas as pd

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GeneCatalogBuilder:
    def __init__(self, call_orfs_output_dir, output_dir, threads=1, 
                 evalue_cutoff=0.01, identity_threshold=0.9, coverage_threshold=0.8):
        """
        Initialize the gene catalog builder.
        
        Args:
            call_orfs_output_dir: Directory containing call-orfs output (annotations and fasta files)
            output_dir: Output directory for the gene catalog
            threads: Number of threads for MMseqs2
            evalue_cutoff: E-value cutoff for HMM annotations
            identity_threshold: Identity threshold for MMseqs2 clustering
            coverage_threshold: Coverage threshold for MMseqs2 clustering
        """
        self.call_orfs_output_dir = Path(call_orfs_output_dir)
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.evalue_cutoff = evalue_cutoff
        self.identity_threshold = identity_threshold
        self.coverage_threshold = coverage_threshold
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if MMseqs2 is available
        self._check_mmseqs2()
    
    def _check_mmseqs2(self):
        """Check if MMseqs2 is available in PATH."""
        try:
            subprocess.run(['mmseqs', 'version'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error("MMseqs2 not found in PATH. Please install MMseqs2.")
            sys.exit(1)
    
    def _get_annotation_files(self):
        """Get all annotation files from call-orfs output."""
        annotation_files = []
        
        # Look for summary files in the main output directory
        for summary_file in self.call_orfs_output_dir.glob("*_orf_summary*.tsv"):
            annotation_files.append(summary_file)
        
        # Also look in subdirectories for individual annotation files
        for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
            annot_dir = self.call_orfs_output_dir / subdir / 'annot'
            if annot_dir.exists():
                for hmm_file in annot_dir.glob("*.hmm.tsv"):
                    annotation_files.append(hmm_file)
        
        return annotation_files
    
    def _parse_annotations(self):
        """Parse HMM annotations and create gene-to-annotation mapping."""
        # Dictionary to store all annotations per gene: {gene_id: [(evalue, annotation), ...]}
        gene_all_annotations = defaultdict(list)
        annotated_genes = set()
        unannotated_genes = set()
        
        annotation_files = self._get_annotation_files()
        
        for annotation_file in annotation_files:
            logger.info(f"Processing annotation file: {annotation_file}")
            
            try:
                if annotation_file.name.endswith('_orf_summary.tsv'):
                    # This is a comprehensive summary file
                    df = pd.read_csv(annotation_file, sep='\t')
                    
                    # Check if this is the eukaryotic format or bacterial/viral format
                    if 'target_name' in df.columns:
                        # Eukaryotic format
                        for _, row in df.iterrows():
                            sample_id = row['sample_id']
                            gene_id = row['target_name']
                            query_name = row.get('query_name', '')
                            full_evalue = row.get('full_evalue', '')
                            
                            if pd.notna(query_name) and query_name != '' and pd.notna(full_evalue):
                                try:
                                    evalue = float(full_evalue)
                                    if evalue <= self.evalue_cutoff:
                                        full_gene_id = f"{sample_id}-----{gene_id}"
                                        gene_all_annotations[full_gene_id].append((evalue, query_name))
                                        annotated_genes.add(full_gene_id)
                                    else:
                                        full_gene_id = f"{sample_id}-----{gene_id}"
                                        unannotated_genes.add(full_gene_id)
                                except (ValueError, TypeError):
                                    # If evalue can't be parsed, treat as unannotated
                                    full_gene_id = f"{sample_id}-----{gene_id}"
                                    unannotated_genes.add(full_gene_id)
                            else:
                                full_gene_id = f"{sample_id}-----{gene_id}"
                                unannotated_genes.add(full_gene_id)
                    else:
                        # Bacterial/viral format - need to check for HMM data
                        # This would require additional parsing logic
                        logger.warning(f"Skipping {annotation_file} - bacterial/viral format not yet implemented")
                        continue
                        
                elif annotation_file.name.endswith('.hmm.tsv'):
                    # This is an individual HMM annotation file
                    sample_id = annotation_file.stem.replace('.hmm', '')
                    
                    with open(annotation_file, 'r') as f:
                        for line in f:
                            if line.startswith('#'):
                                continue
                            parts = line.strip().split()
                            if len(parts) >= 5:
                                query_name = parts[0]  # HMM name
                                gene_id = parts[2]    # Gene name
                                try:
                                    evalue = float(parts[4])
                                    if evalue <= self.evalue_cutoff:
                                        full_gene_id = f"{sample_id}-----{gene_id}"
                                        gene_all_annotations[full_gene_id].append((evalue, query_name))
                                        annotated_genes.add(full_gene_id)
                                    else:
                                        full_gene_id = f"{sample_id}-----{gene_id}"
                                        unannotated_genes.add(full_gene_id)
                                except (ValueError, TypeError):
                                    full_gene_id = f"{sample_id}-----{gene_id}"
                                    unannotated_genes.add(full_gene_id)
                            else:
                                full_gene_id = f"{sample_id}-----{gene_id}"
                                unannotated_genes.add(full_gene_id)
            
            except Exception as e:
                logger.warning(f"Error processing {annotation_file}: {e}")
                continue
        
        # Now find the BEST annotation per gene (lowest e-value)
        gene_annotations = {}
        for gene_id, annotations in gene_all_annotations.items():
            if annotations:
                # Sort by e-value (ascending) and take the first (lowest e-value)
                best_annotation = min(annotations, key=lambda x: x[0])
                gene_annotations[gene_id] = best_annotation[1]
        
        logger.info(f"Found {len(annotated_genes)} annotated genes and {len(unannotated_genes)} unannotated genes")
        return gene_annotations, annotated_genes, unannotated_genes
    
    def _get_fasta_files(self):
        """Get all fasta files from the call-orfs output directory."""
        fasta_files = []
        
        # Look for .faa files (protein sequences) in manicure and annot directories
        for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
            manicure_dir = self.call_orfs_output_dir / subdir / 'manicure'
            annot_dir = self.call_orfs_output_dir / subdir / 'annot'
            
            if manicure_dir.exists():
                for fasta_file in manicure_dir.glob("*.faa"):
                    fasta_files.append(fasta_file)
            
            if annot_dir.exists():
                for fasta_file in annot_dir.glob("*.fas"):
                    if not fasta_file.name.endswith('.codon.fas'):
                        fasta_files.append(fasta_file)
        
        return fasta_files
    
    def _create_unannotated_fasta(self, unannotated_genes, fasta_files):
        """Create a fasta file containing only unannotated genes."""
        unannotated_fasta = self.output_dir / "unannotated_genes.faa"
        
        gene_to_file = {}
        for fasta_file in fasta_files:
            sample_id = fasta_file.parent.name if fasta_file.parent.name in ['manicure', 'annot'] else fasta_file.stem
            
            with open(fasta_file, 'r') as f:
                current_gene = None
                current_sequence = []
                
                for line in f:
                    if line.startswith('>'):
                        if current_gene and current_sequence:
                            gene_to_file[current_gene] = (current_gene, ''.join(current_sequence))
                        
                        current_gene = line.strip()[1:]  # Remove '>'
                        current_sequence = []
                    else:
                        current_sequence.append(line.strip())
                
                if current_gene and current_sequence:
                    gene_to_file[current_gene] = (current_gene, ''.join(current_sequence))
        
        # Write unannotated genes to fasta
        with open(unannotated_fasta, 'w') as out_f:
            for gene_id in unannotated_genes:
                if gene_id in gene_to_file:
                    header, sequence = gene_to_file[gene_id]
                    out_f.write(f">{header}\n{sequence}\n")
        
        return unannotated_fasta
    
    def _cluster_unannotated_genes(self, unannotated_fasta):
        """Cluster unannotated genes using MMseqs2 easy-cluster."""
        if not unannotated_fasta.exists() or unannotated_fasta.stat().st_size == 0:
            logger.info("No unannotated genes to cluster")
            return {}
        
        logger.info("Clustering unannotated genes using MMseqs2")
        
        # Create temporary directory for MMseqs2
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            
            # Run MMseqs2 easy-cluster
            cluster_prefix = temp_dir / "cluster"
            cmd = [
                'mmseqs', 'easy-cluster',
                str(unannotated_fasta),
                str(cluster_prefix),
                temp_dir / "tmp",
                '--min-seq-id', str(self.identity_threshold),
                '--c', str(self.coverage_threshold),
                '--threads', str(self.threads)
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"MMseqs2 clustering failed: {e}")
                return {}
            
            # Parse cluster results
            cluster_file = cluster_prefix.with_suffix('.cluster')
            if not cluster_file.exists():
                logger.warning("No cluster file generated by MMseqs2")
                return {}
            
            gene_clusters = {}
            with open(cluster_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        representative = parts[0]
                        gene_id = parts[1]
                        gene_clusters[gene_id] = representative
            
            logger.info(f"Created {len(set(gene_clusters.values()))} clusters from unannotated genes")
            return gene_clusters
    
    def _create_gene_mappings(self, gene_annotations, annotated_genes, unannotated_genes, gene_clusters):
        """Create the final gene mapping files."""
        # Main cluster mapping
        cluster_mapping = self.output_dir / "gene_clusters.tsv"
        singleton_mapping = self.output_dir / "gene_singletons.tsv"
        
        # Track which genes have been processed
        processed_genes = set()
        
        with open(cluster_mapping, 'w') as cluster_f, open(singleton_mapping, 'w') as singleton_f:
            # Write headers
            cluster_f.write("sampleid\tgene\tclusterid\n")
            singleton_f.write("sampleid\tgene\tclusterid\n")
            
            # Process annotated genes
            for gene_id in annotated_genes:
                if gene_id in gene_annotations:
                    sample_id, gene = gene_id.split('-----', 1)
                    cluster_id = gene_annotations[gene_id]
                    cluster_f.write(f"{sample_id}\t{gene}\t{cluster_id}\n")
                    processed_genes.add(gene_id)
            
            # Process unannotated genes
            for gene_id in unannotated_genes:
                if gene_id in gene_clusters:
                    sample_id, gene = gene_id.split('-----', 1)
                    cluster_id = gene_clusters[gene_id]
                    cluster_f.write(f"{sample_id}\t{gene}\t{cluster_id}\n")
                    processed_genes.add(gene_id)
                else:
                    # This is a singleton
                    sample_id, gene = gene_id.split('-----', 1)
                    singleton_f.write(f"{sample_id}\t{gene}\t{sample_id}\n")
                    processed_genes.add(gene_id)
        
        logger.info(f"Created gene cluster mapping: {cluster_mapping}")
        logger.info(f"Created gene singleton mapping: {singleton_mapping}")
        
        # Count statistics
        cluster_count = len(set(gene_annotations.values())) + len(set(gene_clusters.values()))
        singleton_count = len(unannotated_genes - set(gene_clusters.keys()))
        
        logger.info(f"Total clusters: {cluster_count}")
        logger.info(f"Total singletons: {singleton_count}")
        
        return cluster_mapping, singleton_mapping
    
    def build_catalog(self):
        """Main method to build the gene catalog."""
        logger.info("Starting gene catalog construction")
        
        # Parse annotations
        gene_annotations, annotated_genes, unannotated_genes = self._parse_annotations()
        
        # Get fasta files
        fasta_files = self._get_fasta_files()
        
        # Create unannotated genes fasta
        unannotated_fasta = self._create_unannotated_fasta(unannotated_genes, fasta_files)
        
        # Cluster unannotated genes
        gene_clusters = self._cluster_unannotated_genes(unannotated_fasta)
        
        # Create final mappings
        cluster_mapping, singleton_mapping = self._create_gene_mappings(
            gene_annotations, annotated_genes, unannotated_genes, gene_clusters
        )
        
        logger.info("Gene catalog construction completed successfully")
        return cluster_mapping, singleton_mapping

def main():
    parser = argparse.ArgumentParser(
        description="Build gene catalog by clustering genes based on annotations and sequence similarity"
    )
    
    parser.add_argument(
        '--call-orfs-output',
        required=True,
        help='Directory containing call-orfs output (annotations and fasta files)'
    )
    

    
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for the gene catalog'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Number of threads for MMseqs2 (default: 1)'
    )
    
    parser.add_argument(
        '--evalue-cutoff',
        type=float,
        default=0.01,
        help='E-value cutoff for HMM annotations (default: 0.01)'
    )
    
    parser.add_argument(
        '--identity-threshold',
        type=float,
        default=0.9,
        help='Identity threshold for MMseqs2 clustering (default: 0.9)'
    )
    
    parser.add_argument(
        '--coverage-threshold',
        type=float,
        default=0.8,
        help='Coverage threshold for MMseqs2 clustering (default: 0.8)'
    )
    
    args = parser.parse_args()
    
    # Create and run the gene catalog builder
    builder = GeneCatalogBuilder(
        call_orfs_output_dir=args.call_orfs_output,
        output_dir=args.output_dir,
        threads=args.threads,
        evalue_cutoff=args.evalue_cutoff,
        identity_threshold=args.identity_threshold,
        coverage_threshold=args.coverage_threshold
    )
    
    try:
        builder.build_catalog()
    except Exception as e:
        logger.error(f"Gene catalog construction failed: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
