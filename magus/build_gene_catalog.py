#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict
import pandas as pd

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GeneCatalogBuilder:
    def __init__(self, summary_file, faa_dir, output_dir, threads=1, 
                 evalue_cutoff=0.01, identity_threshold=0.3, coverage_threshold=0.8,
                 identity_only=False):
        """
        Initialize the gene catalog builder.
        
        Args:
            summary_file: Path to call_orfs summary file with annotation data
            faa_dir: Directory containing .faa files where filenames match IDs in summary
            output_dir: Output directory for the gene catalog
            threads: Number of threads for MMseqs2
            evalue_cutoff: E-value cutoff for HMM annotations
            identity_threshold: Identity threshold for MMseqs2 clustering
            coverage_threshold: Coverage threshold for MMseqs2 clustering
            identity_only: If True, skip functional annotation and only use MMseqs2 clustering
        """
        self.summary_file = Path(summary_file)
        self.faa_dir = Path(faa_dir)
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.evalue_cutoff = evalue_cutoff
        self.identity_threshold = identity_threshold
        self.coverage_threshold = coverage_threshold
        self.identity_only = identity_only
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if MMseqs2 is available
        self._check_mmseqs2()
        
        # Store gene sequences
        self.gene_sequences = {}
    
    def _check_mmseqs2(self):
        """Check if MMseqs2 is available in PATH."""
        try:
            subprocess.run(['mmseqs', 'version'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error("MMseqs2 not found in PATH. Please install MMseqs2.")
            sys.exit(1)
    
    def _parse_summary_file(self):
        """Parse the call_orfs summary file to get annotated vs unannotated genes."""
        logger.info(f"Parsing summary file: {self.summary_file}")
        
        # Read the summary file and store as instance variable
        self.summary_df = pd.read_csv(self.summary_file, sep='\t')
        df = self.summary_df
        
        # Determine the format based on columns
        if 'target_name' in df.columns and 'query_name' in df.columns:
            # Eukaryotic format
            annotated_genes = set()
            unannotated_genes = set()
            gene_annotations = {}
            
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
                            gene_annotations[full_gene_id] = query_name
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
        
        elif 'sequence_id' in df.columns and ('query_name' in df.columns or 'query_accession' in df.columns):
            # Bacterial/viral format
            annotated_genes = set()
            unannotated_genes = set()
            gene_annotations = {}
            self.evalue_mapping = {}  # Build efficient E-value lookup
            
            # Determine which column to use for query name
            query_col = 'query_name' if 'query_name' in df.columns else 'query_accession'
            
            for _, row in df.iterrows():
                sample_id = row['sample_id']
                gene_id = row['sequence_id']
                query_name = row.get(query_col, '')
                full_evalue = row.get('full_evalue', '')
                
                # Store E-value for efficient lookup (if available)
                if pd.notna(full_evalue):
                    try:
                        evalue = float(full_evalue)
                        self.evalue_mapping[(sample_id, gene_id)] = evalue
                    except (ValueError, TypeError):
                        self.evalue_mapping[(sample_id, gene_id)] = float('inf')
                else:
                    self.evalue_mapping[(sample_id, gene_id)] = float('inf')
                
                if pd.notna(query_name) and query_name != '' and pd.notna(full_evalue):
                    try:
                        evalue = float(full_evalue)
                        if evalue <= self.evalue_cutoff:
                            full_gene_id = f"{sample_id}-----{gene_id}"
                            gene_annotations[full_gene_id] = query_name
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
        
        else:
            logger.error("Unrecognized summary file format")
            sys.exit(1)
        
        logger.info(f"Found {len(annotated_genes)} annotated genes and {len(unannotated_genes)} unannotated genes")
        return gene_annotations, annotated_genes, unannotated_genes
    
    def _load_faa_sequences(self):
        """Load all .faa sequences from the specified directory."""
        logger.info(f"Loading .faa sequences from: {self.faa_dir}")
        
        for faa_file in self.faa_dir.glob("*.faa"):
            sample_id = faa_file.stem  # filename without .faa extension
            
            with open(faa_file, 'r') as f:
                current_gene = None
                current_sequence = []
                
                for line in f:
                    if line.startswith('>'):
                        if current_gene and current_sequence:
                            # Store the full gene ID and sequence
                            full_gene_id = f"{sample_id}-----{current_gene}"
                            self.gene_sequences[full_gene_id] = ''.join(current_sequence)
                        
                        # Extract just the gene ID part (before any # or space)
                        header = line.strip()[1:]  # Remove '>'
                        current_gene = header.split('#')[0].split()[0]  # Take first part before # or space
                        current_sequence = []
                    else:
                        current_sequence.append(line.strip())
                
                if current_gene and current_sequence:
                    full_gene_id = f"{sample_id}-----{current_gene}"
                    self.gene_sequences[full_gene_id] = ''.join(current_sequence)
        
        logger.info(f"Loaded {len(self.gene_sequences)} gene sequences")
    
    def _create_unannotated_fasta(self, unannotated_genes):
        """Create a fasta file containing only unannotated genes."""
        unannotated_fasta = self.output_dir / "unannotated_genes.faa"
        
        with open(unannotated_fasta, 'w') as out_f:
            for gene_id in unannotated_genes:
                if gene_id in self.gene_sequences:
                    sequence = self.gene_sequences[gene_id]
                    out_f.write(f">{gene_id}\n{sequence}\n")
        
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
    
    def _cluster_all_genes_identity_only(self):
        """Cluster all genes using MMseqs2 when identity_only=True."""
        logger.info("Running identity-only clustering for all genes")
        
        # Create a combined FASTA with all genes
        all_genes_fasta = self.output_dir / "all_genes.faa"
        
        with open(all_genes_fasta, 'w') as out_f:
            for gene_id, sequence in self.gene_sequences.items():
                out_f.write(f">{gene_id}\n{sequence}\n")
        
        # Create temporary directory for MMseqs2
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            
            # Run MMseqs2 easy-cluster
            cluster_prefix = temp_dir / "cluster"
            cmd = [
                'mmseqs', 'easy-cluster',
                str(all_genes_fasta),
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
            
            logger.info(f"Created {len(set(gene_clusters.values()))} clusters from all genes")
            return gene_clusters
    
    def _create_gene_mappings(self, gene_annotations, annotated_genes, unannotated_genes, gene_clusters):
        """Create the final gene mapping files."""
        # Main cluster mapping
        cluster_mapping = self.output_dir / "gene_clusters.tsv"
        singleton_mapping = self.output_dir / "gene_singletons.tsv"
        
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
            
            # Process unannotated genes
            for gene_id in unannotated_genes:
                if gene_id in gene_clusters:
                    sample_id, gene = gene_id.split('-----', 1)
                    cluster_id = gene_clusters[gene_id]
                    cluster_f.write(f"{sample_id}\t{gene}\t{cluster_id}\n")
                else:
                    # This is a singleton
                    sample_id, gene = gene_id.split('-----', 1)
                    singleton_f.write(f"{sample_id}\t{gene}\t{sample_id}\n")
        
        logger.info(f"Created gene cluster mapping: {cluster_mapping}")
        logger.info(f"Created gene singleton mapping: {singleton_mapping}")
        
        return cluster_mapping, singleton_mapping
    
    def _create_identity_only_mappings(self, gene_clusters):
        """Create gene mappings for identity-only mode."""
        # Main cluster mapping
        cluster_mapping = self.output_dir / "gene_clusters.tsv"
        singleton_mapping = self.output_dir / "gene_singletons.tsv"
        
        with open(cluster_mapping, 'w') as cluster_f, open(singleton_mapping, 'w') as singleton_f:
            # Write headers
            cluster_f.write("sampleid\tgene\tclusterid\n")
            singleton_f.write("sampleid\tgene\tclusterid\n")
            
            # Process all genes
            all_genes = set(self.gene_sequences.keys())
            clustered_genes = set(gene_clusters.keys())
            singleton_genes = all_genes - clustered_genes
            
            # Write clustered genes
            for gene_id in clustered_genes:
                sample_id, gene = gene_id.split('-----', 1)
                cluster_id = gene_clusters[gene_id]
                cluster_f.write(f"{sample_id}\t{gene}\t{cluster_id}\n")
            
            # Write singleton genes
            for gene_id in singleton_genes:
                sample_id, gene = gene_id.split('-----', 1)
                singleton_f.write(f"{sample_id}\t{gene}\t{sample_id}\n")
        
        logger.info(f"Created identity-only gene cluster mapping: {cluster_mapping}")
        logger.info(f"Created identity-only gene singleton mapping: {singleton_mapping}")
        
        return cluster_mapping, singleton_mapping
    
    def _create_identity_only_fastas(self, gene_clusters):
        """Create non-redundant FASTA files for identity-only mode."""
        logger.info("Creating identity-only non-redundant FASTA files")
        
        # Files to create
        cluster_faa = self.output_dir / "gene_catalog_clusters.faa"
        singleton_faa = self.output_dir / "gene_catalog_singletons.faa"
        
        # Group genes by their cluster
        cluster_groups = defaultdict(list)
        all_genes = set(self.gene_sequences.keys())
        clustered_genes = set(gene_clusters.keys())
        singleton_genes = all_genes - clustered_genes
        
        # Process clustered genes
        for gene_id in clustered_genes:
            cluster_id = gene_clusters[gene_id]
            cluster_groups[cluster_id].append(gene_id)
        
        # Write cluster FASTA with representative sequences
        with open(cluster_faa, 'w') as out_f:
            for cluster_id, gene_list in cluster_groups.items():
                if gene_list:
                    # Choose longest sequence as representative
                    representative_gene = max(gene_list, key=lambda g: len(self.gene_sequences.get(g, '')))
                    sequence = self.gene_sequences.get(representative_gene, '')
                    if sequence:
                        out_f.write(f">{cluster_id}\n{sequence}\n")
        
        # Write singleton FASTA
        with open(singleton_faa, 'w') as out_f:
            for gene_id in singleton_genes:
                sequence = self.gene_sequences.get(gene_id, '')
                if sequence:
                    out_f.write(f">{gene_id}\n{sequence}\n")
        
        logger.info(f"Created identity-only cluster FASTA: {cluster_faa}")
        logger.info(f"Created identity-only singleton FASTA: {singleton_faa}")
        
        return cluster_faa, singleton_faa
    
    def _create_nonredundant_fastas(self, gene_annotations, annotated_genes, unannotated_genes, gene_clusters):
        """Create non-redundant FASTA files for clusters and singletons."""
        logger.info("Creating non-redundant FASTA files")
        
        # Files to create
        cluster_faa = self.output_dir / "gene_catalog_clusters.faa"
        singleton_faa = self.output_dir / "gene_catalog_singletons.faa"
        
        # Group genes by their cluster/annotation
        cluster_groups = defaultdict(list)
        singleton_genes = set()
        
        # Process annotated genes (functional groups)
        for gene_id in annotated_genes:
            if gene_id in gene_annotations:
                cluster_id = gene_annotations[gene_id]
                cluster_groups[cluster_id].append(gene_id)
        
        # Process unannotated genes
        for gene_id in unannotated_genes:
            if gene_id in gene_clusters:
                cluster_id = gene_clusters[gene_id]
                cluster_groups[cluster_id].append(gene_id)
            else:
                singleton_genes.add(gene_id)
        
        # Write cluster FASTA with representative sequences
        with open(cluster_faa, 'w') as out_f:
            for cluster_id, gene_list in cluster_groups.items():
                if gene_list:
                    # For functional clusters, choose representative based on best E-value
                    # For sequence clusters, choose longest sequence as fallback
                    if cluster_id in [gene_annotations.get(g, '') for g in gene_list if g in gene_annotations]:
                        # This is a functional cluster - choose best E-value
                        best_gene = None
                        best_evalue = float('inf')
                        
                        for gene_id in gene_list:
                            if gene_id in gene_annotations:
                                # Use the pre-built mapping to get E-value efficiently
                                sample_id, seq_id = gene_id.split('-----', 1)
                                # Create the lookup key that matches the summary file format
                                lookup_key = (sample_id, seq_id)
                                
                                if lookup_key in self.evalue_mapping:
                                    evalue = self.evalue_mapping[lookup_key]
                                    if evalue < best_evalue:
                                        best_evalue = evalue
                                        best_gene = gene_id
                        
                        representative_gene = best_gene if best_gene else gene_list[0]
                    else:
                        # This is a sequence cluster - choose longest sequence
                        representative_gene = max(gene_list, key=lambda g: len(self.gene_sequences.get(g, '')))
                    
                    sequence = self.gene_sequences.get(representative_gene, '')
                    if sequence:
                        out_f.write(f">{cluster_id}\n{sequence}\n")
        
        # Write singleton FASTA
        with open(singleton_faa, 'w') as out_f:
            for gene_id in singleton_genes:
                sequence = self.gene_sequences.get(gene_id, '')
                if sequence:
                    out_f.write(f">{gene_id}\n{sequence}\n")
        
        logger.info(f"Created cluster FASTA: {cluster_faa}")
        logger.info(f"Created singleton FASTA: {singleton_faa}")
        
        return cluster_faa, singleton_faa
    
    def build_catalog(self):
        """Main method to build the gene catalog."""
        logger.info("Starting gene catalog construction")
        
        if self.identity_only:
            logger.info("Running in identity-only mode - skipping functional annotation")
            # Load .faa sequences
            self._load_faa_sequences()
            
            # Cluster all genes using MMseqs2
            gene_clusters = self._cluster_all_genes_identity_only()
            
            # Create identity-only mappings
            cluster_mapping, singleton_mapping = self._create_identity_only_mappings(gene_clusters)
            
            # Create non-redundant FASTA files
            cluster_faa, singleton_faa = self._create_identity_only_fastas(gene_clusters)
            
            logger.info("Identity-only gene catalog construction completed successfully")
            return cluster_mapping, singleton_mapping, cluster_faa, singleton_faa
        else:
            # Standard functional annotation mode
            # Parse summary file
            gene_annotations, annotated_genes, unannotated_genes = self._parse_summary_file()
            
            # Load .faa sequences
            self._load_faa_sequences()
            
            # Create unannotated genes fasta
            unannotated_fasta = self._create_unannotated_fasta(unannotated_genes)
            
            # Cluster unannotated genes
            gene_clusters = self._cluster_unannotated_genes(unannotated_fasta)
            
            # Create final mappings
            cluster_mapping, singleton_mapping = self._create_gene_mappings(
                gene_annotations, annotated_genes, unannotated_genes, gene_clusters
            )
            
            # Create non-redundant FASTA files
            cluster_faa, singleton_faa = self._create_nonredundant_fastas(
                gene_annotations, annotated_genes, unannotated_genes, gene_clusters
            )
            
            logger.info("Gene catalog construction completed successfully")
            return cluster_mapping, singleton_mapping, cluster_faa, singleton_faa

def main():
    parser = argparse.ArgumentParser(
        description="Build gene catalog by clustering genes based on annotations and sequence similarity"
    )
    
    parser.add_argument(
        '--summary-file',
        required=True,
        help='Path to call_orfs summary file with annotation data'
    )
    
    parser.add_argument(
        '--faa-dir',
        required=True,
        help='Directory containing .faa files where filenames match IDs in summary'
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
    
    parser.add_argument(
        '--identity-only',
        action='store_true',
        help='Skip functional annotation and only use MMseqs2 sequence clustering'
    )
    
    args = parser.parse_args()
    
    # Create and run the gene catalog builder
    builder = GeneCatalogBuilder(
        summary_file=args.summary_file,
        faa_dir=args.faa_dir,
        output_dir=args.output_dir,
        threads=args.threads,
        evalue_cutoff=args.evalue_cutoff,
        identity_threshold=args.identity_threshold,
        coverage_threshold=args.coverage_threshold,
        identity_only=args.identity_only
    )
    
    try:
        builder.build_catalog()
    except Exception as e:
        logger.error(f"Gene catalog construction failed: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()