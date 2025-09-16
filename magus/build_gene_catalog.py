#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
except ImportError:
    print("ERROR: Biopython is required but not installed.")
    print("Please install it with: pip install biopython")
    sys.exit(1)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GeneCatalogBuilder:
    def __init__(self, summary_file, faa_dir, output_dir, threads=1, 
                 evalue_cutoff=0.01, identity_threshold=0.3, coverage_threshold=0.8,
                 identity_only=False, multi_sample_single_copy=False, tmp_dir='./tmp/'):
        self.summary_file = Path(summary_file) if summary_file else None
        self.faa_dir = Path(faa_dir)
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.evalue_cutoff = evalue_cutoff
        self.identity_threshold = identity_threshold
        self.coverage_threshold = coverage_threshold
        self.identity_only = identity_only
        self.multi_sample_single_copy = multi_sample_single_copy
        self.tmp_dir = Path(tmp_dir)
        
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
    
    def _parse_summary_file(self):
        """Parse the call_orfs summary file to get annotated vs unannotated genes."""
        logger.info(f"Parsing summary file: {self.summary_file}")
        
        annotated_genes = set()
        unannotated_genes = set()
        gene_annotations = {}
        self.evalue_mapping = {}
        
        with open(self.summary_file, 'r') as f:
            # Read header to determine format
            header_line = f.readline().strip()
            if not header_line:
                logger.error("Empty summary file")
                sys.exit(1)
            
            headers = header_line.split('\t')
            
            # Determine the format based on columns
            if 'target_name' in headers and 'query_name' in headers:
                # Eukaryotic format
                for line_num, line in enumerate(f, start=2):
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    # Pad with empty strings if line has fewer columns than headers
                    while len(parts) < len(headers):
                        parts.append('')
                    
                    row = dict(zip(headers, parts))
                    sample_id = row['sample_id']
                    gene_id = row['target_name']
                    query_name = row.get('query_name', '')
                    full_evalue = row.get('full_evalue', '')
                    
                    if query_name and query_name.strip():
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
        
            elif 'sequence_id' in headers and ('query_name' in headers or 'query_accession' in headers):
                # Bacterial/viral format
                query_col = 'query_name' if 'query_name' in headers else 'query_accession'
                
                for line_num, line in enumerate(f, start=2):
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    # Pad with empty strings if line has fewer columns than headers
                    while len(parts) < len(headers):
                        parts.append('')
                    
                    row = dict(zip(headers, parts))
                    sample_id = row['sample_id']
                    gene_id = row['sequence_id']
                    query_name = row.get(query_col, '')
                    full_evalue = row.get('full_evalue', '')
                    
                    # Store E-value for efficient lookup (if available)
                    if full_evalue and full_evalue.strip():
                        try:
                            evalue = float(full_evalue)
                            self.evalue_mapping[(sample_id, gene_id)] = evalue
                        except (ValueError, TypeError):
                            self.evalue_mapping[(sample_id, gene_id)] = float('inf')
                    else:
                        self.evalue_mapping[(sample_id, gene_id)] = float('inf')
                    
                    if query_name and query_name.strip():
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
    
    def _create_unannotated_fasta(self, unannotated_genes):
        """Create a fasta file containing only unannotated genes using Biopython."""
        if self.identity_only:
            # In identity-only mode, combine all FASTA files without filtering
            logger.info("Identity-only mode: combining all FASTA files")
            return self._combine_all_fasta_files()
        
        unannotated_fasta = self.output_dir / "unannotated_genes.faa"
        
        # Create a set for faster lookup
        unannotated_set = set(unannotated_genes)
        
        # Get unique sample IDs from the summary file
        sample_ids = set()
        for gene_id in unannotated_genes:
            sample_id = gene_id.split('-----')[0]
            sample_ids.add(sample_id)
        
        with open(unannotated_fasta, 'w') as out_f:
            for sample_id in sample_ids:
                faa_file = self.faa_dir / f"{sample_id}.faa"
                if not faa_file.exists():
                    logger.warning(f"FASTA file not found for sample {sample_id}: {faa_file}")
                    continue
                    
                try:
                    # Parse FASTA file with Biopython
                    for record in SeqIO.parse(faa_file, "fasta"):
                        # Extract sample_id from filename (remove .faa extension)
                        sample_id = faa_file.stem
                        # The record.id already contains the full format, extract just the basic gene ID
                        if '-----' in record.id:
                            # Format: sample_id-----gene_id-----additional_info
                            full_gene_id = record.id.split('-----', 2)[0] + '-----' + record.id.split('-----', 2)[1]
                        else:
                            # Fallback for basic format
                            basic_gene_id = record.id.split()[0] if ' ' in record.id else record.id
                            full_gene_id = f"{sample_id}-----{basic_gene_id}"
                        
                        # Check if this gene is unannotated
                        if full_gene_id in unannotated_set:
                            # Write the sequence to output
                            SeqIO.write(record, out_f, "fasta")
                            
                except Exception as e:
                    logger.error(f"Error processing {faa_file}: {e}")
                    continue
        return unannotated_fasta
    
    def _combine_all_fasta_files(self):
        """Combine all FASTA files into one for identity-only mode."""
        combined_fasta = self.output_dir / "tocluster.faa"
        
        with open(combined_fasta, 'w') as out_f:
            for faa_file in self.faa_dir.glob("*.faa"):
                try:
                    # Parse and write all sequences from this file
                    for record in SeqIO.parse(faa_file, "fasta"):
                        SeqIO.write(record, out_f, "fasta")
                    
                except Exception as e:
                    logger.error(f"Error processing {faa_file}: {e}")
                    continue
        return combined_fasta
    
    def _cluster_unannotated_genes(self, unannotated_fasta):
        """Cluster unannotated genes using MMseqs2 easy-cluster."""
        if not unannotated_fasta.exists() or unannotated_fasta.stat().st_size == 0:
            logger.info("No unannotated genes to cluster")
            return {}
        
        logger.info("Clustering unannotated genes using MMseqs2")
        
        # Create tmp directory if it doesn't exist
        self.tmp_dir.mkdir(parents=True, exist_ok=True)
        
        # Run MMseqs2 easy-cluster - output goes to output dir
        cluster_prefix = self.output_dir / "mmseqs_genecat"
        cmd = [
            'mmseqs', 'easy-cluster',
            str(unannotated_fasta),
            str(cluster_prefix),
            str(self.tmp_dir),
            '--min-seq-id', str(self.identity_threshold),
            '-c', str(self.coverage_threshold),
            '--threads', str(self.threads)
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"MMseqs2 ERROR: {e}")
            print(f"MMseqs2 STDERR: {e.stderr}")
            print(f"MMseqs2 STDOUT: {e.stdout}")
            logger.error(f"MMseqs2 clustering failed: {e}")
            logger.error(f"MMseqs2 stderr: {e.stderr}")
            logger.error(f"MMseqs2 stdout: {e.stdout}")
            return {}
        
        # Parse cluster results - MMseqs2 creates genecat_cluster.tsv
        cluster_file = cluster_prefix.parent / f"{cluster_prefix.name}_cluster.tsv"
        if not cluster_file.exists():
            logger.error("MMseqs2 failed to generate cluster file")
            logger.error(f"Expected cluster file: {cluster_file}")
            logger.error("This is a critical failure - MMseqs2 clustering is required for the pipeline")
            raise RuntimeError("MMseqs2 clustering failed - no cluster file generated")
        
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
        """Create the final gene mapping files with proper frequency-based clustering."""
        # Main cluster mapping
        cluster_mapping = self.output_dir / "gene_clusters.tsv"
        singleton_mapping = self.output_dir / "gene_singletons.tsv"
        single_copy_mapping = self.output_dir / "single_copy_genes.tsv"
        
        # Count annotation frequencies across samples
        annotation_frequency = defaultdict(int)
        annotation_samples = defaultdict(set)
        annotation_genes = defaultdict(list)
        
        # Process annotated genes to count frequencies
        for gene_id in annotated_genes:
            if gene_id in gene_annotations:
                sample_id, gene = gene_id.split('-----', 1)
                annotation = gene_annotations[gene_id]
                annotation_frequency[annotation] += 1
                annotation_samples[annotation].add(sample_id)
                annotation_genes[annotation].append(gene_id)
        
        # Separate annotations into clusters vs singletons based on frequency
        cluster_annotations = set()
        singleton_annotations = set()
        multi_sample_single_copy_annotations = set()
        
        for annotation, count in annotation_frequency.items():
            if count > 1:
                # Appears multiple times - goes to clusters
                cluster_annotations.add(annotation)
            else:
                # Appears only once - check if multi-sample single-copy
                if self.multi_sample_single_copy and len(annotation_samples[annotation]) > 1:
                    # Appears once per sample but across multiple samples
                    multi_sample_single_copy_annotations.add(annotation)
                    cluster_annotations.add(annotation)  # Also goes to clusters
                else:
                    # True singleton - appears only once total
                    singleton_annotations.add(annotation)
        
        with open(cluster_mapping, 'w') as cluster_f, open(singleton_mapping, 'w') as singleton_f, open(single_copy_mapping, 'w') as single_copy_f:
            # Write headers - 3 columns only
            cluster_f.write("sampleid\tgene\tcluster_rep\n")
            singleton_f.write("sampleid\tgene\tcluster_rep\n")
            single_copy_f.write("sampleid\tgene\tcluster_rep\n")
            
            # Process annotated genes
            for gene_id in annotated_genes:
                if gene_id in gene_annotations:
                    sample_id, gene = gene_id.split('-----', 1)
                    annotation = gene_annotations[gene_id]
                    
                    if annotation in cluster_annotations:
                        # Non-singleton functional annotation - goes to clusters
                        cluster_f.write(f"{sample_id}\t{gene}\t{annotation}\n")
                        if annotation in multi_sample_single_copy_annotations:
                            single_copy_f.write(f"{sample_id}\t{gene}\t{annotation}\n")
                    elif annotation in singleton_annotations:
                        # Singleton functional annotation - goes to singletons
                        singleton_f.write(f"{sample_id}\t{gene}\t{annotation}\n")
            
            # Process unannotated genes
            for gene_id in unannotated_genes:
                sample_id, gene = gene_id.split('-----', 1)
                if gene_id in gene_clusters:
                    # Non-singleton sequence cluster - goes to clusters
                    cluster_rep = gene_clusters[gene_id]
                    cluster_f.write(f"{sample_id}\t{gene}\t{cluster_rep}\n")
                else:
                    # Singleton sequence - goes to singletons
                    singleton_f.write(f"{sample_id}\t{gene}\t{gene_id}\n")
        
        logger.info(f"Created gene cluster mapping: {cluster_mapping}")
        logger.info(f"Created gene singleton mapping: {singleton_mapping}")
        if self.multi_sample_single_copy:
            logger.info(f"Created single copy genes mapping: {single_copy_mapping}")
        
        return cluster_mapping, singleton_mapping, single_copy_mapping
    
    def _create_identity_only_mappings(self, gene_clusters):
        """Create simple mappings for identity-only mode."""
        cluster_mapping = self.output_dir / "gene_clusters.tsv"
        singleton_mapping = self.output_dir / "gene_singletons.tsv"
        
        with open(cluster_mapping, 'w') as cluster_f, open(singleton_mapping, 'w') as singleton_f:
            # Write headers
            cluster_f.write("sampleid\tgene\tcluster_rep\n")
            singleton_f.write("sampleid\tgene\tcluster_rep\n")
            
            # Process all genes from FASTA files
            for faa_file in self.faa_dir.glob("*.faa"):
                sample_id = faa_file.stem
                
                try:
                    for record in SeqIO.parse(faa_file, "fasta"):
                        # Extract just the basic gene ID from the full format
                        if '-----' in record.id:
                            gene_id = record.id.split('-----', 2)[1]  # Get the middle part
                            full_gene_id = record.id.split('-----', 2)[0] + '-----' + record.id.split('-----', 2)[1]
                        else:
                            gene_id = record.id
                            full_gene_id = f"{sample_id}-----{gene_id}"
                        
                        if full_gene_id in gene_clusters:
                            # Non-singleton sequence cluster
                            cluster_rep = gene_clusters[full_gene_id]
                            cluster_f.write(f"{sample_id}\t{gene_id}\t{cluster_rep}\n")
                        else:
                            # Singleton sequence
                            singleton_f.write(f"{sample_id}\t{gene_id}\t{full_gene_id}\n")
                            
                except Exception as e:
                    logger.error(f"Error processing {faa_file}: {e}")
                    continue
        
        return cluster_mapping, singleton_mapping
    
    def build_catalog(self):
        """Main method to build the gene catalog."""
        logger.info("Starting gene catalog construction")
        
        if self.identity_only:
            # In identity-only mode, skip annotation parsing and just cluster all sequences
            logger.info("Running in identity-only mode - clustering all sequences")
            
            # Create combined fasta file
            combined_fasta = self._combine_all_fasta_files()
            
            # Cluster all sequences
            gene_clusters = self._cluster_unannotated_genes(combined_fasta)
            
            # Create simple mappings for identity-only mode
            cluster_mapping, singleton_mapping = self._create_identity_only_mappings(gene_clusters)
            
        else:
            # Parse summary file
            gene_annotations, annotated_genes, unannotated_genes = self._parse_summary_file()
            
            # Create unannotated genes fasta
            unannotated_fasta = self._create_unannotated_fasta(unannotated_genes)
            
            # Cluster unannotated genes
            gene_clusters = self._cluster_unannotated_genes(unannotated_fasta)
            
            # Create final mappings
            cluster_mapping, singleton_mapping, single_copy_mapping = self._create_gene_mappings(
                gene_annotations, annotated_genes, unannotated_genes, gene_clusters
            )
        
        logger.info("Gene catalog construction completed successfully")
        return cluster_mapping, singleton_mapping

def main():
    parser = argparse.ArgumentParser(
        description="Build gene catalog by clustering genes based on annotations and sequence similarity"
    )
    
    parser.add_argument(
        '--summary-file',
        required=False,
        help='Path to call_orfs summary file with annotation data (not required for identity-only mode)'
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
    
    parser.add_argument(
        '--multi-sample-single-copy',
        action='store_true',
        help='Identify genes that appear once per sample but across multiple samples'
    )
    
    parser.add_argument(
        '--tmp-dir',
        default='./tmp/',
        help='Temporary directory for MMseqs2 (default: ./tmp/)'
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if not args.identity_only and not args.summary_file:
        parser.error("--summary-file is required when not using --identity-only mode")
    
    # Create and run the gene catalog builder
    builder = GeneCatalogBuilder(
        summary_file=args.summary_file,
        faa_dir=args.faa_dir,
        output_dir=args.output_dir,
        threads=args.threads,
        evalue_cutoff=args.evalue_cutoff,
        identity_threshold=args.identity_threshold,
        coverage_threshold=args.coverage_threshold,
        identity_only=args.identity_only,
        multi_sample_single_copy=args.multi_sample_single_copy,
        tmp_dir=args.tmp_dir
    )
    
    try:
        builder.build_catalog()
    except Exception as e:
        logger.error(f"Gene catalog construction failed: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()