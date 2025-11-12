#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import subprocess
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GeneCatalogBuilder:
    def __init__(self, sequence_dir=None, sequence_file=None, output_dir=None, threads=1, 
                 identity_threshold=0.9, identity_thresholds=None, 
                 coverage_threshold=0.8, extension='faa', tmpdir='./tmp/', split_singletons=False):
        if sequence_dir and sequence_file:
            logger.error("Cannot specify both --sequence-dir and --sequence-file")
            sys.exit(1)
        if not sequence_dir and not sequence_file:
            logger.error("Must specify either --sequence-dir or --sequence-file")
            sys.exit(1)
        
        self.sequence_dir = Path(sequence_dir) if sequence_dir else None
        self.sequence_file = Path(sequence_file) if sequence_file else None
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.identity_threshold = identity_threshold
        self.identity_thresholds = identity_thresholds or []
        self.coverage_threshold = coverage_threshold
        self.extension = extension
        self.tmp_dir = Path(tmpdir)
        self.split_singletons = split_singletons
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._check_mmseqs2()
    
    def _check_mmseqs2(self):
        try:
            subprocess.run(['mmseqs', 'version'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error("MMseqs2 not found in PATH. Please install MMseqs2.")
            sys.exit(1)
    
    def _get_input_files(self):
        if self.sequence_file:
            if not self.sequence_file.is_file():
                logger.error(f"Sequence file does not exist: {self.sequence_file}")
                sys.exit(1)
            return [self.sequence_file]
        elif self.sequence_dir:
            if not self.sequence_dir.is_dir():
                logger.error(f"Sequence directory does not exist: {self.sequence_dir}")
                sys.exit(1)
            ext = f".{self.extension}" if not self.extension.startswith('.') else self.extension
            files = list(self.sequence_dir.glob(f"*{ext}"))
            if not files:
                logger.error(f"No {ext} files found in {self.sequence_dir}")
                sys.exit(1)
            return files
    
    def _combine_fasta_files(self, input_files):
        combined_fasta = self.output_dir / "tocluster.faa"
        with open(combined_fasta, 'w') as out_f:
            for faa_file in input_files:
                with open(faa_file, 'r') as in_f:
                    out_f.write(in_f.read())
        return combined_fasta
    
    def _cluster_with_threshold(self, input_fasta, min_seq_id, tag):
        self.tmp_dir.mkdir(parents=True, exist_ok=True)
        cluster_prefix = self.output_dir / f"mmseqs_iter_{tag}"
        cmd = [
            'mmseqs', 'easy-cluster',
            str(input_fasta),
            str(cluster_prefix),
            str(self.tmp_dir),
            '--min-seq-id', str(min_seq_id),
            '-c', str(self.coverage_threshold),
            '--threads', str(self.threads)
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"MMseqs2 clustering failed at threshold {min_seq_id}: {e.stderr}")
            raise
        cluster_file = cluster_prefix.parent / f"{cluster_prefix.name}_cluster.tsv"
        rep_fasta = cluster_prefix.parent / f"{cluster_prefix.name}_rep_seq.fasta"
        if not cluster_file.exists() or not rep_fasta.exists():
            raise RuntimeError(f"MMseqs2 outputs missing for tag {tag}")
        mapping = {}
        with open(cluster_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    representative = parts[0]
                    member = parts[1]
                    mapping[member] = representative
        return mapping, rep_fasta
    
    def _cluster_single_threshold(self, input_fasta):
        self.tmp_dir.mkdir(parents=True, exist_ok=True)
        cluster_prefix = self.output_dir / "mmseqs_genecat"
        cmd = [
            'mmseqs', 'easy-cluster',
            str(input_fasta),
            str(cluster_prefix),
            str(self.tmp_dir),
            '--min-seq-id', str(self.identity_threshold),
            '-c', str(self.coverage_threshold),
            '--threads', str(self.threads)
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"MMseqs2 clustering failed: {e.stderr}")
            raise
        cluster_file = cluster_prefix.parent / f"{cluster_prefix.name}_cluster.tsv"
        if not cluster_file.exists():
            raise RuntimeError("MMseqs2 clustering failed - no cluster file generated")
        gene_clusters = {}
        with open(cluster_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    representative = parts[0]
                    gene_id = parts[1]
                    gene_clusters[gene_id] = representative
        return gene_clusters
    
    def _iterative_cluster_chain(self, input_fasta):
        current_fasta = input_fasta
        all_mappings = []
        for idx, thr in enumerate(self.identity_thresholds, start=1):
            tag = f"t{idx}_{str(thr).replace('.','p')}"
            logger.info(f"Clustering at identity {thr}")
            mapping, rep_fasta = self._cluster_with_threshold(current_fasta, thr, tag)
            all_mappings.append((thr, mapping))
            current_fasta = rep_fasta
        return all_mappings
    
    def _extract_sample_gene(self, gene_id):
        if '-----' in gene_id:
            parts = gene_id.split('-----', 2)
            sample = parts[0]
            gene = parts[1] if len(parts) > 1 else gene_id
        else:
            sample = 'unknown'
            gene = gene_id
        return sample, gene
    
    def _create_gene_catalog(self, clustering_result, input_files, is_iterative=False):
        all_genes = {}
        
        for faa_file in input_files:
            with open(faa_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        gene_id = line.strip()[1:].split()[0]
                        sample, gene = self._extract_sample_gene(gene_id)
                        if gene_id not in all_genes:
                            all_genes[gene_id] = (sample, gene)
        
        final_catalog = self.output_dir / "gene_catalog.tsv"
        
        if is_iterative:
            threshold_mappings = clustering_result
            thresholds = [str(thr) for thr, _ in threshold_mappings]
            header = "sample\tgene\t" + "\t".join(thresholds) + "\n"
            
            with open(final_catalog, 'w') as out_f:
                out_f.write(header)
                for gene_id, (sample, gene) in all_genes.items():
                    row = [sample, gene]
                    current_id = gene_id
                    for thr, mapping in threshold_mappings:
                        if current_id in mapping:
                            current_id = mapping[current_id]
                        row.append(current_id)
                    out_f.write("\t".join(row) + "\n")
        else:
            gene_clusters = clustering_result
            header = f"sample\tgene\t{self.identity_threshold}\n"
            
            with open(final_catalog, 'w') as out_f:
                out_f.write(header)
                for gene_id, (sample, gene) in all_genes.items():
                    if gene_id in gene_clusters:
                        rep = gene_clusters[gene_id]
                    else:
                        rep = gene_id
                    out_f.write(f"{sample}\t{gene}\t{rep}\n")
        
        return final_catalog
    
    def _split_singletons(self, catalog_file):
        logger.info("Identifying and splitting singletons")
        
        with open(catalog_file, 'r') as f:
            header = f.readline().strip()
            lines = [line.strip() for line in f]
        
        header_parts = header.split('\t')
        if len(header_parts) < 3:
            logger.error("Invalid catalog format")
            return
        
        rep_columns = header_parts[2:]
        
        rep_counts = {}
        gene_rows = []
        
        for line in lines:
            parts = line.split('\t')
            if len(parts) < len(header_parts):
                continue
            gene_rows.append(parts)
            
            for idx, rep_col in enumerate(rep_columns):
                rep_val = parts[2 + idx]
                if rep_val not in rep_counts:
                    rep_counts[rep_val] = {}
                if idx not in rep_counts[rep_val]:
                    rep_counts[rep_val][idx] = 0
                rep_counts[rep_val][idx] += 1
        
        singleton_genes = set()
        for parts in gene_rows:
            is_singleton = True
            for idx in range(len(rep_columns)):
                rep_val = parts[2 + idx]
                if rep_counts[rep_val][idx] > 1:
                    is_singleton = False
                    break
            if is_singleton:
                gene_key = f"{parts[0]}\t{parts[1]}"
                singleton_genes.add(gene_key)
        
        nonsingleton_file = self.output_dir / "gene_catalog_nonsingletons.tsv"
        singleton_file = self.output_dir / "gene_catalog_singletons.tsv"
        
        with open(nonsingleton_file, 'w') as ns_f, open(singleton_file, 'w') as s_f:
            ns_f.write(header + "\n")
            s_f.write(header + "\n")
            
            for parts in gene_rows:
                gene_key = f"{parts[0]}\t{parts[1]}"
                line = "\t".join(parts) + "\n"
                if gene_key in singleton_genes:
                    s_f.write(line)
                else:
                    ns_f.write(line)
        
        logger.info(f"Split catalog: {len(singleton_genes)} singletons, {len(gene_rows) - len(singleton_genes)} non-singletons")
        return nonsingleton_file, singleton_file
    
    def _track_singletons_per_threshold(self, catalog_file):
        tracking_file = self.output_dir / "singleton_counts.tsv"
        
        with open(catalog_file, 'r') as f:
            header = f.readline().strip()
            lines = [line.strip() for line in f]
        
        header_parts = header.split('\t')
        if len(header_parts) < 3:
            return
        
        rep_columns = header_parts[2:]
        thresholds = rep_columns
        
        rep_counts_per_threshold = []
        for idx in range(len(thresholds)):
            rep_counts = {}
            for line in lines:
                parts = line.split('\t')
                if len(parts) < len(header_parts):
                    continue
                rep_val = parts[2 + idx]
                rep_counts[rep_val] = rep_counts.get(rep_val, 0) + 1
            rep_counts_per_threshold.append(rep_counts)
        
        singleton_counts = []
        for idx, rep_counts in enumerate(rep_counts_per_threshold):
            singleton_count = sum(1 for count in rep_counts.values() if count == 1)
            singleton_counts.append(singleton_count)
        
        file_exists = tracking_file.exists()
        with open(tracking_file, 'a') as f:
            if not file_exists:
                f.write("singleton_count\n")
            for count in singleton_counts:
                f.write(f"{count}\n")
        
        logger.info(f"Recorded singleton counts per threshold: {singleton_counts}")
    
    def build_catalog(self):
        logger.info("Starting gene catalog construction")
        
        input_files = self._get_input_files()
        logger.info(f"Found {len(input_files)} input file(s)")
        
        combined_fasta = self._combine_fasta_files(input_files)
        
        if self.identity_thresholds:
            logger.info(f"Iterative clustering at thresholds: {self.identity_thresholds}")
            clustering_result = self._iterative_cluster_chain(combined_fasta)
            final_catalog = self._create_gene_catalog(clustering_result, input_files, is_iterative=True)
        else:
            logger.info(f"Single-threshold clustering at identity {self.identity_threshold}")
            clustering_result = self._cluster_single_threshold(combined_fasta)
            final_catalog = self._create_gene_catalog(clustering_result, input_files, is_iterative=False)
        
        self._track_singletons_per_threshold(final_catalog)
        
        if self.split_singletons:
            self._split_singletons(final_catalog)
        
        logger.info(f"Gene catalog construction completed: {final_catalog}")
        return final_catalog

def main():
    parser = argparse.ArgumentParser(
        description="Build gene catalog by clustering protein sequences with MMseqs2"
    )
    
    parser.add_argument(
        '--sequence-dir',
        dest='sequence_dir',
        default=None,
        help='Directory containing FASTA files'
    )
    
    parser.add_argument(
        '--sequence-file',
        dest='sequence_file',
        default=None,
        help='Single FASTA file to cluster'
    )
    
    parser.add_argument(
        '-x', '--extension',
        default='faa',
        help='File extension to match when using --sequence-dir (default: faa)'
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
        '--identity-threshold',
        type=float,
        default=0.9,
        help='Identity threshold for MMseqs2 clustering (default: 0.9)'
    )
    
    parser.add_argument(
        '--identity-thresholds',
        type=str,
        default=None,
        help='Comma-separated identity thresholds for iterative clustering, e.g., 0.9,0.3'
    )
    
    parser.add_argument(
        '--coverage-threshold',
        type=float,
        default=0.8,
        help='Coverage threshold for MMseqs2 clustering (default: 0.8)'
    )
    
    parser.add_argument(
        '--tmpdir',
        default='./tmp/',
        help='Temporary directory for MMseqs2 (default: ./tmp/)'
    )
    
    parser.add_argument(
        '--split-singletons',
        action='store_true',
        help='Split singletons into separate file (gene_catalog_singletons.tsv)'
    )
    
    args = parser.parse_args()
    
    identity_thresholds = None
    if args.identity_thresholds:
        identity_thresholds = [float(x.strip()) for x in args.identity_thresholds.split(',')]
    
    builder = GeneCatalogBuilder(
        sequence_dir=args.sequence_dir,
        sequence_file=args.sequence_file,
        output_dir=args.output_dir,
        threads=args.threads,
        identity_threshold=args.identity_threshold,
        identity_thresholds=identity_thresholds,
        coverage_threshold=args.coverage_threshold,
        extension=args.extension,
        tmpdir=args.tmpdir,
        split_singletons=args.split_singletons
    )
    
    try:
        builder.build_catalog()
    except Exception as e:
        logger.error(f"Gene catalog construction failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
