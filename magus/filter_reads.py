#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import csv

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_config(config_file: str) -> dict:
    """Parse the input config file."""
    config_data = {}
    with open(config_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[0].lower() == 'filename':
                continue
            if len(parts) >= 2:
                filename = parts[0]
                pe1 = parts[1]
                pe2 = parts[2] if len(parts) > 2 and parts[2].strip() else None
                config_data[filename] = (pe1, pe2)
    return config_data

def filter_reads(perq_file: str, fastq_file: str, output_file: str, min_kmers: int, threads: int):
    """Filter reads using bash one-liner and seqkit."""
    tofilter = f"{output_file}.tofilter"
    # Remove .gz extension for intermediate file
    uncompressed_output = output_file.replace('.gz', '')
    
    cmd = f"grep -v 'No matches found' {perq_file} | awk -F'\t' '$6 > {min_kmers} {{print $1}}' > {tofilter}"
    subprocess.run(cmd, shell=True, check=True)
    cmd = f"seqkit grep -f {tofilter} -v {fastq_file} -j {threads} -w 0 > {uncompressed_output}"
    subprocess.run(cmd, shell=True, check=True)
    cmd = f"gzip {uncompressed_output}"
    subprocess.run(cmd, shell=True, check=True)

def process_file(filename: str, pe1: str, pe2: str, perq_dir: str, output_dir: str, min_kmers: int, threads: int):
    """Process a single file or file pair."""
    perq_file = os.path.join(perq_dir, f"{filename}.perq")
    if not os.path.exists(perq_file):
        logger.warning(f"No perq file found for {filename}, skipping.")
        return
    
    # Determine output extension based on input file
    pe1_ext = '.fa.gz' if pe1.endswith(('.fa', '.fasta')) else '.fq.gz'
    
    if pe2:
        # Paired-end mode - apply same filtering to both files
        logger.info(f"Processing paired-end files for {filename}")
        pe2_ext = '.fa.gz' if pe2.endswith(('.fa', '.fasta')) else '.fq.gz'
        
        # Use R1/R2 naming convention for output files
        r1_output = os.path.join(output_dir, f"{filename}.R1{pe1_ext}")
        r2_output = os.path.join(output_dir, f"{filename}.R2{pe2_ext}")
        
        # Apply same filtering logic to both files
        filter_reads(perq_file, pe1, r1_output, min_kmers, threads)
        filter_reads(perq_file, pe2, r2_output, min_kmers, threads)
        logger.info(f"Filtered paired-end reads for {filename}: {r1_output}, {r2_output}")
    else:
        # Single-end mode
        logger.info(f"Processing single-end file for {filename}")
        pe1_output = os.path.join(output_dir, f"{filename}_filtered{pe1_ext}")
        
        filter_reads(perq_file, pe1, pe1_output, min_kmers, threads)
        logger.info(f"Filtered single-end reads for {filename}: {pe1_output}")

def main():
    parser = argparse.ArgumentParser(description='Filter reads based on perq output files.')
    parser.add_argument('--config', required=True, help='Input config file (filename, pe1, pe2)')
    parser.add_argument('--perq-dir', '--perq_dir', dest='perq_dir', required=True, help='Directory containing perq output files')
    parser.add_argument('--output-dir', '--output_dir', dest='output_dir', default='filtered_reads', help='Output directory for filtered reads')
    parser.add_argument('--min_kmers', type=int, default=10, help='Minimum number of kmers to keep a read')
    parser.add_argument('--max-workers', '--max_workers', dest='max_workers', type=int, default=1, help='Maximum number of parallel workers')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads for seqkit')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    config_data = parse_config(args.config)

    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = []
        for filename, (pe1, pe2) in config_data.items():
            future = executor.submit(process_file, filename, pe1, pe2, args.perq_dir, args.output_dir, args.min_kmers, args.threads)
            futures.append(future)
        
        # Wait for all tasks to complete
        for future in futures:
            future.result()

if __name__ == '__main__':
    main() 