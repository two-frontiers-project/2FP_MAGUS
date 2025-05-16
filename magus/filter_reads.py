#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

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
                pe2 = parts[2] if len(parts) > 2 else None
                config_data[filename] = (pe1, pe2)
    return config_data

def filter_reads(perq_file: str, fastq_file: str, output_file: str, min_kmers: int, threads: int):
    """Filter reads using bash one-liner and seqkit."""
    tofilter = f"{output_file}.tofilter"
    cmd = f"grep -v 'No matches found' {perq_file} | awk -F'\t' '$6 > {min_kmers} {{print $1}}' > {tofilter}"
    subprocess.run(cmd, shell=True, check=True)
    cmd = f"seqkit grep -f {tofilter} -v {fastq_file} -j {threads} > {output_file}"
    subprocess.run(cmd, shell=True, check=True)
    cmd = f"gzip {output_file}"
    subprocess.run(cmd, shell=True, check=True)

def process_file(filename: str, pe1: str, pe2: str, perq_dir: str, output_dir: str, min_kmers: int, threads: int):
    """Process a single file pair."""
    perq_file = os.path.join(perq_dir, f"{filename}.perq")
    if not os.path.exists(perq_file):
        logger.warning(f"No perq file found for {filename}, skipping.")
        return
    
    pe1_basename = os.path.basename(pe1)
    pe1_output = os.path.join(output_dir, f"{filename}_{pe1_basename}")
    filter_reads(perq_file, pe1, pe1_output, min_kmers, threads)
    
    if pe2:
        pe2_basename = os.path.basename(pe2)
        pe2_output = os.path.join(output_dir, f"{filename}_{pe2_basename}")
        filter_reads(perq_file, pe2, pe2_output, min_kmers, threads)

def main():
    parser = argparse.ArgumentParser(description='Filter reads based on perq output files.')
    parser.add_argument('--config', required=True, help='Input config file (filename, pe1, pe2)')
    parser.add_argument('--perq_dir', required=True, help='Directory containing perq output files')
    parser.add_argument('--output_dir', default='filtered_reads', help='Output directory for filtered reads')
    parser.add_argument('--min_kmers', type=int, default=10, help='Minimum number of kmers to keep a read')
    parser.add_argument('--max_workers', type=int, default=1, help='Maximum number of parallel workers')
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