#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import subprocess
import gzip
from pathlib import Path
from typing import Dict, Set, Tuple, Optional, TextIO
from concurrent.futures import ThreadPoolExecutor

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ReadFilter:
    def __init__(self, config_file: str, perq_dir: str, output_dir: str = 'filtered_reads',
                 min_kmers: int = 1, output_config: str = 'configs/filtered_reads_config',
                 max_workers: int = 1):
        self.config_file = config_file
        self.perq_dir = perq_dir
        self.output_dir = output_dir
        self.min_kmers = min_kmers
        self.output_config = output_config
        self.max_workers = max_workers
        self.reads_to_filter: Dict[str, Set[str]] = {}  # filename -> set of reads to filter
        
    def parse_config(self) -> Dict[str, Tuple[str, Optional[str]]]:
        """Parse the input config file."""
        config_data = {}
        with open(self.config_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                # Skip header line
                if parts[0].lower() == 'filename':
                    continue
                if len(parts) >= 2:
                    filename = parts[0]
                    pe1 = parts[1]
                    pe2 = parts[2] if len(parts) > 2 else None
                    config_data[filename] = (pe1, pe2)
        return config_data

    def prefilter_perq_file(self, perq_file: str) -> str:
        """Pre-filter perq file using grep to remove 'No matches found' lines and awk to filter by min_kmers."""
        logger.info(f"Pre-filtering {perq_file}")
        filtered_file = f"{perq_file}.filtered.tmp"
        try:
            # Use grep to remove 'No matches found' lines, then awk to filter by min_kmers
            with open(filtered_file, 'w') as outfile:
                subprocess.run(f"grep -v 'No matches found' {perq_file} | awk -F'\t' '$6 >= {self.min_kmers}'", shell=True, stdout=outfile, check=True)
            return filtered_file
        except subprocess.CalledProcessError as e:
            logger.error(f"Error pre-filtering {perq_file}: {e}")
            return perq_file

    def parse_perq_file(self, perq_file: str) -> Set[str]:
        """Parse a perq output file and return set of reads to filter."""
        reads_to_filter = set()
        logger.debug(f"Parsing perq file: {perq_file}")
        
        try:
            with open(perq_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:  # We need 6 columns
                        read_id = parts[0]
                        # Use column index 5 for kmer count
                        total_kmers = int(parts[5]) if parts[5].isdigit() else 0
                        if total_kmers >= self.min_kmers:
                            reads_to_filter.add(read_id)
            
            # Clean up temporary filtered file
            if perq_file.endswith('.filtered.tmp'):
                os.remove(perq_file)
                
        except Exception as e:
            logger.error(f"Error parsing {perq_file}: {e}")
            
        return reads_to_filter

    def load_perq_files(self):
        """Load all perq files and store reads to filter."""
        perq_files = [os.path.join(self.perq_dir, f) for f in os.listdir(self.perq_dir) if f.endswith('.perq')]
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_file = {executor.submit(self.prefilter_perq_file, perq_file): perq_file for perq_file in perq_files}
            for future in future_to_file:
                perq_file = future_to_file[future]
                try:
                    filtered_file = future.result()
                    base_name = os.path.splitext(os.path.basename(perq_file))[0]
                    logger.info(f"Processing perq file: {os.path.basename(perq_file)}")
                    self.reads_to_filter[base_name] = self.parse_perq_file(filtered_file)
                    logger.info(f"Found {len(self.reads_to_filter[base_name])} reads to filter in {os.path.basename(perq_file)}")
                except Exception as e:
                    logger.error(f"Error processing {perq_file}: {e}")

    def open_file(self, file_path: str, mode: str = 'r') -> TextIO:
        """Open a file, handling gzip compression if needed."""
        if file_path.endswith('.gz'):
            return gzip.open(file_path, mode + 't')
        return open(file_path, mode)

    def filter_fastq_file(self, fastq_file: str, reads_to_remove: Set[str], output_file: str):
        """Filter a FASTQ/FASTA file using seqkit grep."""
        try:
            # Add .gz extension to output file if not already present
            if not output_file.endswith('.gz'):
                output_file += '.gz'
            
            # Create a temporary file with read IDs to remove
            remove_file = f"{output_file}.remove"
            with open(remove_file, 'w') as f:
                for read_id in reads_to_remove:
                    f.write(f"{read_id}\n")
            
            # Use seqkit grep to filter the file
            cmd = f"seqkit grep -f {remove_file} -p -v {fastq_file} | gzip > {output_file}"
            logger.info(f"Filtering {fastq_file} to {output_file}")
            subprocess.run(cmd, shell=True, check=True)
            
            # Clean up temporary file
            os.remove(remove_file)
            
        except Exception as e:
            logging.error(f"Error filtering {fastq_file}: {str(e)}")
            raise

    def process_single_file(self, filename: str, pe1: str, pe2: Optional[str]) -> Tuple[str, str, Optional[str]]:
        """Process a single file completely in serial and return output paths."""
        logger.info(f"Processing file: {filename}")
        
        # Step 1: Process perq file for this sample
        perq_file = os.path.join(self.perq_dir, f"{filename}.perq")
        if not os.path.exists(perq_file):
            logger.warning(f"No perq file found for {filename}, skipping.")
            return filename, pe1, pe2
            
        # Pre-filter and parse perq file
        filtered_perq = self.prefilter_perq_file(perq_file)
        reads_to_filter = self.parse_perq_file(filtered_perq)
        logger.info(f"Found {len(reads_to_filter)} reads to filter in {filename}")
        
        # Step 2: Process fastq/fasta files for this sample
        pe1_basename = os.path.basename(pe1)
        pe1_output = os.path.join(self.output_dir, f"{filename}_{pe1_basename}")
        self.filter_fastq_file(pe1, reads_to_filter, pe1_output)
        
        # Process PE2 if it exists
        pe2_output = None
        if pe2:
            pe2_basename = os.path.basename(pe2)
            pe2_output = os.path.join(self.output_dir, f"{filename}_{pe2_basename}")
            self.filter_fastq_file(pe2, reads_to_filter, pe2_output)
        
        # Add .gz extension if not already present
        pe1_output = pe1_output if pe1_output.endswith('.gz') else pe1_output + '.gz'
        pe2_output = pe2_output if pe2_output and pe2_output.endswith('.gz') else (pe2_output + '.gz' if pe2_output else None)
        
        return filename, pe1_output, pe2_output

    def process_files(self):
        """Process all files according to config using parallel processing."""
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(os.path.dirname(self.output_config), exist_ok=True)
        
        # Get config data
        config_data = self.parse_config()
        
        # Process files in parallel using ThreadPoolExecutor
        results = []
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            future_to_file = {
                executor.submit(self.process_single_file, filename, pe1, pe2): filename
                for filename, (pe1, pe2) in config_data.items()
            }
            
            # Collect results as they complete
            for future in future_to_file:
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    logger.error(f"Error processing file: {e}")
        
        # Write results to config file
        with open(self.output_config, 'w') as out_config:
            # Write headers
            out_config.write("filename\tpe1\tpe2\n")
            for filename, pe1_output, pe2_output in results:
                if pe2_output:
                    out_config.write(f"{filename}\t{pe1_output}\t{pe2_output}\n")
                else:
                    out_config.write(f"{filename}\t{pe1_output}\n")

def main():
    parser = argparse.ArgumentParser(description='Filter reads based on perq output files.')
    parser.add_argument('--config', required=True, help='Input config file (filename, pe1, pe2)')
    parser.add_argument('--perq_dir', required=True, help='Directory containing perq output files')
    parser.add_argument('--output_dir', default='filtered_reads', help='Output directory for filtered reads')
    parser.add_argument('--min_kmers', type=int, default=1, help='Minimum number of kmers to keep a read')
    parser.add_argument('--output_config', default='configs/filtered_reads_config',
                       help='Output config file path')
    parser.add_argument('--max_workers', type=int, default=1,
                       help='Maximum number of parallel workers for processing files')
    args = parser.parse_args()

    # Initialize and run filter
    filter = ReadFilter(args.config, args.perq_dir, args.output_dir, args.min_kmers, 
                       args.output_config, args.max_workers)
    filter.process_files()
    logger.info("Read filtering completed successfully.")

if __name__ == '__main__':
    main() 