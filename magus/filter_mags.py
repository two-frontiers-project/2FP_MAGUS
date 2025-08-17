#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import subprocess
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MAGFilter:
    def __init__(self, output_dir, perq_dir, mag_dir, kmer_threshold=10):
        """
        Initialize the MAG filter.
        
        Args:
            output_dir: Output directory for filtered MAGs
            perq_dir: Directory containing perq files (raw_alignments)
            mag_dir: Directory containing MAG files to filter
            kmer_threshold: K-mer count threshold for filtering (default: 10)
        """
        self.output_dir = Path(output_dir)
        self.perq_dir = Path(perq_dir)
        self.mag_dir = Path(mag_dir)
        self.kmer_threshold = kmer_threshold
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Validate input directories
        self._validate_directories()
    
    def _validate_directories(self):
        """Validate that input directories exist."""
        if not self.perq_dir.exists():
            raise FileNotFoundError(f"Perq directory not found: {self.perq_dir}")
        if not self.mag_dir.exists():
            raise FileNotFoundError(f"MAG directory not found: {self.mag_dir}")
    
    def _get_perq_files(self):
        """Get all perq files from the perq directory."""
        perq_files = list(self.perq_dir.glob("*.perq"))
        if not perq_files:
            raise FileNotFoundError(f"No .perq files found in {self.perq_dir}")
        logger.info(f"Found {len(perq_files)} perq files")
        return perq_files
    
    def _create_all_perq_data(self, perq_files):
        """Combine all perq files and add sample ID column."""
        all_perq_data_path = self.output_dir / "all_perq_data"
        
        logger.info("Creating combined perq data file")
        
        with open(all_perq_data_path, 'w') as outfile:
            for perq_file in perq_files:
                sample_id = perq_file.stem  # filename without .perq extension
                
                try:
                    with open(perq_file, 'r') as infile:
                        for line in infile:
                            line = line.strip()
                            if line and "No matches found" not in line:
                                # Split the line and extract columns 1, 2, and 6 (0-indexed: 0, 1, 5)
                                parts = line.split('\t')
                                if len(parts) >= 6:
                                    # Format: sample_id, contig, reference, kmer_count
                                    contig = parts[0]  # Column 1: contig ID
                                    reference = parts[1]  # Column 2: reference name
                                    kmer_count = parts[5]  # Column 6: kmer count
                                    outfile.write(f"{sample_id}\t{contig}\t{reference}\t{kmer_count}\n")
                except Exception as e:
                    logger.warning(f"Error processing {perq_file}: {e}")
                    continue
        
        logger.info(f"Created combined perq data file: {all_perq_data_path}")
        return all_perq_data_path
    
    def _parse_perq_data(self, all_perq_data_path):
        """Parse the combined perq data to identify contigs to remove."""
        to_remove = defaultdict(set)
        
        logger.info("Parsing perq data to identify contigs for removal")
        
        try:
            with open(all_perq_data_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 4:
                        logger.warning(f"Line {line_num}: insufficient columns, skipping")
                        continue
                    
                    sample_id, contig, query, kmer_count = parts
                    
                    try:
                        kmer_count_int = int(kmer_count)
                        if kmer_count_int > self.kmer_threshold:
                            to_remove[sample_id].add(contig)
                    except ValueError:
                        logger.warning(f"Line {line_num}: invalid kmer count '{kmer_count}', skipping")
                        continue
            
            logger.info(f"Identified contigs to remove for {len(to_remove)} samples")
            for sample_id, contigs in to_remove.items():
                logger.info(f"  {sample_id}: {len(contigs)} contigs to remove")
                
        except Exception as e:
            logger.error(f"Error parsing perq data: {e}")
            raise
        
        return to_remove
    
    def _filter_mags(self, to_remove):
        """Filter each MAG file by removing identified contigs."""
        logger.info("Filtering MAG files")
        
        filtered_count = 0
        total_contigs_removed = 0
        
        for sample_id, contigs_to_remove in to_remove.items():
            # Look for the corresponding MAG file - it should have the exact same name
            mag_file = self.mag_dir / sample_id
            
            if not mag_file.exists():
                logger.warning(f"MAG file not found for sample {sample_id}")
                continue
            
            try:
                # Read and filter the MAG file
                records = list(SeqIO.parse(mag_file, "fasta"))
                kept = [r for r in records if r.id not in contigs_to_remove]
                removed = len(records) - len(kept)
                
                # Write filtered output
                output_file = self.output_dir / f"filtered_{mag_file.name}"
                SeqIO.write(kept, output_file, "fasta")
                
                filtered_count += 1
                total_contigs_removed += removed
                
                logger.info(f"{sample_id}: removed {removed} contigs, kept {len(kept)} contigs")
                
            except Exception as e:
                logger.error(f"Error filtering {sample_id}: {e}")
                continue
        
        logger.info(f"Filtered {filtered_count} MAG files, removed {total_contigs_removed} total contigs")
    
    def filter_mags(self):
        """Main method to filter MAGs based on perq data."""
        logger.info("Starting MAG filtering process")
        
        try:
            # Get perq files
            perq_files = self._get_perq_files()
            
            # Create combined perq data
            all_perq_data_path = self._create_all_perq_data(perq_files)
            
            # Parse perq data to identify contigs to remove
            to_remove = self._parse_perq_data(all_perq_data_path)
            
            # Filter MAGs
            self._filter_mags(to_remove)
            
            logger.info("MAG filtering completed successfully")
            
        except Exception as e:
            logger.error(f"MAG filtering failed: {e}")
            raise

def main():
    parser = argparse.ArgumentParser(
        description="Filter MAGs based on alignment data from taxonomy output"
    )
    
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for filtered MAGs'
    )
    
    parser.add_argument(
        '--perq-dir',
        required=True,
        help='Directory containing perq files (raw_alignments)'
    )
    
    parser.add_argument(
        '--mag-dir',
        required=True,
        help='Directory containing MAG files to filter'
    )
    
    parser.add_argument(
        '--kmer-threshold',
        type=int,
        default=10,
        help='K-mer count threshold for filtering (default: 10)'
    )
    
    args = parser.parse_args()
    
    # Create and run the MAG filter
    mag_filter = MAGFilter(
        output_dir=args.output_dir,
        perq_dir=args.perq_dir,
        mag_dir=args.mag_dir,
        kmer_threshold=args.kmer_threshold
    )
    
    try:
        mag_filter.filter_mags()
    except Exception as e:
        logger.error(f"MAG filtering failed: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
