#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import csv
import glob
import re
from typing import Optional
# pandas is not required here

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ORFCaller:
    def __init__(self, output_dir, extension, args):
        self.output_dir = output_dir
        # Ensure extension starts with a dot for consistent usage
        self.extension = f".{extension}" if not extension.startswith('.') else extension
        self.args = args

    def call_bacterial_orfs(self, genome_file, sample_id=None):
        """Call ORFs in bacterial genomes using prodigal."""
        annot_dir = os.path.join(self.output_dir, 'bacteria', 'annot')
        manicure_dir = os.path.join(self.output_dir, 'bacteria', 'manicure')
        os.makedirs(annot_dir, exist_ok=True)
        os.makedirs(manicure_dir, exist_ok=True)

        FN = sample_id if sample_id else os.path.basename(genome_file).replace(self.extension, '')
        
        # Check if the output files already exist and have nonzero length
        faa_file = os.path.join(annot_dir, f"{FN}.faa")
        ffn_file = os.path.join(annot_dir, f"{FN}.ffn")
        gff_file = os.path.join(annot_dir, f"{FN}.gff")
        
        if (os.path.exists(faa_file) and os.path.getsize(faa_file) > 0 and 
            os.path.exists(ffn_file) and os.path.getsize(ffn_file) > 0 and 
            os.path.exists(gff_file) and os.path.getsize(gff_file) > 0 and 
            not self.args.force):
            logger.info(f"Skipping ORF calling for {genome_file} as output already exists.")
        else:
            log_file = os.path.join(self.output_dir, 'bacteria', f"{FN}_prodigal.log")
            cmd = ['prodigal', '-i', genome_file, '-d', ffn_file, '-a', faa_file, '-o', gff_file, '-f', 'gff']
            logger.info(f"Calling bacterial ORFs for {genome_file} using prodigal")
            with open(log_file, 'w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)

            # Manicure the output files
            manicure_file = os.path.join(manicure_dir, f"{FN}.faa")
            with open(os.path.join(annot_dir, f"{FN}.faa"), 'r') as infile, open(manicure_file, 'w') as outfile:
                for line in infile:
                    if line.startswith('>'):
                        outfile.write(line.replace('>',f'>{FN}-----').replace(' # ', '-----', 1).replace('*', 'X').replace(' # ', '+').replace(' # ', '-'))
                    else:
                        outfile.write(line)

            manicure_ffn = os.path.join(manicure_dir, f"{FN}.ffn")
            with open(os.path.join(annot_dir, f"{FN}.ffn"), 'r') as infile, open(manicure_ffn, 'w') as outfile:
                for line in infile:
                    if line.startswith('>'):
                        outfile.write(line.replace('>',f'>{FN}-----').replace('+', '-----+').replace('*', 'X').replace(' # ', '+').replace(' # ', '-'))
                    else:
                        outfile.write(line)

            # Keep the annot files for the comprehensive summary
            # (Files are kept in annot directory for summary generation)

        # HMM annotation (put results in annot directory for summary)
        manicure_file = os.path.join(manicure_dir, f"{FN}.faa")
        self.run_hmm_annotation(manicure_file, annot_dir, FN, self.args.hmmfile)

        # Remove MetaEuk temporary directory named after the sample within annot_dir
        sample_tmp_dir = os.path.join(annot_dir, FN)
        if os.path.isdir(sample_tmp_dir):
            import shutil
            shutil.rmtree(sample_tmp_dir, ignore_errors=True)
        return manicure_file

    def call_viral_orfs(self, genome_file, sample_id=None):
        """Call ORFs in viral genomes using prodigal-gv."""
        annot_dir = os.path.join(self.output_dir, 'viruses', 'annot')
        manicure_dir = os.path.join(self.output_dir, 'viruses', 'manicure')
        os.makedirs(annot_dir, exist_ok=True)
        os.makedirs(manicure_dir, exist_ok=True)

        FN = sample_id if sample_id else os.path.basename(genome_file).replace(self.extension, '')
        
        # Check if the output files already exist and have nonzero length
        faa_file = os.path.join(annot_dir, f"{FN}.faa")
        ffn_file = os.path.join(annot_dir, f"{FN}.ffn")
        gff_file = os.path.join(annot_dir, f"{FN}.gff")
        
        if (os.path.exists(faa_file) and os.path.getsize(faa_file) > 0 and 
            os.path.exists(ffn_file) and os.path.getsize(ffn_file) > 0 and 
            os.path.exists(gff_file) and os.path.getsize(gff_file) > 0 and 
            not self.args.force):
            logger.info(f"Skipping ORF calling for {genome_file} as output already exists.")
        else:
            log_file = os.path.join(self.output_dir, 'viruses', f"{FN}_prodigal.log")
            cmd = ['prodigal-gv', '-p', '-q', '-i', genome_file, '-d', ffn_file, '-a', faa_file, '-o', gff_file, '-f', 'gff']
            logger.info(f"Calling viral ORFs for {genome_file} using prodigal-gv")
            with open(log_file, 'w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)
        
        # Manicure the output files (same as bacteria since it's also Prodigal)
        manicure_file = os.path.join(manicure_dir, f"{FN}.faa")
        with open(faa_file, 'r') as infile, open(manicure_file, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    outfile.write(line.replace('>',f'>{FN}-----').replace(' # ', '-----', 1).replace('*', 'X').replace(' # ', '+').replace(' # ', '-'))
                else:
                    outfile.write(line)

        manicure_ffn = os.path.join(manicure_dir, f"{FN}.ffn")
        with open(ffn_file, 'r') as infile, open(manicure_ffn, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    outfile.write(line.replace('>',f'>{FN}-----').replace('+', '-----+').replace('*', 'X').replace(' # ', '+').replace(' # ', '-'))
                else:
                    outfile.write(line)

        # HMM annotation (put results in annot directory for summary)
        if self.args.hmmfile:
            self.run_hmm_annotation(faa_file, annot_dir, FN, self.args.hmmfile)
        
        return faa_file, ffn_file, gff_file

    def call_eukaryotic_orfs(self, genome_file, sample_id=None):
        """Call ORFs in eukaryotic genomes using MetaEuk."""
        annot_dir = os.path.join(self.output_dir, 'eukaryotes', 'annot')
        manicure_dir = os.path.join(self.output_dir, 'eukaryotes', 'manicure')
        os.makedirs(annot_dir, exist_ok=True)
        os.makedirs(manicure_dir, exist_ok=True)

        FN = sample_id if sample_id else os.path.basename(genome_file).replace(self.extension, '')
        
        # Check if the output files already exist and have nonzero length
        faa_file = os.path.join(annot_dir, f"{FN}.fas")
        ffn_file = os.path.join(annot_dir, f"{FN}.codon.fas")
        
        if (os.path.exists(faa_file) and os.path.getsize(faa_file) > 0 and 
            os.path.exists(ffn_file) and os.path.getsize(ffn_file) > 0 and 
            not self.args.force):
            logger.info(f"Skipping ORF calling for {genome_file} as output already exists.")
        else:
            log_file = os.path.join(self.output_dir, 'eukaryotes', f"{FN}_metaeuk.log")
            cmd = ['metaeuk', 'easy-predict', genome_file, self.args.eukdb, os.path.join(annot_dir, f"{FN}"), os.path.join(annot_dir, f"{FN}"), '--threads', str(self.args.threads)]
            logger.info(f"Calling eukaryotic ORFs for {genome_file} using MetaEuk with {self.args.threads} threads")
            with open(log_file, 'w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)
        
        # Symlink MetaEuk output to manicure directory
        manicure_faa = os.path.join(manicure_dir, f"{FN}.faa")
        manicure_ffn = os.path.join(manicure_dir, f"{FN}.ffn")
        
        # Create symlinks (or copy if symlink fails)
        try:
            if os.path.exists(manicure_faa):
                os.remove(manicure_faa)
            if os.path.exists(manicure_ffn):
                os.remove(manicure_ffn)
            os.symlink(faa_file, manicure_faa)
            os.symlink(ffn_file, manicure_ffn)
        except OSError:
            # Fallback to copying if symlink fails
            import shutil
            shutil.copy2(faa_file, manicure_faa)
            shutil.copy2(ffn_file, manicure_ffn)
        
        # HMM annotation (put results in annot directory for summary)
        if self.args.hmmfile:
            hmm_out = self.run_hmm_annotation(faa_file, annot_dir, FN, self.args.hmmfile)
            
            # Parse HMM output to clean CSV and merge with FASTA data
            if hmm_out and os.path.exists(hmm_out):
                # Parse HMM tblout to CSV with e-value filtering
                hmm_csv = os.path.join(annot_dir, f"{FN}.hmm_clean.csv")
                self.parse_hmm_tblout_to_csv(hmm_out, hmm_csv, self.args.annotation_domain_evalue)
                
                # Read FASTA headers and merge with HMM data
                fasta_headers = self.read_fasta_headers(faa_file)
                if fasta_headers:
                    summary_file = os.path.join(annot_dir, f"{FN}_summary.tsv")
                    self.merge_eukaryotic_hmm_data(FN, fasta_headers, hmm_csv, summary_file)
        
        return faa_file, ffn_file, None  # MetaEuk doesn't produce GFF

    def call_metagenome_orfs(self, genome_file, sample_id=None):
        """Call ORFs in metagenome mode using prodigal."""
        annot_dir = os.path.join(self.output_dir, 'metagenomes', 'annot')
        manicure_dir = os.path.join(self.output_dir, 'metagenomes', 'manicure')
        os.makedirs(annot_dir, exist_ok=True)
        os.makedirs(manicure_dir, exist_ok=True)

        FN = sample_id if sample_id else os.path.basename(genome_file).replace(self.extension, '')
        
        # Check if the output files already exist and have nonzero length
        faa_file = os.path.join(annot_dir, f"{FN}.faa")
        ffn_file = os.path.join(annot_dir, f"{FN}.ffn")
        gff_file = os.path.join(annot_dir, f"{FN}.gff")
        
        if (os.path.exists(faa_file) and os.path.getsize(faa_file) > 0 and 
            os.path.exists(ffn_file) and os.path.getsize(ffn_file) > 0 and 
            os.path.exists(gff_file) and os.path.getsize(gff_file) > 0 and 
            not self.args.force):
            logger.info(f"Skipping ORF calling for {genome_file} as output already exists.")
        else:
            log_file = os.path.join(self.output_dir, 'metagenomes', f"{FN}_prodigal.log")
            cmd = ['prodigal', '-p', 'meta', '-q', '-i', genome_file, '-d', ffn_file, '-a', faa_file, '-o', gff_file, '-f', 'gbk']
            logger.info(f"Calling metagenome ORFs for {genome_file} using prodigal (metagenome mode)")
            with open(log_file, 'w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)
        
        # Manicure the output files (same as bacteria since it's also Prodigal)
        manicure_file = os.path.join(manicure_dir, f"{FN}.faa")
        with open(faa_file, 'r') as infile, open(manicure_file, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    outfile.write(line.replace('>',f'>{FN}-----').replace(' # ', '-----', 1).replace('*', 'X').replace(' # ', '+').replace(' # ', '-'))
                else:
                    outfile.write(line)

        manicure_ffn = os.path.join(manicure_dir, f"{FN}.ffn")
        with open(ffn_file, 'r') as infile, open(manicure_ffn, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    outfile.write(line.replace('>',f'>{FN}-----').replace('+', '-----+').replace('*', 'X').replace(' # ', '+').replace(' # ', '-'))
                else:
                    outfile.write(line)

        # HMM annotation (put results in annot directory for summary)
        if self.args.hmmfile:
            self.run_hmm_annotation(faa_file, annot_dir, FN, self.args.hmmfile)
        
        return faa_file, ffn_file, gff_file

    def run_hmm_annotation(self, faa_file, output_dir, sample_id, hmmfile, hmm_suffix=None):
        """Run HMM annotation on a protein file.

        If hmm_suffix is provided, write tblout to
        "{sample_id}.hmm.{suffix}.tsv" so multiple annotation sets can coexist.
        """
        if not hmmfile:
            logger.warning(f"No HMM file provided for {sample_id}. Skipping annotation.")
            return None
            
        if not os.path.exists(faa_file):
            logger.warning(f"Protein file {faa_file} does not exist. Skipping annotation.")
            return None
            
        if hmm_suffix:
            hmm_out = os.path.join(output_dir, f"{sample_id}.hmm.{hmm_suffix}.tsv")
        else:
            hmm_out = os.path.join(output_dir, f"{sample_id}.hmm.tsv")
        cmd = [
            'hmmsearch',
            '--tblout', hmm_out,
            '--noali',
            '--notextw',
            '--cpu', str(self.args.threads),
            hmmfile,
            faa_file
        ]
        logger.info(f"Running HMM annotation for {sample_id}")
        with open(os.path.join(output_dir, f"{sample_id}_hmmsearch.log"), 'w') as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=log)
        
        return hmm_out

    def parse_hmm_tblout_to_csv(self, hmm_file, output_csv, evalue_cutoff=0.01):
        """Parse HMM tblout output to clean CSV with e-value filtering.
        
        This function works for all domains (bacteria, viruses, eukaryotes, metagenomes).
        It filters on both full-sequence and domain e-values immediately upon loading.
        For each ORF, it keeps only the hit with the lowest full-sequence e-value.
        """
        if not os.path.exists(hmm_file):
            logger.warning(f"HMM file {hmm_file} does not exist. Skipping parsing.")
            return None
            
        # HMM tblout columns (removing target_accession which is just '-')
        columns = [
            'target_name', 'query_name', 'query_accession',
            'full_evalue', 'full_score', 'full_bias', 'dom_evalue', 'dom_score', 
            'dom_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description'
        ]
        
        # Dictionary to store best hit per ORF (target_name)
        best_hits = {}
        total_lines = 0
        filtered_lines = 0
        
        with open(hmm_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                    
                total_lines += 1
                parts = line.split()
                
                if len(parts) < 19:  # HMM tblout has 19 columns
                    continue
                
                # Extract evalue for filtering
                try:
                    full_evalue = float(parts[4])
                    dom_evalue = float(parts[7])
                except (ValueError, IndexError):
                    continue
                
                # Apply evalue filter immediately
                if full_evalue > evalue_cutoff and dom_evalue > evalue_cutoff:
                    continue
                    
                filtered_lines += 1
                
                # Create row data (skip target_accession which is parts[1] and always '-')
                row = {}
                # Map columns to parts, skipping target_accession (parts[1])
                row['target_name'] = parts[0]      # parts[0]
                row['query_name'] = parts[2]       # parts[2] (skip parts[1] which is target_accession)
                row['query_accession'] = parts[3]  # parts[3]
                row['full_evalue'] = parts[4]      # parts[4]
                row['full_score'] = parts[5]       # parts[5]
                row['full_bias'] = parts[6]        # parts[6]
                row['dom_evalue'] = parts[7]       # parts[7]
                row['dom_score'] = parts[8]        # parts[8]
                row['dom_bias'] = parts[9]         # parts[9]
                row['exp'] = parts[10]             # parts[10]
                row['reg'] = parts[11]             # parts[11]
                row['clu'] = parts[12]             # parts[12]
                row['ov'] = parts[13]              # parts[13]
                row['env'] = parts[14]             # parts[14]
                row['dom'] = parts[15]             # parts[15]
                row['rep'] = parts[16]             # parts[16]
                row['inc'] = parts[17]             # parts[17]
                row['description'] = ' '.join(parts[18:])  # parts[18:] for description
                
                # Keep only the best hit (lowest full_evalue) for each ORF
                target_name = row['target_name']
                if target_name not in best_hits or float(row['full_evalue']) < float(best_hits[target_name]['full_evalue']):
                    best_hits[target_name] = row
        
        # Convert best hits to list
        rows = list(best_hits.values())
        
        # Write to CSV
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
        
        return rows

    def read_fasta_headers(self, fas_file):
        """Read FASTA headers from .fas file."""
        if not os.path.exists(fas_file):
            logger.warning(f"FASTA file {fas_file} does not exist.")
            return set()
            
        headers = set()
        with open(fas_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line.strip()[1:]  # Remove '>'
                    headers.add(header)
        
        return headers

    def merge_eukaryotic_hmm_data(self, sample_id, fasta_headers, hmm_csv_file, output_file):
        """Merge eukaryotic FASTA headers with HMM data and write to output file.
        
        This creates the proper eukaryotic summary format:
        sample_id \t target_name \t hmm_annotation_data
        """
        if not os.path.exists(hmm_csv_file):
            # Create summary with just FASTA headers and no HMM data
            with open(output_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                header = ['sample_id', 'target_name', 'query_name', 'query_accession',
                         'full_evalue', 'full_score', 'full_bias', 'dom_evalue', 'dom_score', 
                         'dom_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description']
                writer.writerow(header)
                
                for target_name in fasta_headers:
                    row = [sample_id, target_name] + [''] * (len(header) - 2)
                    writer.writerow(row)
            
            logger.info(f"Created summary for {sample_id}: {len(fasta_headers)} ORFs (no HMM data)")
            return len(fasta_headers), 0, len(fasta_headers)
        
        # Read HMM data from CSV
        hmm_data = {}
        with open(hmm_csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                target_name = row['target_name']
                hmm_data[target_name] = row
        
        # HMM annotation columns (excluding target_name since it's in the merge)
        hmm_cols = [
            'query_name', 'query_accession', 'full_evalue', 'full_score', 'full_bias',
            'dom_evalue', 'dom_score', 'dom_bias', 'exp', 'reg', 'clu', 'ov', 'env',
            'dom', 'rep', 'inc', 'description'
        ]
        
        # Write merged output
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            header = ['sample_id', 'target_name'] + hmm_cols
            writer.writerow(header)
            
            # Write data rows
            rows_written = 0
            rows_with_hmm = 0
            rows_without_hmm = 0
            
            for target_name in fasta_headers:
                if target_name in hmm_data:
                    # This target has HMM hits
                    hmm_info = hmm_data[target_name]
                    row = [sample_id, target_name]
                    for col in hmm_cols:
                        row.append(hmm_info.get(col, ''))
                    writer.writerow(row)
                    rows_with_hmm += 1
                else:
                    # This target has no HMM hits
                    row = [sample_id, target_name] + [''] * len(hmm_cols)
                    writer.writerow(row)
                    rows_without_hmm += 1
                rows_written += 1
        
        logger.info(f"Created summary for {sample_id}: {rows_written} ORFs ({rows_with_hmm} with HMM hits, {rows_without_hmm} without)")
        
        return rows_written, rows_with_hmm, rows_without_hmm

def find_genome_files(mag_dir, extension, wildcard):
    """Find genome files using directory and wildcard pattern, similar to dereplicate_genomes.py"""
    input_paths = [Path(p).resolve() for p in glob.glob(mag_dir, recursive=True)]
    all_matches = []
    
    for base_path in input_paths:
        if base_path.is_dir():
            suffix = f".{extension}" if not extension.startswith('.') else extension
            for path in base_path.rglob(f"*{suffix}"):
                if wildcard in str(path):
                    all_matches.append(path.resolve())
    
    if not all_matches:
        raise RuntimeError("No matching genome files found.")
    
    return all_matches

def process_genomes(orf_caller, genomes_data, max_workers=1):
    """Process a list of genomes data (tuples of sample_id, genome_path, domain)"""
    
    def process_single_genome(genome_data):
        sample_id, genome_path, domain = genome_data
        domain = domain.lower()
        try:
            if domain == 'bacterial':
                orf_caller.call_bacterial_orfs(genome_path, sample_id=sample_id)
            elif domain == 'viral':
                orf_caller.call_viral_orfs(genome_path, sample_id=sample_id)
            elif domain == 'eukaryotic':
                orf_caller.call_eukaryotic_orfs(genome_path, sample_id=sample_id)
            elif domain == 'metagenomic':
                orf_caller.call_metagenome_orfs(genome_path, sample_id=sample_id)
            else:
                logger.error(f"Unknown domain '{domain}' for sample {sample_id}. Skipping.")
                return False
            return True
        except Exception as e:
            logger.error(f"Error processing {sample_id}: {e}")
            return False
    
    if max_workers > 1:
        logger.info(f"Processing {len(genomes_data)} genomes with {max_workers} parallel workers")
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(process_single_genome, genomes_data))
        successful = sum(results)
        logger.info(f"Successfully processed {successful}/{len(genomes_data)} genomes")
    else:
        logger.info(f"Processing {len(genomes_data)} genomes sequentially")
        for genome_data in genomes_data:
            process_single_genome(genome_data)

def create_comprehensive_summary(output_dir, hmmfile, suffix=None,
                                 hmm_fullseq_evalue_cutoff: Optional[float] = None,
                                 hmm_domain_evalue_cutoff: Optional[float] = None):
    """Create ONE comprehensive summary file per domain with ALL information.

    If suffix is provided, domain summary is named with _{suffix} and HMM inputs
    are read from .hmm.{suffix}.tsv files.
    """
    for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
        annot_dir = os.path.join(output_dir, subdir, 'annot')
        if not os.path.exists(annot_dir):
            continue
            
        summary_file = os.path.join(
            output_dir,
            f"{subdir}_orf_summary_{suffix}.tsv" if suffix else f"{subdir}_orf_summary.tsv"
        )
        logger.info(f"Creating comprehensive summary for {subdir}: {summary_file}")
        
        with open(summary_file, 'w') as summary:
            if subdir == 'eukaryotes':
                logger.info("Processing eukaryotic samples")
                
                # Create temporary ORFCaller instance to use the methods
                class TempArgs:
                    def __init__(self):
                        self.annotation_domain_evalue = hmm_domain_evalue_cutoff or 0.01
                        self.threads = 1
                
                temp_args = TempArgs()
                temp_caller = ORFCaller(output_dir, "fas", temp_args)
                
                # Write header once
                header = ['sample_id', 'target_name', 'query_name', 'query_accession',
                         'full_evalue', 'full_score', 'full_bias', 'dom_evalue', 'dom_score', 
                         'dom_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description']
                summary.write('\t'.join(header) + '\n')
                
                # Process each sample: read FASTA, clean HMM, merge, write to final summary
                for file in os.listdir(annot_dir):
                    if file.endswith('.fas') and not file.endswith('.codon.fas'):
                        sample_id = file.replace('.fas', '')
                        fas_file = os.path.join(annot_dir, file)
                        
                        # Read FASTA headers for this sample
                        fasta_headers = temp_caller.read_fasta_headers(fas_file)
                        if not fasta_headers:
                            continue
                        
                        # Check if HMM file exists and clean it
                        hmm_file = os.path.join(annot_dir, f"{sample_id}.hmm.tsv")
                        sample_hmm_data = {}
                        
                        if os.path.exists(hmm_file):
                            hmm_csv = os.path.join(annot_dir, f"{sample_id}.hmm_clean.csv")
                            rows = temp_caller.parse_hmm_tblout_to_csv(hmm_file, hmm_csv, hmm_domain_evalue_cutoff)
                            if rows:
                                for row in rows:
                                    target_name = row['target_name']
                                    sample_hmm_data[target_name] = row
                        
                        # Write all rows for this sample to the final summary
                        for target_name in fasta_headers:
                            # ALWAYS write a row for every FASTA header (left join)
                            row_data = [sample_id, target_name]
                            
                            # Add HMM data if it exists, otherwise empty strings
                            if target_name in sample_hmm_data:
                                hmm_row = sample_hmm_data[target_name]
                                row_data.extend([
                                    hmm_row.get('query_name', ''),
                                    hmm_row.get('query_accession', ''),
                                    hmm_row.get('full_evalue', ''),
                                    hmm_row.get('full_score', ''),
                                    hmm_row.get('full_bias', ''),
                                    hmm_row.get('dom_evalue', ''),
                                    hmm_row.get('dom_score', ''),
                                    hmm_row.get('dom_bias', ''),
                                    hmm_row.get('exp', ''),
                                    hmm_row.get('reg', ''),
                                    hmm_row.get('clu', ''),
                                    hmm_row.get('ov', ''),
                                    hmm_row.get('env', ''),
                                    hmm_row.get('dom', ''),
                                    hmm_row.get('rep', ''),
                                    hmm_row.get('inc', ''),
                                    hmm_row.get('description', '')
                                ])
                            else:
                                # No HMM data for this target - add empty strings
                                row_data.extend([''] * 17)
                            
                            # WRITE THE ROW - this should happen for EVERY FASTA header
                            summary.write('\t'.join(str(field) for field in row_data) + '\n')
            else:
                # Bacteria/Viruses/Metagenomes: simple GFF + HMM merge
                logger.info(f"Processing {subdir} samples")
                
                import pandas as pd
                
                # Create temporary ORFCaller instance to use the methods
                class TempArgs:
                    def __init__(self):
                        self.annotation_domain_evalue = hmm_domain_evalue_cutoff or 0.01
                        self.threads = 1
                
                temp_args = TempArgs()
                
                # Collect all data first, then write at once
                all_merged_data = []
                
                manicure_dir = os.path.join(output_dir, subdir, 'manicure')
                
                for file in os.listdir(manicure_dir):
                    if file.endswith('.faa'):
                        sample_id = file.replace('.faa', '')
                        
                        # 1. Load FAA file and parse headers for annotation data
                        faa_file = os.path.join(manicure_dir, f"{sample_id}.faa")
                        if not os.path.exists(faa_file):
    
                            continue
                            

                        faa_data = []
                        
                        with open(faa_file, 'r') as f:
                            for line in f:
                                if line.startswith('>'):
                                    # Parse header: >GCF_964277365.1_iyApoUnic1.Wolbachia_sp_1.1_genomic-----NZ_OZ195518.1_1-----3+1016+1+ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.365
                                    header = line.strip()[1:]  # Remove '>'
                                    parts = header.split('-----')
                                    
                                    if len(parts) >= 3:
                                        sequence_id = parts[1]  # NZ_OZ195518.1_1
                                        coord_part = parts[2]  # 3+1016+1+ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.365
                                        coord_parts = coord_part.split('+')
                                        if len(coord_parts) >= 3:
                                            start = coord_parts[0]
                                            end = coord_parts[1]
                                            strand = '+' if coord_parts[2] == '1' else '-'
                                        else:
                                            start = ''
                                            end = ''
                                            strand = ''
                                        
                                        # Parse annotation attributes
                                        annotation_str = coord_part.split('+', 3)[-1] if len(coord_parts) > 3 else ''
                                        
                                        # Parse annotation attributes
                                        annotation = {}
                                        if ';' in annotation_str:
                                            note_parts = annotation_str.split(';')
                                            for part in note_parts:
                                                if '=' in part:
                                                    key, value = part.split('=', 1)
                                                    annotation[key.strip()] = value.strip()
                                        
                                        faa_data.append({
                                            'sample_id': sample_id,
                                            'sequence_id': sequence_id,
                                            'start': start,
                                            'end': end,
                                            'strand': strand,
                                            'ID': annotation.get('ID', ''),
                                            'partial': annotation.get('partial', ''),
                                            'start_type': annotation.get('start_type', ''),
                                            'rbs_motif': annotation.get('rbs_motif', ''),
                                            'rbs_spacer': annotation.get('rbs_spacer', ''),
                                            'gc_cont': annotation.get('gc_cont', '')
                                        })
                        
                        # 2. Convert to DataFrame
                        faa_df = pd.DataFrame(faa_data)
                        
                        # 3. Check if HMM annotation exists and merge if available
                        hmm_file = os.path.join(annot_dir, f"{sample_id}.hmm.tsv")
                        if os.path.exists(hmm_file):
                            # HMM annotation exists - parse and merge
                            hmm_csv = os.path.join(annot_dir, f"{sample_id}.hmm_clean.csv")
                            # Create temporary ORFCaller instance to use the parse method
                            temp_caller = ORFCaller(output_dir, "fas", temp_args)
                            hmm_rows = temp_caller.parse_hmm_tblout_to_csv(hmm_file, hmm_csv, hmm_domain_evalue_cutoff)
                            if hmm_rows:
                                # Create HMM DataFrame from parsed data
                                hmm_df = pd.DataFrame(hmm_rows)
                                # Extract sequence ID from target_name (first part between dashes)
                                hmm_df['sequence_id'] = hmm_df['target_name'].str.split('-----').str[1]
                                
                                # Ensure both sequence_id columns are strings and handle NaN values
                                faa_df['sequence_id'] = faa_df['sequence_id'].astype(str)
                                hmm_df['sequence_id'] = hmm_df['sequence_id'].astype(str)
                                
                                # 4. Left join on sequence_id
                                # Create a mapping from FAA sequence_id to HMM sequence_id for joining
                                merged_df = pd.merge(faa_df, hmm_df, left_on='sequence_id', right_on='sequence_id', how='left')
                                
                                # Ensure all required HMM columns exist
                                required_hmm_cols = ['target_name', 'query_accession', 'full_evalue', 'full_score', 'full_bias', 
                                                   'dom_evalue', 'dom_score', 'dom_bias', 'exp', 'reg', 'clu', 'ov', 
                                                   'env', 'dom', 'rep', 'inc', 'description']
                                for col in required_hmm_cols:
                                    if col not in merged_df.columns:
                                        merged_df[col] = ''
                            else:
                                # HMM parsing failed, use FAA data with empty HMM columns
                                merged_df = faa_df.copy()
                                empty_hmm_cols = ['query_accession', 'full_evalue', 'full_score', 'full_bias', 
                                                 'dom_evalue', 'dom_score', 'dom_bias', 'exp', 'reg', 'clu', 'ov', 
                                                 'env', 'dom', 'rep', 'inc', 'description']
                                for col in empty_hmm_cols:
                                    merged_df[col] = ''
                        else:
                            # No HMM annotation - just use FAA data with empty HMM columns
                            merged_df = faa_df.copy()
                            empty_hmm_cols = ['query_accession', 'full_evalue', 'full_score', 'full_bias', 
                                             'dom_evalue', 'dom_score', 'dom_bias', 'exp', 'reg', 'clu', 'ov', 
                                             'env', 'dom', 'rep', 'inc', 'description']
                            for col in empty_hmm_cols:
                                merged_df[col] = ''
                        
                        # Drop redundant query_name column before collecting
                        if 'query_name' in merged_df.columns:
                            merged_df = merged_df.drop('query_name', axis=1)
                        
                        # Collect data for later writing
                        if not faa_df.empty:
                            all_merged_data.append(merged_df)
                            logger.info(f"Collected {len(merged_df)} rows for sample {sample_id}")
                
                # Write all data at once with proper header
                if all_merged_data:
                    final_df = pd.concat(all_merged_data, ignore_index=True)
                    final_df.to_csv(summary_file, index=False, sep='\t')
                    logger.info(f"Wrote {len(final_df)} total rows to {summary_file}")
                else:
                    logger.warning(f"No data collected for {subdir}")
        
        # HMM data for bacteria/viruses/metagenomes was already merged inline above
        # No additional processing needed here
        logger.info(f"HMM data for {subdir} was merged inline during summary generation")
        
        logger.info(f"Created comprehensive {subdir} summary: {summary_file}")

def main():
    parser = argparse.ArgumentParser(description='Call ORFs in genomes using a config file or directory search.')
    
    # Config file mode arguments
    parser.add_argument('--config', type=str, default=None, help='Tab-delimited file: sample_id<TAB>genome_path<TAB>domain')
    
    # Directory mode arguments (modeled after dereplicate_genomes.py)
    parser.add_argument('-m', '--mag_dir', type=str, default=None, help='Path or glob to genome files (e.g. asm/*/bins).')
    parser.add_argument('-w', '--wildcard', type=str, default='', help='Pattern to match anywhere in genome file path (can be pipe-separated for multiple patterns).')
    parser.add_argument('--domain', type=str, choices=['bacterial', 'viral', 'eukaryotic', 'metagenomic'], 
                        help='Domain type for all genomes when using directory mode.')
    
    # Common arguments
    parser.add_argument('--output_directory', type=str, default='magus_output/orf_calling', help='Directory to store ORF output files (default: magus_output/orf_calling).')
    parser.add_argument('--max_workers', type=int, default=1, help='Number of ORF calling jobs to run in parallel.')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads for tools (default: 4).')
    parser.add_argument('--extension', type=str, default='fa', help='Extension of genome files (default: fa).')
    parser.add_argument('--force', action='store_true', help='Force rewriting of output files even if they already exist.')
    parser.add_argument('--hmmfile', type=str, default=None, help='Path to HMM file for annotation (optional).')

    parser.add_argument('--annotation-fullseq-evalue', '--afe', type=float, default=1e-2, dest='annotation_fullseq_evalue', help='Filter out HMM hits with full-sequence E-value exceeding this cutoff. Default: 0.01')
    parser.add_argument('--annotation-domain-evalue', '--ade', type=float, default=1e-2, dest='annotation_domain_evalue', help='Filter out HMM hits with best-domain E-value exceeding this cutoff. Default: 0.01')
    parser.add_argument('--suffix', type=str, default=None, help='Suffix to tag outputs (e.g., kegg). Summaries and HMM tblout files include this suffix.')
    parser.add_argument('--eukdb', type=str, default='data/uniref90', help='Path to UniRef90 database for MetaEuk (default: data/uniref90).')
    parser.add_argument('--cleanup', action='store_true', help='Clean up annotation directories after processing.')
    
    # Restart functionality
    parser.add_argument('--restart', type=str, choices=['orf-calling', 'annotation', 'hmm-parsing', 'create-summary'], 
                        help='Restart from specific stage: orf-calling, annotation, hmm-parsing, or create-summary')
    
    args = parser.parse_args()

    # Validate arguments (skip validation for restart modes that don't need genome input)
    if not args.restart:
        if not args.config and not args.mag_dir:
            parser.error("Either --config or --mag_dir must be provided.")
        
        if args.config and args.mag_dir:
            parser.error("Cannot use both --config and --mag_dir. Choose one mode.")
        
        if args.mag_dir and not args.domain:
            parser.error("--domain must be specified when using --mag_dir.")

    os.makedirs(args.output_directory, exist_ok=True)
    orf_caller = ORFCaller(args.output_directory, args.extension, args)

    # Handle restart modes
    if args.restart == 'orf-calling':
        logger.info("Running in orf-calling restart mode - forcing overwrite of existing ORF files")
        args.force = True  # Force overwrite of existing files
        # Continue to normal processing below
        
    elif args.restart == 'annotation':
        logger.info("Running in annotation restart mode - running HMM annotations on existing ORF calls")
        if not args.hmmfile:
            logger.error("--hmmfile must be provided when using --restart annotation")
            return
        
        # Collect all HMM annotation tasks
        annotation_tasks = []
        for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
            annot_dir = os.path.join(args.output_directory, subdir, 'annot')
            if not os.path.exists(annot_dir):
                continue
            
            logger.info(f"Collecting HMM annotation tasks for {subdir}")
            for file in os.listdir(annot_dir):
                # For eukaryotes, use .fas files (protein sequences) but not .codon.fas
                # For others, use .faa files (protein sequences)
                if subdir == 'eukaryotes' and file.endswith('.fas') and not file.endswith('.codon.fas'):
                    sample_id = file.replace('.fas', '')
                    faa_file = os.path.join(annot_dir, file)
                elif subdir != 'eukaryotes' and file.endswith('.faa'):
                    sample_id = file.replace('.faa', '')
                    faa_file = os.path.join(annot_dir, file)
                else:
                    continue
                
                # Check if HMM results already exist
                hmm_file = os.path.join(annot_dir, f"{sample_id}.hmm.tsv")
                if os.path.exists(hmm_file) and not args.force:
                    logger.info(f"HMM results already exist for {sample_id} in {subdir}, skipping")
                    continue
                
                annotation_tasks.append((faa_file, annot_dir, sample_id, args.hmmfile))
        
        # Run HMM annotations in parallel
        if annotation_tasks:
            logger.info(f"Running HMM annotations on {len(annotation_tasks)} files with {args.max_workers} parallel workers")

            def run_single_annotation(task):
                faa_file, annot_dir, sample_id, hmmfile = task
                try:
                    orf_caller = ORFCaller(args.output_directory, args.extension, args)
                    orf_caller.run_hmm_annotation(faa_file, annot_dir, sample_id, hmmfile, args.suffix)
                    return True
                except Exception as e:
                    logger.error(f"Error running HMM annotation for {sample_id}: {e}")
                    return False
            
            if args.max_workers > 1 and len(annotation_tasks) > 1:
                with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
                    results = list(executor.map(run_single_annotation, annotation_tasks))
                successful = sum(results)
                logger.info(f"Successfully annotated {successful}/{len(annotation_tasks)} files")
            else:
                logger.info(f"Running HMM annotations on {len(annotation_tasks)} files sequentially")
                for task in annotation_tasks:
                    run_single_annotation(task)
        else:
            logger.info("No files found that need HMM annotation")
        
        logger.info("Annotation restart completed successfully.")
        return
        
    elif args.restart == 'hmm-parsing':
        logger.info("Running in hmm-parsing restart mode - regenerating clean HMM files")
        
        # Regenerate all HMM clean CSV files
        for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
            annot_dir = os.path.join(args.output_directory, subdir, 'annot')
            if not os.path.exists(annot_dir):
                continue
            
            logger.info(f"Regenerating HMM clean files for {subdir}")
            for file in os.listdir(annot_dir):
                if file.endswith('.hmm.tsv'):
                    sample_id = file.replace('.hmm.tsv', '')
                    hmm_file = os.path.join(annot_dir, file)
                    hmm_csv = os.path.join(annot_dir, f"{sample_id}.hmm_clean.csv")
                    
                    # Create temporary ORFCaller instance
                    temp_caller = ORFCaller(args.output_directory, args.extension, args)
                    logger.info(f"Regenerating clean HMM file for {sample_id}")
                    temp_caller.parse_hmm_tblout_to_csv(hmm_file, hmm_csv, args.annotation_domain_evalue)
        
        logger.info("HMM parsing restart completed successfully.")
        # Continue to create summaries below instead of returning
        
    elif args.restart == 'create-summary':
        logger.info("Running in create-summary restart mode - regenerating comprehensive summaries")
        create_comprehensive_summary(
            args.output_directory, args.hmmfile, args.suffix,
            hmm_fullseq_evalue_cutoff=args.annotation_fullseq_evalue,
            hmm_domain_evalue_cutoff=args.annotation_domain_evalue,
        )
        logger.info("Comprehensive summary creation completed successfully.")
        return

    # Only process genomes if not in restart mode (except orf-calling restart)
    if not args.restart or args.restart == 'orf-calling':
        genomes_data = []

        if args.config:
            # Config file mode
            logger.info(f"Using config file mode: {args.config}")
            with open(args.config, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if not row or row[0].startswith('#') or len(row) < 3:
                        continue
                    sample_id, genome_path, domain = row[0], row[1], row[2]
                    genomes_data.append((sample_id, genome_path, domain))
        
        else:
            # Directory mode
            logger.info(f"Using directory mode: {args.mag_dir}")
            
            # Handle multiple wildcard patterns separated by pipes
            wildcards = args.wildcard.split('|') if args.wildcard else ['']
            
            all_genome_files = []
            for wildcard in wildcards:
                try:
                    genome_files = find_genome_files(args.mag_dir, args.extension, wildcard.strip())
                    all_genome_files.extend(genome_files)
                    logger.info(f"Found {len(genome_files)} files matching wildcard '{wildcard.strip()}'")
                except RuntimeError as e:
                    if len(wildcards) == 1:  # Only one wildcard, so this is an error
                        raise e
                    else:  # Multiple wildcards, so this one just didn't match anything
                        logger.warning(f"No files found for wildcard '{wildcard.strip()}'")
            
            # Remove duplicates while preserving order
            seen = set()
            unique_genome_files = []
            for path in all_genome_files:
                if path not in seen:
                    seen.add(path)
                    unique_genome_files.append(path)
            
            if not unique_genome_files:
                raise RuntimeError("No matching genome files found with any of the provided wildcards.")
            
            logger.info(f"Total unique genome files found: {len(unique_genome_files)}")
            
            # Create genome data tuples - ALL files get the same domain specified at command line
            for genome_path in unique_genome_files:
                # Use the filename (without extension) as sample_id
                sample_id = genome_path.name
                # Remove the extension properly
                extension = f".{args.extension}" if not args.extension.startswith('.') else args.extension
                if sample_id.endswith(extension):
                    sample_id = sample_id[:-len(extension)]
                genomes_data.append((sample_id, str(genome_path), args.domain))

        # Process all genomes (Stage 1: ORF calling)
        logger.info("Stage 1: Calling ORFs")
        process_genomes(orf_caller, genomes_data, args.max_workers)

    # Stage 2: Create comprehensive summaries with HMM results
    logger.info("Stage 2: Creating comprehensive summaries")
    create_comprehensive_summary(
        args.output_directory, args.hmmfile, args.suffix,
        hmm_fullseq_evalue_cutoff=args.annotation_fullseq_evalue,
        hmm_domain_evalue_cutoff=args.annotation_domain_evalue,
    )

    # Clean up annot directories
    for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
        annot_dir = os.path.join(args.output_directory, subdir, 'annot')
        if os.path.exists(annot_dir) and args.cleanup:
            import shutil
            shutil.rmtree(annot_dir)

    logger.info("ORF calling completed successfully.")

if __name__ == '__main__':
    main() 
