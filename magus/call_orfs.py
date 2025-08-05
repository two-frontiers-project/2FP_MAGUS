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
import pandas as pd

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
        os.makedirs(annot_dir, exist_ok=True)

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
            cmd = ['prodigal', '-i', genome_file, '-d', ffn_file, '-a', faa_file, '-o', gff_file]
            logger.info(f"Calling bacterial ORFs for {genome_file} using prodigal")
            with open(log_file, 'w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)
        
        return faa_file, ffn_file, gff_file

    def call_viral_orfs(self, genome_file, sample_id=None):
        """Call ORFs in viral genomes using prodigal-gv."""
        annot_dir = os.path.join(self.output_dir, 'viruses', 'annot')
        os.makedirs(annot_dir, exist_ok=True)

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
            cmd = ['prodigal-gv', '-p', '-q', '-i', genome_file, '-d', ffn_file, '-a', faa_file, '-o', gff_file]
            logger.info(f"Calling viral ORFs for {genome_file} using prodigal-gv")
            with open(log_file, 'w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)
        
        return faa_file, ffn_file, gff_file

    def call_eukaryotic_orfs(self, genome_file, sample_id=None):
        """Call ORFs in eukaryotic genomes using MetaEuk."""
        annot_dir = os.path.join(self.output_dir, 'eukaryotes', 'annot')
        os.makedirs(annot_dir, exist_ok=True)

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
        
        return faa_file, ffn_file, None  # MetaEuk doesn't produce GFF

    def call_metagenome_orfs(self, genome_file, sample_id=None):
        """Call ORFs in metagenome mode using prodigal."""
        annot_dir = os.path.join(self.output_dir, 'metagenomes', 'annot')
        os.makedirs(annot_dir, exist_ok=True)

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
            cmd = ['prodigal', '-p', 'meta', '-q', '-i', genome_file, '-d', ffn_file, '-a', faa_file, '-o', gff_file]
            logger.info(f"Calling metagenome ORFs for {genome_file} using prodigal (metagenome mode)")
            with open(log_file, 'w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)
        
        return faa_file, ffn_file, gff_file

    def run_hmm_annotation(self, faa_file, output_dir, sample_id, hmmfile):
        """Run HMM annotation on a protein file."""
        if not hmmfile:
            logger.warning(f"No HMM file provided for {sample_id}. Skipping annotation.")
            return None
            
        if not os.path.exists(faa_file):
            logger.warning(f"Protein file {faa_file} does not exist. Skipping annotation.")
            return None
            
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

def process_genomes(orf_caller, genomes_data):
    """Process a list of genomes data (tuples of sample_id, genome_path, domain)"""
    for sample_id, genome_path, domain in genomes_data:
        domain = domain.lower()
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

def summarize_calls(output_dir):
    """Summarize ORF calls into a single TSV per domain type."""
    for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
        annot_dir = os.path.join(output_dir, subdir, 'annot')
        if not os.path.exists(annot_dir):
            continue
            
        summary_file = os.path.join(output_dir, f'{subdir}_calls_summary.tsv')
        all_calls = []
        
        for file in os.listdir(annot_dir):
            if file.endswith('.faa') or file.endswith('.fas'):
                sample_id = file.replace('.faa', '').replace('.fas', '')
                faa_file = os.path.join(annot_dir, file)
                
                with open(faa_file, 'r') as infile:
                    for line in infile:
                        if line.startswith('>'):
                            # Parse header and extract information
                            header = line.strip()
                            # Add sample_id as first column
                            call_info = [sample_id, header]
                            all_calls.append(call_info)
        
        # Write summary file
        with open(summary_file, 'w') as summary:
            summary.write('sample_id\torf_header\n')
            for call in all_calls:
                summary.write(f'{call[0]}\t{call[1]}\n')
        
        logger.info(f"Created {subdir} calls summary: {summary_file}")

def run_annotations(output_dir, hmmfile):
    """Run HMM annotations on all protein files."""
    if not hmmfile:
        logger.warning("No HMM file provided. Skipping annotation stage.")
        return
        
    for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
        annot_dir = os.path.join(output_dir, subdir, 'annot')
        if not os.path.exists(annot_dir):
            continue
            
        for file in os.listdir(annot_dir):
            if file.endswith('.faa') or file.endswith('.fas'):
                sample_id = file.replace('.faa', '').replace('.fas', '')
                faa_file = os.path.join(annot_dir, file)
                
                orf_caller = ORFCaller(output_dir, 'fa', type('Args', (), {'threads': 4, 'force': False})())
                orf_caller.run_hmm_annotation(faa_file, annot_dir, sample_id, hmmfile)

def summarize_annotations(output_dir):
    """Create manicured files and final summary with annotations."""
    for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
        annot_dir = os.path.join(output_dir, subdir, 'annot')
        manicure_dir = os.path.join(output_dir, subdir, 'manicure')
        os.makedirs(manicure_dir, exist_ok=True)
        
        if not os.path.exists(annot_dir):
            continue
            
        for file in os.listdir(annot_dir):
            if file.endswith('.faa') or file.endswith('.fas'):
                sample_id = file.replace('.faa', '').replace('.fas', '')
                faa_file = os.path.join(annot_dir, file)
                manicure_file = os.path.join(manicure_dir, f"{sample_id}.faa")
                
                # Manicure the protein file
                with open(faa_file, 'r') as infile, open(manicure_file, 'w') as outfile:
                    for line in infile:
                        if line.startswith('>'):
                            if subdir == 'eukaryotes':
                                # Handle MetaEuk header format: >UniRef90_ID|sample|strand|score|evalue|num_exons|start|end|exon_info
                                parts = line.strip().split('|')
                                if len(parts) >= 7:
                                    uniref_id = parts[0].replace('>', '')
                                    contig_id = parts[1]
                                    strand = parts[2]
                                    start = parts[6]
                                    end = parts[7]
                                    # Create simplified header format for MetaEuk
                                    new_header = f">{sample_id}-----{contig_id}-----{start}+{end}+{strand}+ID={uniref_id};partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.500\n"
                                    outfile.write(new_header)
                                else:
                                    outfile.write(line)
                            else:
                                # Handle Prodigal header format
                                outfile.write(line.replace('>',f'>{sample_id}-----').replace(' # ', '-----', 1).replace('*', 'X').replace(' # ', '+').replace(' # ', '-'))
                        else:
                            outfile.write(line)
                
                # Also manicure the nucleotide file if it exists
                if file.endswith('.faa'):
                    ffn_file = faa_file.replace('.faa', '.ffn')
                    if os.path.exists(ffn_file):
                        manicure_ffn = os.path.join(manicure_dir, f"{sample_id}.ffn")
                        with open(ffn_file, 'r') as infile, open(manicure_ffn, 'w') as outfile:
                            for line in infile:
                                if line.startswith('>'):
                                    outfile.write(line.replace('>',f'>{sample_id}-----').replace('+', '-----+').replace('*', 'X').replace(' # ', '+').replace(' # ', '-'))
                                else:
                                    outfile.write(line)
                elif file.endswith('.fas'):
                    ffn_file = faa_file.replace('.fas', '.codon.fas')
                    if os.path.exists(ffn_file):
                        manicure_ffn = os.path.join(manicure_dir, f"{sample_id}.ffn")
                        with open(ffn_file, 'r') as infile, open(manicure_ffn, 'w') as outfile:
                            for line in infile:
                                if line.startswith('>'):
                                    if subdir == 'eukaryotes':
                                        # Handle MetaEuk header format for nucleotide sequences
                                        parts = line.strip().split('|')
                                        if len(parts) >= 7:
                                            uniref_id = parts[0].replace('>', '')
                                            contig_id = parts[1]
                                            strand = parts[2]
                                            start = parts[6]
                                            end = parts[7]
                                            # Create simplified header format for MetaEuk
                                            new_header = f">{sample_id}-----{contig_id}-----{start}+{end}+{strand}+ID={uniref_id};partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.500\n"
                                            outfile.write(new_header)
                                        else:
                                            outfile.write(line)
                                    else:
                                        outfile.write(line.replace('>',f'>{sample_id}-----').replace('+', '-----+').replace('*', 'X').replace(' # ', '+').replace(' # ', '-'))
                                else:
                                    outfile.write(line)
        
        # Create final summary for this domain
        summary_file = os.path.join(output_dir, f'{subdir}_final_summary.tsv')
        with open(summary_file, 'w') as summary:
            summary.write('sample_id\tcontig_id\tstart\tend\tstrand\tID\tpartial\tstart_type\trbs_motif\trbs_spacer\tgc_cont\tmode\n')
            manicure_dir = os.path.join(output_dir, subdir, 'manicure')
            if os.path.exists(manicure_dir):
                for file in os.listdir(manicure_dir):
                    if file.endswith('.faa'):
                        with open(os.path.join(manicure_dir, file), 'r') as infile:
                            for line in infile:
                                if line.startswith('>'):
                                    parts = line.strip().split('-----')
                                    if len(parts) >= 3:
                                        sample_id = parts[0].replace('>', '')
                                        contig_id = parts[1]
                                        metadata = parts[2]
                                        metadata_parts = metadata.split('+')
                                        if len(metadata_parts) >= 4:
                                            start = metadata_parts[0]
                                            end = metadata_parts[1]
                                            strand = metadata_parts[2]
                                            additional_info = metadata_parts[3].split(';')
                                            id_part = additional_info[0].split('=')[1]
                                            partial_part = additional_info[1].split('=')[1]
                                            start_type_part = additional_info[2].split('=')[1]
                                            rbs_motif_part = additional_info[3].split('=')[1]
                                            rbs_spacer_part = additional_info[4].split('=')[1]
                                            gc_cont_part = additional_info[5].split('=')[1]
                                            mode = subdir
                                            summary.write(f'{sample_id}\t{contig_id}\t{start}\t{end}\t{strand}\t{id_part}\t{partial_part}\t{start_type_part}\t{rbs_motif_part}\t{rbs_spacer_part}\t{gc_cont_part}\t{mode}\n')

def main():
    parser = argparse.ArgumentParser(description='Call ORFs in genomes using a config file or directory search.')
    
    # Config file mode arguments
    parser.add_argument('--config_file', type=str, default=None, help='Tab-delimited file: sample_id<TAB>genome_path<TAB>domain')
    
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
    parser.add_argument('--eukdb', type=str, default='data/uniref90', help='Path to UniRef90 database for MetaEuk (default: data/uniref90).')
    parser.add_argument('--cleanup', action='store_true', help='Clean up annotation directories after processing.')
    
    # Restart functionality
    parser.add_argument('--restart', type=str, choices=['summarize-calls', 'annotation', 'summarize-annotations'], 
                        help='Restart from specific stage: summarize-calls, annotation, or summarize-annotations')
    
    args = parser.parse_args()

    # Validate arguments
    if not args.config_file and not args.mag_dir:
        parser.error("Either --config_file or --mag_dir must be provided.")
    
    if args.config_file and args.mag_dir:
        parser.error("Cannot use both --config_file and --mag_dir. Choose one mode.")
    
    if args.mag_dir and not args.domain:
        parser.error("--domain must be specified when using --mag_dir.")

    os.makedirs(args.output_directory, exist_ok=True)
    orf_caller = ORFCaller(args.output_directory, args.extension, args)

    # Handle restart modes
    if args.restart == 'summarize-calls':
        logger.info("Running in summarize-calls mode - generating call summaries")
        summarize_calls(args.output_directory)
        logger.info("Call summarization completed successfully.")
        return
    elif args.restart == 'annotation':
        logger.info("Running in annotation mode - running HMM annotations")
        run_annotations(args.output_directory, args.hmmfile)
        logger.info("Annotation completed successfully.")
        return
    elif args.restart == 'summarize-annotations':
        logger.info("Running in summarize-annotations mode - creating manicured files and final summaries")
        summarize_annotations(args.output_directory)
        logger.info("Annotation summarization completed successfully.")
        return

    genomes_data = []

    if args.config_file:
        # Config file mode
        logger.info(f"Using config file mode: {args.config_file}")
        with open(args.config_file, 'r') as f:
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
    process_genomes(orf_caller, genomes_data)

    # Stage 2: Summarize calls
    logger.info("Stage 2: Summarizing calls")
    summarize_calls(args.output_directory)

    # Stage 3: Run annotations (if HMM file provided)
    logger.info("Stage 3: Running annotations")
    run_annotations(args.output_directory, args.hmmfile)

    # Stage 4: Summarize annotations
    logger.info("Stage 4: Summarizing annotations")
    summarize_annotations(args.output_directory)

    # Clean up annot directories
    for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
        annot_dir = os.path.join(args.output_directory, subdir, 'annot')
        if os.path.exists(annot_dir) and args.cleanup:
            import shutil
            shutil.rmtree(annot_dir)

    logger.info("ORF calling completed successfully.")

if __name__ == '__main__':
    main() 
