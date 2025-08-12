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
            cmd = ['prodigal', '-i', genome_file, '-d', ffn_file, '-a', faa_file, '-o', gff_file]
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
            cmd = ['prodigal-gv', '-p', '-q', '-i', genome_file, '-d', ffn_file, '-a', faa_file, '-o', gff_file]
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
            self.run_hmm_annotation(faa_file, annot_dir, FN, self.args.hmmfile)
        
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
            cmd = ['prodigal', '-p', 'meta', '-q', '-i', genome_file, '-d', ffn_file, '-a', faa_file, '-o', gff_file]
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

def create_comprehensive_summary(output_dir, hmmfile):
    """Create ONE comprehensive summary file per domain with ALL information."""
    for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
        annot_dir = os.path.join(output_dir, subdir, 'annot')
        if not os.path.exists(annot_dir):
            continue
            
        summary_file = os.path.join(output_dir, f'{subdir}_orf_summary.tsv')
        logger.info(f"Creating comprehensive summary for {subdir}: {summary_file}")
        
        with open(summary_file, 'w') as summary:
            if subdir == 'eukaryotes':
                # Eukaryotes: use headersMap.tsv files
                header_written = False
                for file in os.listdir(annot_dir):
                    if file.endswith('.headersMap.tsv'):
                        sample_id = file.replace('.headersMap.tsv', '')
                        file_path = os.path.join(annot_dir, file)
                        
                        with open(file_path, 'r') as infile:
                            for i, line in enumerate(infile):
                                line = line.rstrip('\n')
                                if not line.strip():
                                    continue
                                
                                # Write header only once
                                if i == 0 and not header_written:
                                    summary.write(f'sample_id\t{line}\n')
                                    header_written = True
                                elif i > 0 or header_written:
                                    summary.write(f'{sample_id}\t{line}\n')
            else:
                # Bacteria/Viruses/Metagenomes: parse Prodigal output
                summary.write('sample_id\tcontig_id\tstart\tend\tstrand\tpartial\tstart_type\trbs_motif\trbs_spacer\tgc_cont\n')
                
                for file in os.listdir(annot_dir):
                    if file.endswith('.faa'):
                        sample_id = file.replace('.faa', '')
                        faa_file = os.path.join(annot_dir, file)
                        
                        with open(faa_file, 'r') as infile:
                            for line in infile:
                                if line.startswith('>'):
                                    # Parse Prodigal header: >contig_1 # 1 # 279 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.500
                                    parts = line.strip().split(' # ')
                                    if len(parts) >= 4:
                                        contig_id = parts[0].replace('>', '')
                                        start = parts[1]
                                        end = parts[2]
                                        strand = parts[3]
                                        
                                        # Parse the metadata part
                                        metadata = parts[4] if len(parts) > 4 else ''
                                        partial = '00'
                                        start_type = 'ATG'
                                        rbs_motif = 'None'
                                        rbs_spacer = 'None'
                                        gc_cont = '0.500'
                                        
                                        if 'partial=' in metadata:
                                            partial = metadata.split('partial=')[1].split(';')[0]
                                        if 'start_type=' in metadata:
                                            start_type = metadata.split('start_type=')[1].split(';')[0]
                                        if 'rbs_motif=' in metadata:
                                            rbs_motif = metadata.split('rbs_motif=')[1].split(';')[0]
                                        if 'rbs_spacer=' in metadata:
                                            rbs_spacer = metadata.split('rbs_spacer=')[1].split(';')[0]
                                        if 'gc_cont=' in metadata:
                                            gc_cont = metadata.split('gc_cont=')[1].split(';')[0]
                                        
                                        summary.write(f'{sample_id}\t{contig_id}\t{start}\t{end}\t{strand}\t{partial}\t{start_type}\t{rbs_motif}\t{rbs_spacer}\t{gc_cont}\n')
        
        # Now add HMM results if available
        if hmmfile:
            logger.info(f"Adding HMM results to {subdir} summary")
            hmm_results = {}
            
            # Collect HMM results from all samples
            for file in os.listdir(annot_dir):
                if file.endswith('.hmm.tsv'):
                    sample_id = file.replace('.hmm.tsv', '')
                    hmm_file = os.path.join(annot_dir, file)
                    
                    with open(hmm_file, 'r') as infile:
                        for line in infile:
                            if line.startswith('#'):
                                continue
                            parts = line.strip().split()
                            if len(parts) >= 5:
                                target_name = parts[0]
                                query_name = parts[2]
                                evalue = parts[4]
                                score = parts[5]
                                
                                if query_name not in hmm_results:
                                    hmm_results[query_name] = []
                                hmm_results[query_name].append(f"{target_name}:{evalue}:{score}")
            
            # Add HMM column to summary
            if hmm_results:
                # Read existing summary and add HMM column
                temp_summary = summary_file + '.tmp'
                with open(summary_file, 'r') as infile, open(temp_summary, 'w') as outfile:
                    header = infile.readline().strip()
                    outfile.write(f'{header}\thmm_hits\n')
                    
                    for line in infile:
                        line = line.strip()
                        if not line:
                            continue
                        
                        # Extract query name from line to match with HMM results
                        parts = line.split('\t')
                        if subdir == 'eukaryotes':
                            # For eukaryotes, need to extract query name from headersMap format
                            # This depends on the actual format of your headersMap.tsv
                            query_name = parts[1] if len(parts) > 1 else 'unknown'
                        else:
                            # For bacteria/viruses, query is contig_id
                            query_name = parts[1] if len(parts) > 1 else 'unknown'
                        
                        hmm_hits = ';'.join(hmm_results.get(query_name, []))
                        outfile.write(f'{line}\t{hmm_hits}\n')
                
                # Replace original with enhanced version
                os.replace(temp_summary, summary_file)
        
        logger.info(f"Created comprehensive {subdir} summary: {summary_file}")

# This function has been removed as it was redundant with create_comprehensive_summary

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
    parser.add_argument('--eukdb', type=str, default='data/uniref90', help='Path to UniRef90 database for MetaEuk (default: data/uniref90).')
    parser.add_argument('--cleanup', action='store_true', help='Clean up annotation directories after processing.')
    
    # Restart functionality
    parser.add_argument('--restart', type=str, choices=['create-summary', 'annotations'], 
                        help='Restart from specific stage: create-summary or annotations')
    
    args = parser.parse_args()

    # Validate arguments
    if not args.config and not args.mag_dir:
        parser.error("Either --config or --mag_dir must be provided.")
    
    if args.config and args.mag_dir:
        parser.error("Cannot use both --config and --mag_dir. Choose one mode.")
    
    if args.mag_dir and not args.domain:
        parser.error("--domain must be specified when using --mag_dir.")

    os.makedirs(args.output_directory, exist_ok=True)
    orf_caller = ORFCaller(args.output_directory, args.extension, args)

    # Handle restart modes
    if args.restart == 'create-summary':
        logger.info("Running in create-summary mode - generating comprehensive summaries")
        create_comprehensive_summary(args.output_directory, args.hmmfile)
        logger.info("Comprehensive summary creation completed successfully.")
        return
    elif args.restart == 'annotations':
        logger.info("Running in annotations mode - running HMM annotations on existing ORF calls")
        if not args.hmmfile:
            logger.error("--hmmfile must be provided when using --restart annotations")
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
                    orf_caller.run_hmm_annotation(faa_file, annot_dir, sample_id, hmmfile)
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
        
        # Create comprehensive summaries with new HMM results
        logger.info("Creating comprehensive summaries with updated HMM results")
        create_comprehensive_summary(args.output_directory, args.hmmfile)
        logger.info("Annotations restart completed successfully.")
        return

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
    create_comprehensive_summary(args.output_directory, args.hmmfile)

    # Clean up annot directories
    for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
        annot_dir = os.path.join(args.output_directory, subdir, 'annot')
        if os.path.exists(annot_dir) and args.cleanup:
            import shutil
            shutil.rmtree(annot_dir)

    logger.info("ORF calling completed successfully.")

if __name__ == '__main__':
    main() 
