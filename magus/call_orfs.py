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

class ORFCaller:
    def __init__(self, output_dir, extension, args):
        self.output_dir = output_dir
        self.extension = extension
        self.args = args

    def run_hmm_annotation(self, faa_file, output_dir, sample_id, hmmfile):
        if hmmfile and os.path.exists(faa_file):
            hmm_out = os.path.join(output_dir, f"{sample_id}.hmm.tsv")
            cmd = [
                'hmmsearch',
                '--tblout', hmm_out,
                '--noali',
                '--notextw',
                hmmfile,
                faa_file
            ]
            logger.info(f"Running HMM annotation for {sample_id}")
            with open(os.path.join(output_dir, f"{sample_id}_hmmsearch.log"), 'w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)

    def call_bacterial_orfs(self, genome_file, sample_id=None, hmmfile=None):
        """Call ORFs in bacterial genomes using prodigal."""
        annot_dir = os.path.join(self.output_dir, 'bacteria', 'annot')
        manicure_dir = os.path.join(self.output_dir, 'bacteria', 'manicure')
        os.makedirs(annot_dir, exist_ok=True)
        os.makedirs(manicure_dir, exist_ok=True)

        FN = sample_id if sample_id else os.path.basename(genome_file).replace(self.extension, '')
        manicure_file = os.path.join(manicure_dir, f"{FN}.faa")

        # Check if the manicure file already exists and has nonzero length
        if os.path.exists(manicure_file) and os.path.getsize(manicure_file) > 0 and not self.args.force:
            logger.info(f"Skipping ORF calling for {genome_file} as output already exists.")
        else:
            log_file = os.path.join(self.output_dir, 'bacteria', f"{FN}_prodigal.log")
            cmd = ['prodigal', '-i', genome_file, '-d', os.path.join(annot_dir, f"{FN}.ffn"), '-a', os.path.join(annot_dir, f"{FN}.faa"), '-o', '/dev/null']
            logger.info(f"Calling bacterial ORFs for {genome_file} using prodigal")
            with open(log_file, 'w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)

            # Manicure the output files
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

            # Clean up the annot files
            os.remove(os.path.join(annot_dir, f"{FN}.ffn"))
            os.remove(os.path.join(annot_dir, f"{FN}.faa"))

        # HMM annotation (standalone)
        self.run_hmm_annotation(manicure_file, manicure_dir, FN, hmmfile)
        return manicure_file

    def call_viral_orfs(self, genome_file, sample_id=None, hmmfile=None):
        """Call ORFs in viral genomes using prodigal-gv."""
        annot_dir = os.path.join(self.output_dir, 'viruses', 'annot')
        manicure_dir = os.path.join(self.output_dir, 'viruses', 'manicure')
        os.makedirs(annot_dir, exist_ok=True)
        os.makedirs(manicure_dir, exist_ok=True)

        FN = sample_id if sample_id else os.path.basename(genome_file).replace(self.extension, '')
        manicure_file = os.path.join(manicure_dir, f"{FN}.faa")

        # Check if the manicure file already exists and has nonzero length
        if os.path.exists(manicure_file) and os.path.getsize(manicure_file) > 0 and not self.args.force:
            logger.info(f"Skipping ORF calling for {genome_file} as output already exists.")
            return manicure_file

        log_file = os.path.join(self.output_dir, 'viruses', f"{FN}_prodigal.log")
        cmd = ['prodigal-gv', '-p', '-q', '-i', genome_file, '-d', os.path.join(annot_dir, f"{FN}.ffn"), '-a', os.path.join(annot_dir, f"{FN}.faa"), '-o', '/dev/null']
        logger.info(f"Calling viral ORFs for {genome_file} using prodigal-gv")
        with open(log_file, 'w') as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=log)

        # Manicure the output files
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

        # Clean up the annot files
        os.remove(os.path.join(annot_dir, f"{FN}.ffn"))
        os.remove(os.path.join(annot_dir, f"{FN}.faa"))

        # HMM annotation (standalone)
        self.run_hmm_annotation(manicure_file, manicure_dir, FN, hmmfile)
        return manicure_file

    def call_eukaryotic_orfs(self, genome_file, sample_id=None, hmmfile=None):
        """Call ORFs in eukaryotic genomes using MetaEuk."""
        annot_dir = os.path.join(self.output_dir, 'eukaryotes', 'annot')
        manicure_dir = os.path.join(self.output_dir, 'eukaryotes', 'manicure')
        os.makedirs(annot_dir, exist_ok=True)
        os.makedirs(manicure_dir, exist_ok=True)

        FN = sample_id if sample_id else os.path.basename(genome_file).replace(self.extension, '')
        manicure_file = os.path.join(manicure_dir, f"{FN}.faa")

        # Check if the manicure file already exists and has nonzero length
        if os.path.exists(manicure_file) and os.path.getsize(manicure_file) > 0 and not self.args.force:
            logger.info(f"Skipping ORF calling for {genome_file} as output already exists.")
            return manicure_file

        log_file = os.path.join(self.output_dir, 'eukaryotes', f"{FN}_metaeuk.log")
        cmd = ['metaeuk', 'easy-predict', genome_file, 'data/uniref90', os.path.join(annot_dir, f"{FN}"), os.path.join(annot_dir, f"{FN}")]
        logger.info(f"Calling eukaryotic ORFs for {genome_file} using MetaEuk")
        with open(log_file, 'w') as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=log)

        # Manicure the output files
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

        # Clean up the annot files
        os.remove(os.path.join(annot_dir, f"{FN}.ffn"))
        os.remove(os.path.join(annot_dir, f"{FN}.faa"))

        # HMM annotation (standalone)
        self.run_hmm_annotation(manicure_file, manicure_dir, FN, hmmfile)
        return manicure_file

    def call_metagenome_orfs(self, genome_file, sample_id=None, hmmfile=None):
        """Call ORFs in metagenome mode using prodigal-gv."""
        annot_dir = os.path.join(self.output_dir, 'metagenomes', 'annot')
        manicure_dir = os.path.join(self.output_dir, 'metagenomes', 'manicure')
        os.makedirs(annot_dir, exist_ok=True)
        os.makedirs(manicure_dir, exist_ok=True)

        FN = sample_id if sample_id else os.path.basename(genome_file).replace(self.extension, '')
        manicure_file = os.path.join(manicure_dir, f"{FN}.faa")

        # Check if the manicure file already exists and has nonzero length
        if os.path.exists(manicure_file) and os.path.getsize(manicure_file) > 0 and not self.args.force:
            logger.info(f"Skipping ORF calling for {genome_file} as output already exists.")
            return manicure_file

        log_file = os.path.join(self.output_dir, 'metagenomes', f"{FN}_prodigal.log")
        cmd = ['prodigal', '-p', 'meta', '-q', '-i', genome_file, '-d', os.path.join(annot_dir, f"{FN}.ffn"), '-a', os.path.join(annot_dir, f"{FN}.faa"), '-o', '/dev/null']
        logger.info(f"Calling metagenome ORFs for {genome_file} using prodigal-gv")
        with open(log_file, 'w') as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=log)

        # Manicure the output files
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

        # Clean up the annot files
        os.remove(os.path.join(annot_dir, f"{FN}.ffn"))
        os.remove(os.path.join(annot_dir, f"{FN}.faa"))

        # HMM annotation (standalone)
        self.run_hmm_annotation(manicure_file, manicure_dir, FN, hmmfile)
        return manicure_file

def summarize_output(output_dir):
    """Summarize the output in a table with specified columns."""
    summary_file = os.path.join(output_dir, 'summary.txt')
    with open(summary_file, 'w') as summary:
        summary.write('sample_id\tcontig_id\tstart\tend\tstrand\tID\tpartial\tstart_type\trbs_motif\trbs_spacer\tgc_cont\tmode\n')
        for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
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
    parser = argparse.ArgumentParser(description='Call ORFs in bacterial, viral, and eukaryotic genomes.')
    parser.add_argument('--bacterial_genomes', type=str, help='Path to folder with bacterial genomes.')
    parser.add_argument('--viral_genomes', type=str, help='Path to folder with viral genomes.')
    parser.add_argument('--eukaryotic_genomes', type=str, help='Path to folder with eukaryotic genomes.')
    parser.add_argument('--metagenomes', type=str, help='Path to folder with metagenomic contigs.')
    parser.add_argument('--max_workers', type=int, default=1, help='Number of ORF calling jobs to run in parallel.')
    parser.add_argument('--extension', type=str, default='.fa', help='Extension of genome files (default: .fa).')
    parser.add_argument('--output_directory', type=str, default='magus_output/orf_calling', help='Directory to store ORF output files (default: magus_output/orf_calling).')
    parser.add_argument('--force', action='store_true', help='Force rewriting of output files even if they already exist.')
    parser.add_argument('--hmmfile', type=str, default=None, help='Path to HMM file for annotation (optional).')
    parser.add_argument('--config_file', type=str, default=None, help='Tab-delimited file: sample_id<TAB>genome_path')
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_directory, exist_ok=True)

    # Initialize ORFCaller
    orf_caller = ORFCaller(args.output_directory, args.extension, args)

    if args.config_file:
        with open(args.config_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if not row or row[0].startswith('#') or len(row) < 2:
                    continue
                sample_id, genome_path = row[0], row[1]
                orf_caller.call_bacterial_orfs(genome_path, sample_id=sample_id, hmmfile=args.hmmfile)
    else:
        # Get list of genome files from the input directories
        bacterial_files = [os.path.join(args.bacterial_genomes, f) for f in os.listdir(args.bacterial_genomes) if f.endswith(args.extension)] if args.bacterial_genomes else []
        viral_files = [os.path.join(args.viral_genomes, f) for f in os.listdir(args.viral_genomes) if f.endswith(args.extension)] if args.viral_genomes else []
        eukaryotic_files = [os.path.join(args.eukaryotic_genomes, f) for f in os.listdir(args.eukaryotic_genomes) if f.endswith(args.extension)] if args.eukaryotic_genomes else []
        metagenome_files = [os.path.join(args.metagenomes, f) for f in os.listdir(args.metagenomes) if f.endswith(args.extension)] if args.metagenomes else []

        # Process genomes based on type
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            for genome_file in bacterial_files:
                executor.submit(orf_caller.call_bacterial_orfs, genome_file, None, args.hmmfile)
            for genome_file in viral_files:
                executor.submit(orf_caller.call_viral_orfs, genome_file, None, args.hmmfile)
            for genome_file in eukaryotic_files:
                executor.submit(orf_caller.call_eukaryotic_orfs, genome_file, None, args.hmmfile)
            for genome_file in metagenome_files:
                executor.submit(orf_caller.call_metagenome_orfs, genome_file, None, args.hmmfile)

    # Remove the annot directory after the script is complete
    for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
        annot_dir = os.path.join(args.output_directory, subdir, 'annot')
        if os.path.exists(annot_dir):
            import shutil
            shutil.rmtree(annot_dir)

    # Summarize the output
    summarize_output(args.output_directory)

    logger.info("ORF calling completed successfully.")

if __name__ == '__main__':
    main() 
