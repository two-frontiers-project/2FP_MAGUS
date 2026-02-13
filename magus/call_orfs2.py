#!/usr/bin/env python3

import os
import argparse
import logging
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import csv
import glob

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ORFCaller:
    def __init__(self, output_dir, extension, threads, force, eukdb):
        self.output_dir = output_dir
        self.extension = f".{extension}" if not extension.startswith('.') else extension
        self.threads = threads
        self.force = force
        self.eukdb = eukdb

    def call_bacterial_orfs(self, genome_file, sample_id=None):
        annot_dir = os.path.join(self.output_dir, 'bacteria', 'annot')
        manicure_dir = os.path.join(self.output_dir, 'bacteria', 'manicure')
        os.makedirs(annot_dir, exist_ok=True)
        os.makedirs(manicure_dir, exist_ok=True)

        FN = sample_id if sample_id else os.path.basename(genome_file).replace(self.extension, '')

        faa_file = os.path.join(annot_dir, f"{FN}.faa")
        ffn_file = os.path.join(annot_dir, f"{FN}.ffn")
        gff_file = os.path.join(annot_dir, f"{FN}.gff")

        if (os.path.exists(faa_file) and os.path.getsize(faa_file) > 0 and
            os.path.exists(ffn_file) and os.path.getsize(ffn_file) > 0 and
            os.path.exists(gff_file) and os.path.getsize(gff_file) > 0 and
            not self.force):
            logger.info(f"Skipping ORF calling for {genome_file} as output already exists.")
        else:
            log_file = os.path.join(self.output_dir, 'bacteria', f"{FN}_prodigal.log")
            cmd = ['prodigal', '-i', genome_file, '-d', ffn_file, '-a', faa_file, '-o', gff_file, '-f', 'gff']
            with open(log_file, 'w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)

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

        sample_tmp_dir = os.path.join(annot_dir, FN)
        if os.path.isdir(sample_tmp_dir):
            import shutil
            shutil.rmtree(sample_tmp_dir, ignore_errors=True)
        return True

    def call_viral_orfs(self, genome_file, sample_id=None):
        annot_dir = os.path.join(self.output_dir, 'viruses', 'annot')
        manicure_dir = os.path.join(self.output_dir, 'viruses', 'manicure')
        os.makedirs(annot_dir, exist_ok=True)
        os.makedirs(manicure_dir, exist_ok=True)

        FN = sample_id if sample_id else os.path.basename(genome_file).replace(self.extension, '')

        faa_file = os.path.join(annot_dir, f"{FN}.faa")
        ffn_file = os.path.join(annot_dir, f"{FN}.ffn")
        gff_file = os.path.join(annot_dir, f"{FN}.gff")

        if (os.path.exists(faa_file) and os.path.getsize(faa_file) > 0 and
            os.path.exists(ffn_file) and os.path.getsize(ffn_file) > 0 and
            os.path.exists(gff_file) and os.path.getsize(gff_file) > 0 and
            not self.force):
            logger.info(f"Skipping ORF calling for {genome_file} as output already exists.")
        else:
            log_file = os.path.join(self.output_dir, 'viruses', f"{FN}_prodigal.log")
            cmd = ['prodigal-gv', '-p', 'meta', '-q', '-i', genome_file, '-d', ffn_file, '-a', faa_file, '-o', gff_file, '-f', 'gff']
            with open(log_file, 'w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)

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
        return True

    def call_eukaryotic_orfs(self, genome_file, sample_id=None):
        annot_dir = os.path.join(self.output_dir, 'eukaryotes', 'annot')
        manicure_dir = os.path.join(self.output_dir, 'eukaryotes', 'manicure')
        os.makedirs(annot_dir, exist_ok=True)
        os.makedirs(manicure_dir, exist_ok=True)

        FN = sample_id if sample_id else os.path.basename(genome_file).replace(self.extension, '')

        faa_file = os.path.join(annot_dir, f"{FN}.fas")
        ffn_file = os.path.join(annot_dir, f"{FN}.codon.fas")

        if (os.path.exists(faa_file) and os.path.getsize(faa_file) > 0 and
            os.path.exists(ffn_file) and os.path.getsize(ffn_file) > 0 and
            not self.force):
            logger.info(f"Skipping ORF calling for {genome_file} as output already exists.")
        else:
            log_file = os.path.join(self.output_dir, 'eukaryotes', f"{FN}_metaeuk.log")
            cmd = ['metaeuk', 'easy-predict', genome_file, self.eukdb, os.path.join(annot_dir, f"{FN}"), os.path.join(annot_dir, f"{FN}"), '--threads', str(self.threads)]
            with open(log_file, 'w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)

        manicure_faa = os.path.join(manicure_dir, f"{FN}.faa")
        manicure_ffn = os.path.join(manicure_dir, f"{FN}.ffn")
        try:
            if os.path.exists(manicure_faa):
                os.remove(manicure_faa)
            if os.path.exists(manicure_ffn):
                os.remove(manicure_ffn)
            os.symlink(faa_file, manicure_faa)
            os.symlink(ffn_file, manicure_ffn)
        except OSError:
            import shutil
            shutil.copy2(faa_file, manicure_faa)
            shutil.copy2(ffn_file, manicure_ffn)
        return True

    def call_metagenome_orfs(self, genome_file, sample_id=None):
        annot_dir = os.path.join(self.output_dir, 'metagenomes', 'annot')
        manicure_dir = os.path.join(self.output_dir, 'metagenomes', 'manicure')
        os.makedirs(annot_dir, exist_ok=True)
        os.makedirs(manicure_dir, exist_ok=True)

        FN = sample_id if sample_id else os.path.basename(genome_file).replace(self.extension, '')

        faa_file = os.path.join(annot_dir, f"{FN}.faa")
        ffn_file = os.path.join(annot_dir, f"{FN}.ffn")
        gff_file = os.path.join(annot_dir, f"{FN}.gff")

        if (os.path.exists(faa_file) and os.path.getsize(faa_file) > 0 and
            os.path.exists(ffn_file) and os.path.getsize(ffn_file) > 0 and
            os.path.exists(gff_file) and os.path.getsize(gff_file) > 0 and
            not self.force):
            logger.info(f"Skipping ORF calling for {genome_file} as output already exists.")
        else:
            log_file = os.path.join(self.output_dir, 'metagenomes', f"{FN}_prodigal.log")
            cmd = ['prodigal', '-p', 'meta', '-q', '-i', genome_file, '-d', ffn_file, '-a', faa_file, '-o', gff_file, '-f', 'gbk']
            with open(log_file, 'w') as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)

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
        return True

def find_genome_files(mag_dir, extension, wildcard):
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
    def process_single_genome(genome_data):
        sample_id, genome_path, domain = genome_data
        domain = domain.lower()
        try:
            if domain == 'bacterial':
                return orf_caller.call_bacterial_orfs(genome_path, sample_id=sample_id)
            elif domain == 'viral':
                return orf_caller.call_viral_orfs(genome_path, sample_id=sample_id)
            elif domain == 'eukaryotic':
                return orf_caller.call_eukaryotic_orfs(genome_path, sample_id=sample_id)
            elif domain == 'metagenomic':
                return orf_caller.call_metagenome_orfs(genome_path, sample_id=sample_id)
            else:
                logger.error(f"Unknown domain '{domain}' for sample {sample_id}. Skipping.")
                return False
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

def main():
    parser = argparse.ArgumentParser(description='Call ORFs for bacterial, viral, eukaryotic, and/or metagenomic genomes.')
    parser.add_argument('--config', type=str, default=None, help='Config file with 2 or 3 columns: sample_id genome_path [domain]. Supports tab or whitespace delimiters.')
    parser.add_argument('-m', '--mag_dir', type=str, default=None, help='Path or glob to genome files (e.g. asm/*/bins).')
    parser.add_argument('-w', '--wildcard', type=str, default='', help='Pattern to match anywhere in genome file path (pipe-separated for multiple).')
    parser.add_argument('--domain', type=str, choices=['bacterial', 'viral', 'eukaryotic', 'metagenomic'], help='Domain type when using directory mode.')

    parser.add_argument('--output_directory', type=str, default='magus_output/orf_calling', help='Directory to store ORF outputs.')
    parser.add_argument('--max_workers', type=int, default=1, help='Parallel ORF calling jobs.')
    parser.add_argument('--threads', type=int, default=4, help='Threads for tools.')
    parser.add_argument('--extension', type=str, default='fa', help='Extension of genome files.')
    parser.add_argument('--force', action='store_true', help='Force rewriting of existing outputs.')
    parser.add_argument('--eukdb', type=str, default='data/uniref90', help='UniRef90 DB path for MetaEuk.')

    args = parser.parse_args()

    if not args.config and not args.mag_dir:
        parser.error("Either --config or --mag_dir must be provided.")
    if args.config and args.mag_dir:
        parser.error("Cannot use both --config and --mag_dir.")
    if args.mag_dir and not args.domain:
        parser.error("--domain must be specified when using --mag_dir.")

    os.makedirs(args.output_directory, exist_ok=True)
    orf_caller = ORFCaller(args.output_directory, args.extension, args.threads, args.force, args.eukdb)

    genomes_data = []
    if args.config:
        logger.info(f"Using config file: {args.config}")
        with open(args.config, 'r') as f:
            for line_number, raw_line in enumerate(f, start=1):
                line = raw_line.strip()
                if not line or line.startswith('#'):
                    continue

                # Prefer tab-delimited parsing, but gracefully support whitespace-delimited files.
                row = line.split('\t') if '\t' in line else line.split()
                if len(row) < 2:
                    logger.warning(f"Skipping malformed config line {line_number}: {raw_line.rstrip()}")
                    continue

                sample_id, genome_path = row[0], row[1]
                if len(row) >= 3 and row[2].strip():
                    domain = row[2].strip()
                elif args.domain:
                    domain = args.domain
                else:
                    raise ValueError(
                        f"Config line {line_number} has no domain column and --domain was not provided: {raw_line.rstrip()}"
                    )

                genomes_data.append((sample_id, genome_path, domain))
    else:
        logger.info(f"Using directory mode: {args.mag_dir}")
        wildcards = args.wildcard.split('|') if args.wildcard else ['']
        all_genome_files = []
        for wildcard in wildcards:
            try:
                genome_files = find_genome_files(args.mag_dir, args.extension, wildcard.strip())
                all_genome_files.extend(genome_files)
                logger.info(f"Found {len(genome_files)} files matching wildcard '{wildcard.strip()}'")
            except RuntimeError as e:
                if len(wildcards) == 1:
                    raise e
        seen = set()
        unique_genome_files = []
        for path in all_genome_files:
            if path not in seen:
                seen.add(path)
                unique_genome_files.append(path)
        if not unique_genome_files:
            raise RuntimeError("No matching genome files found.")
        for genome_path in unique_genome_files:
            sample_id = genome_path.name
            extension = f".{args.extension}" if not args.extension.startswith('.') else args.extension
            if sample_id.endswith(extension):
                sample_id = sample_id[:-len(extension)]
            genomes_data.append((sample_id, str(genome_path), args.domain))

    if not genomes_data:
        raise RuntimeError(
            "No genomes were loaded from the provided inputs. Check config formatting or input paths."
        )

    process_genomes(orf_caller, genomes_data, args.max_workers)
    logger.info("ORF calling (no annotation) completed.")

if __name__ == '__main__':
    main()

