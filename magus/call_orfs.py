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
                                 orfcalling_evalue_cutoff: Optional[float] = None,
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
                # Eukaryotes: parse MetaEuk headersMap.tsv lines and normalize columns.
                # Each data line's last column is a pipe-delimited block:
                # UniRefAcc|sample|strand|score|evalue|num_exons|start|end|exon_info
                header_written = False
                for file in os.listdir(annot_dir):
                    if file.endswith('.headersMap.tsv'):
                        sample_id = file.replace('.headersMap.tsv', '')
                        file_path = os.path.join(annot_dir, file)

                        with open(file_path, 'r') as infile:
                            for raw in infile:
                                raw = raw.rstrip('\n')
                                if not raw.strip():
                                    continue
                                # Prepare header once
                                if not header_written:
                                    summary.write('sample_id\taccession\tsample\tstrand\tscore\tevalue\tnum_exons\tstart\tend\texon_info\n')
                                    header_written = True
                                # Extract last column and parse
                                parts = raw.split('\t')
                                pipe_block = parts[-1] if parts else ''
                                pipe_fields = pipe_block.split('|') if pipe_block else []
                                # Pad to at least 9 fields
                                while len(pipe_fields) < 9:
                                    pipe_fields.append('')
                                accession, sample_name, strand, score, evalue, num_exons, start, end, exon_info = pipe_fields[:9]
                                summary.write('\t'.join([
                                    sample_id,
                                    accession,
                                    sample_name,
                                    strand,
                                    score,
                                    evalue,
                                    num_exons,
                                    start,
                                    end,
                                    exon_info
                                ]) + '\n')
            else:
                # Bacteria/Viruses/Metagenomes: parse Prodigal output
                # Include ID so we can correctly join HMM hits per ORF
                summary.write('sample_id\tcontig_id\tstart\tend\tstrand\tID\tpartial\tstart_type\trbs_motif\trbs_spacer\tgc_cont\n')
                
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
                                        gene_id = 'NA'
                                        partial = '00'
                                        start_type = 'ATG'
                                        rbs_motif = 'None'
                                        rbs_spacer = 'None'
                                        gc_cont = '0.500'
                                        
                                        if 'ID=' in metadata:
                                            gene_id = metadata.split('ID=')[1].split(';')[0]
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
                                        
                                        summary.write(f'{sample_id}\t{contig_id}\t{start}\t{end}\t{strand}\t{gene_id}\t{partial}\t{start_type}\t{rbs_motif}\t{rbs_spacer}\t{gc_cont}\n')
        
        # Now add HMM results from existing files (do not depend on --hmmfile being provided here)
        logger.info(f"Adding HMM results to {subdir} summary (if present)")

        # Parse HMM tblout and build records
        def parse_tblout_row(row: str):
                if not row or row.startswith('#'):
                    return None
                t = row.rstrip('\n').split()
                if len(t) < 19:
                    return None
                rec = {
                    'target_name': t[0],
                    'target_accession': t[1],
                    'query_name': t[2],
                    'query_accession': t[3],
                    'full_evalue': t[4],
                    'full_score': t[5],
                    'full_bias': t[6],
                    'dom_evalue': t[7],
                    'dom_score': t[8],
                    'dom_bias': t[9],
                    'exp': t[10], 'reg': t[11], 'clu': t[12], 'ov': t[13], 'env': t[14], 'dom': t[15], 'rep': t[16], 'inc': t[17],
                    'description': ' '.join(t[18:])
                }
                # Extract possible gene ID from target_name (for non-euks)
                m = re.search(r'ID=([^;\s]+)', rec['target_name'])
                rec['gene_id'] = m.group(1) if m else None
                return rec

        # Collect HMM records grouped by join key
        hmm_by_key = {}
        for file in os.listdir(annot_dir):
            if (suffix and file.endswith(f'.hmm.{suffix}.tsv')) or (not suffix and file.endswith('.hmm.tsv')):
                hmm_path = os.path.join(annot_dir, file)
                with open(hmm_path, 'r') as fin:
                    for raw in fin:
                        rec = parse_tblout_row(raw)
                        if not rec:
                            continue
                        # Choose join key
                        if subdir == 'eukaryotes':
                                # For MetaEuk targets, the accession is the first token before '|'
                                # Example target_name: UniRef90_ID|sample|strand|...
                                tn = rec.get('target_name', '')
                                acc_from_name = tn.split('|', 1)[0] if '|' in tn else tn
                                join_key = acc_from_name
                        else:
                            join_key = rec.get('gene_id') or rec.get('target_name')
                        if not join_key:
                            continue
                        hmm_by_key.setdefault(join_key, []).append(rec)

            if hmm_by_key:
                temp_summary = summary_file + '.tmp'
                with open(summary_file, 'r') as infile, open(temp_summary, 'w') as outfile:
                    original_header = infile.readline().rstrip('\n')
                    hmm_cols = [
                        'hmm_target_name','hmm_target_accession','hmm_query_name','hmm_query_accession',
                        'hmm_full_evalue','hmm_full_score','hmm_full_bias',
                        'hmm_dom_evalue','hmm_dom_score','hmm_dom_bias',
                        'hmm_exp','hmm_reg','hmm_clu','hmm_ov','hmm_env','hmm_dom','hmm_rep','hmm_inc','hmm_description'
                    ]
                    # Write a single header row only once (long format: one row per HMM hit)
                    outfile.write(original_header + '\t' + '\t'.join(hmm_cols) + '\n')

                    header_fields = original_header.split('\t')
                    if subdir == 'eukaryotes' and 'accession' in header_fields:
                        key_idx = header_fields.index('accession')
                    elif 'ID' in header_fields:
                        key_idx = header_fields.index('ID')
                    else:
                        key_idx = header_fields.index('contig_id') if 'contig_id' in header_fields else 1

                    # Determine euk evalue column index if present for ORF-calling evalue filtering
                    euk_evalue_idx = header_fields.index('evalue') if 'evalue' in header_fields else None

                    def parse_float_safe(value: str) -> Optional[float]:
                        try:
                            return float(value)
                        except Exception:
                            return None

                    for row in infile:
                        row = row.rstrip('\n')
                        if not row:
                            continue
                        cols = row.split('\t')
                        join_key = cols[key_idx] if key_idx < len(cols) else None
                        hits = hmm_by_key.get(join_key, []) if join_key else []
                        # Long format: emit one row per HMM hit; also emit a blank-annotation row if there are no hits
                        if hits:
                            for h in hits:
                                # Apply filters
                                # 1) ORF-calling evalue (euk only)
                                if subdir == 'eukaryotes' and orfcalling_evalue_cutoff is not None and euk_evalue_idx is not None:
                                    ev = parse_float_safe(cols[euk_evalue_idx]) if euk_evalue_idx < len(cols) else None
                                    if ev is not None and ev > orfcalling_evalue_cutoff:
                                        continue
                                # 2) HMM full-sequence evalue
                                if hmm_fullseq_evalue_cutoff is not None:
                                    fev = parse_float_safe(h.get('full_evalue', ''))
                                    if fev is not None and fev > hmm_fullseq_evalue_cutoff:
                                        continue
                                # 3) HMM best-domain evalue
                                if hmm_domain_evalue_cutoff is not None:
                                    dev = parse_float_safe(h.get('dom_evalue', ''))
                                    if dev is not None and dev > hmm_domain_evalue_cutoff:
                                        continue
                                values = [
                                    h.get('target_name',''), h.get('target_accession',''), h.get('query_name',''), h.get('query_accession',''),
                                    h.get('full_evalue',''), h.get('full_score',''), h.get('full_bias',''),
                                    h.get('dom_evalue',''), h.get('dom_score',''), h.get('dom_bias',''),
                                    h.get('exp',''), h.get('reg',''), h.get('clu',''), h.get('ov',''), h.get('env',''), h.get('dom',''), h.get('rep',''), h.get('inc',''), h.get('description','')
                                ]
                                outfile.write(row + '\t' + '\t'.join(values) + '\n')
                        else:
                            # No HMMs for this ORF: apply only ORF-calling evalue filter (if any); then keep with empty HMM columns
                            if subdir == 'eukaryotes' and orfcalling_evalue_cutoff is not None and euk_evalue_idx is not None:
                                ev = parse_float_safe(cols[euk_evalue_idx]) if euk_evalue_idx < len(cols) else None
                                if ev is not None and ev > orfcalling_evalue_cutoff:
                                    continue
                            outfile.write(row + '\t' + '\t'.join([''] * len(hmm_cols)) + '\n')

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
    parser.add_argument('--orfcalling-evalue', '--oce', type=float, default=None, dest='orfcalling_evalue', help='Filter out ORFs whose calling evalue exceeds this cutoff (applies to eukaryote MetaEuk evalue column).')
    parser.add_argument('--annotation-fullseq-evalue', '--afe', type=float, default=None, dest='annotation_fullseq_evalue', help='Filter out HMM hits with full-sequence E-value exceeding this cutoff.')
    parser.add_argument('--annotation-domain-evalue', '--ade', type=float, default=None, dest='annotation_domain_evalue', help='Filter out HMM hits with best-domain E-value exceeding this cutoff.')
    parser.add_argument('--suffix', type=str, default=None, help='Suffix to tag outputs (e.g., kegg). Summaries and HMM tblout files include this suffix.')
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
        create_comprehensive_summary(
            args.output_directory, args.hmmfile, args.suffix,
            orfcalling_evalue_cutoff=args.orfcalling_evalue,
            hmm_fullseq_evalue_cutoff=args.annotation_fullseq_evalue,
            hmm_domain_evalue_cutoff=args.annotation_domain_evalue,
        )
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
        
        # Create comprehensive summaries with new HMM results
        logger.info("Creating comprehensive summaries with updated HMM results")
        create_comprehensive_summary(
            args.output_directory, args.hmmfile, args.suffix,
            orfcalling_evalue_cutoff=args.orfcalling_evalue,
            hmm_fullseq_evalue_cutoff=args.annotation_fullseq_evalue,
            hmm_domain_evalue_cutoff=args.annotation_domain_evalue,
        )
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
    create_comprehensive_summary(
        args.output_directory, args.hmmfile, args.suffix,
        orfcalling_evalue_cutoff=args.orfcalling_evalue,
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
