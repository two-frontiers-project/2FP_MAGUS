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

        # Ensure manicured files exist (even if ORF calling was skipped previously)
        source_faa = os.path.join(annot_dir, f"{FN}.faa")
        source_ffn = os.path.join(annot_dir, f"{FN}.ffn")
        manicure_faa = os.path.join(manicure_dir, f"{FN}.faa")
        manicure_ffn = os.path.join(manicure_dir, f"{FN}.ffn")

        if os.path.exists(source_faa) and os.path.getsize(source_faa) > 0:
            if not os.path.exists(manicure_faa) or self.args.force:
                with open(source_faa, 'r') as infile, open(manicure_faa, 'w') as outfile:
                    for line in infile:
                        if line.startswith('>'):
                            outfile.write(line.replace('>',f'>{FN}-----').replace(' # ', '-----', 1).replace('*', 'X').replace(' # ', '+').replace(' # ', '-'))
                        else:
                            outfile.write(line)
        else:
            logger.warning(f"Expected protein file missing for manicure: {source_faa}")

        if os.path.exists(source_ffn) and os.path.getsize(source_ffn) > 0:
            if not os.path.exists(manicure_ffn) or self.args.force:
                with open(source_ffn, 'r') as infile, open(manicure_ffn, 'w') as outfile:
                    for line in infile:
                        if line.startswith('>'):
                            outfile.write(line.replace('>',f'>{FN}-----').replace('+', '-----+').replace('*', 'X').replace(' # ', '+').replace(' # ', '-'))
                        else:
                            outfile.write(line)
        else:
            logger.warning(f"Expected nucleotide file missing for manicure: {source_ffn}")

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
        
        # HMM annotation (put results in annot directory for summary)
        if self.args.hmmfile:
            self.run_hmm_annotation(faa_file, annot_dir, FN, self.args.hmmfile)
        
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
        
        # HMM annotation (put results in annot directory for summary)
        if self.args.hmmfile:
            self.run_hmm_annotation(faa_file, annot_dir, FN, self.args.hmmfile)
        
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
                # Include ID in the summary so we can join HMM hits reliably
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
        
        # Now add HMM results if available
        if hmmfile:
            logger.info(f"Adding full HMM results to {subdir} summary")

            # Helper: parse a HMMER tblout data line into a dict of all columns
            def parse_hmm_tblout_line(data_line: str) -> dict:
                # According to HMMER, first 18 tokens are fixed-width fields; the rest is description
                tokens = data_line.rstrip('\n').split()
                if len(tokens) < 18:
                    return {}
                # Some HMMER versions include accession columns (target acc, query acc) right after names
                # We normalize by checking token count
                # Layout with accessions: target, tacc, query, qacc, fs_e, fs_s, fs_b, bd_e, bd_s, bd_b, exp, reg, clu, ov, env, dom, rep, inc, desc...
                # Layout without accessions: target, query, fs_e, fs_s, fs_b, bd_e, bd_s, bd_b, exp, reg, clu, ov, env, dom, rep, inc, desc...
                has_accessions = False
                # Heuristic: if tokens[1] != tokens[3] and tokens[1] contains characters like '.' or letters but not numeric-only lengths
                if len(tokens) >= 20:  # generous threshold; tblout with accessions typically has >= 20 tokens before desc
                    has_accessions = True
                idx = 0
                result = {}
                result['target_name'] = tokens[idx]; idx += 1
                if has_accessions:
                    result['target_accession'] = tokens[idx]; idx += 1
                result['query_name'] = tokens[idx]; idx += 1
                if has_accessions:
                    result['query_accession'] = tokens[idx]; idx += 1
                # full sequence stats
                result['full_evalue'] = tokens[idx]; idx += 1
                result['full_score'] = tokens[idx]; idx += 1
                result['full_bias'] = tokens[idx]; idx += 1
                # best 1 domain stats
                result['dom_evalue'] = tokens[idx]; idx += 1
                result['dom_score'] = tokens[idx]; idx += 1
                result['dom_bias'] = tokens[idx]; idx += 1
                # domain number estimation
                result['exp'] = tokens[idx]; idx += 1
                result['reg'] = tokens[idx]; idx += 1
                result['clu'] = tokens[idx]; idx += 1
                result['ov'] = tokens[idx]; idx += 1
                result['env'] = tokens[idx]; idx += 1
                result['dom'] = tokens[idx]; idx += 1
                result['rep'] = tokens[idx]; idx += 1
                result['inc'] = tokens[idx]; idx += 1
                # Remaining tokens (if any) form the description
                desc = ' '.join(tokens[idx:]) if idx < len(tokens) else ''
                result['description'] = desc
                # Try to extract ID=... from description for robust joining
                m = re.search(r'ID=([^;\s]+)', desc)
                if m:
                    result['gene_id'] = m.group(1)
                else:
                    # fall back to try extracting from target_name if it contains metadata
                    m2 = re.search(r'ID=([^;\s]+)', result['target_name'])
                    result['gene_id'] = m2.group(1) if m2 else None
                return result

            # Collect HMM results across samples, keyed by gene ID when available, else by target name
            hmm_by_key = {}
            for file in os.listdir(annot_dir):
                if file.endswith('.hmm.tsv'):
                    hmm_path = os.path.join(annot_dir, file)
                    with open(hmm_path, 'r') as fin:
                        for raw in fin:
                            if not raw or raw.startswith('#'):
                                continue
                            rec = parse_hmm_tblout_line(raw)
                            if not rec:
                                continue
                            key = rec.get('gene_id') or rec['target_name']
                            if key not in hmm_by_key:
                                hmm_by_key[key] = []
                            hmm_by_key[key].append(rec)

            if hmm_by_key:
                # Read existing summary and add all HMM columns
                temp_summary = summary_file + '.tmp'
                with open(summary_file, 'r') as infile, open(temp_summary, 'w') as outfile:
                    header = infile.readline().rstrip('\n')
                    hmm_cols = [
                        'hmm_target_name',
                        'hmm_target_accession',
                        'hmm_query_name',
                        'hmm_query_accession',
                        'hmm_full_evalue', 'hmm_full_score', 'hmm_full_bias',
                        'hmm_dom_evalue', 'hmm_dom_score', 'hmm_dom_bias',
                        'hmm_exp', 'hmm_reg', 'hmm_clu', 'hmm_ov', 'hmm_env', 'hmm_dom', 'hmm_rep', 'hmm_inc',
                        'hmm_description'
                    ]
                    outfile.write(f"{header}\t" + "\t".join(hmm_cols) + "\n")

                    # Determine how to fetch the join key for each row
                    # For non-euks we printed ID at column index 5
                    # For euks we don't know; we try to detect 'ID' column by name, else fall back to contig-like column
                    header_fields = header.split('\t')
                    try:
                        id_index = header_fields.index('ID')
                    except ValueError:
                        id_index = None
                    # Common fallbacks
                    contig_index = None
                    for cand in ['contig_id', 'contig', 'sequence', 'seq_id']:
                        if cand in header_fields:
                            contig_index = header_fields.index(cand)
                            break

                    for row in infile:
                        row = row.rstrip('\n')
                        if not row:
                            continue
                        fields = row.split('\t')
                        join_key = None
                        if id_index is not None and id_index < len(fields):
                            join_key = fields[id_index]
                        elif contig_index is not None and contig_index < len(fields):
                            join_key = fields[contig_index]
                        else:
                            # As a last resort, use the second column
                            join_key = fields[1] if len(fields) > 1 else fields[0]

                        hits = hmm_by_key.get(join_key, [])
                        # Aggregate each column across hits using ';' join
                        def agg(col):
                            vals = []
                            for h in hits:
                                v = h.get(col)
                                if v is None:
                                    # ensure positions exist for accessions when missing
                                    if col in ['hmm_target_accession', 'hmm_query_accession']:
                                        vals.append('')
                                    else:
                                        continue
                                else:
                                    vals.append(str(v))
                            return ';'.join(vals)

                        # Prepare output columns in the same order as hmm_cols
                        col_map = {
                            'hmm_target_name': agg('target_name'),
                            'hmm_target_accession': agg('target_accession'),
                            'hmm_query_name': agg('query_name'),
                            'hmm_query_accession': agg('query_accession'),
                            'hmm_full_evalue': agg('full_evalue'),
                            'hmm_full_score': agg('full_score'),
                            'hmm_full_bias': agg('full_bias'),
                            'hmm_dom_evalue': agg('dom_evalue'),
                            'hmm_dom_score': agg('dom_score'),
                            'hmm_dom_bias': agg('dom_bias'),
                            'hmm_exp': agg('exp'),
                            'hmm_reg': agg('reg'),
                            'hmm_clu': agg('clu'),
                            'hmm_ov': agg('ov'),
                            'hmm_env': agg('env'),
                            'hmm_dom': agg('dom'),
                            'hmm_rep': agg('rep'),
                            'hmm_inc': agg('inc'),
                            'hmm_description': agg('description'),
                        }
                        outfile.write(row + '\t' + '\t'.join(col_map[c] for c in hmm_cols) + '\n')

                os.replace(temp_summary, summary_file)
        
        logger.info(f"Created comprehensive {subdir} summary: {summary_file}")

def summarize_annotations(output_dir):
    """Create manicured files and final summary with annotations."""
    for subdir in ['bacteria', 'viruses', 'eukaryotes', 'metagenomes']:
        annot_dir = os.path.join(output_dir, subdir, 'annot')
        manicure_dir = os.path.join(output_dir, subdir, 'manicure')
        os.makedirs(manicure_dir, exist_ok=True)
        
        if not os.path.exists(annot_dir):
            continue
        
        # Special handling for eukaryotes: concatenate MetaEuk headers maps with a sample column
        if subdir == 'eukaryotes':
            summary_file = os.path.join(output_dir, f'{subdir}_final_summary.tsv')
            header_written = False
            with open(summary_file, 'w') as summary:
                for file in os.listdir(annot_dir):
                    if file.endswith('.headersMap.tsv'):
                        sample_id = file.replace('.headersMap.tsv', '')
                        file_path = os.path.join(annot_dir, file)
                        with open(file_path, 'r') as infile:
                            for i, line in enumerate(infile):
                                line = line.rstrip('\n')
                                if not line.strip():
                                    continue
                                
                                # Write header only once (first line of first file)
                                if i == 0 and not header_written:
                                    summary.write(f'sample_id\t{line}\n')
                                    header_written = True
                                elif i > 0 or header_written:  # Skip header lines from subsequent files
                                    summary.write(f'{sample_id}\t{line}\n')
            logger.info(f"Created {subdir} calls summary: {summary_file}")
            
            # Also create manicured FASTA files for eukaryotes
            for file in os.listdir(annot_dir):
                if file.endswith('.fas'):
                    sample_id = file.replace('.fas', '')
                    faa_file = os.path.join(annot_dir, file)
                    manicure_file = os.path.join(manicure_dir, f"{sample_id}.faa")
                    
                    # Manicure the protein file
                    with open(faa_file, 'r') as infile, open(manicure_file, 'w') as outfile:
                        for line in infile:
                            if line.startswith('>'):
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
                                outfile.write(line)
                    
                    # Also manicure the nucleotide file if it exists
                    ffn_file = faa_file.replace('.fas', '.codon.fas')
                    if os.path.exists(ffn_file):
                        manicure_ffn = os.path.join(manicure_dir, f"{sample_id}.ffn")
                        with open(ffn_file, 'r') as infile, open(manicure_ffn, 'w') as outfile:
                            for line in infile:
                                if line.startswith('>'):
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
                                    outfile.write(line)
            
            # Skip the generic manicure/summarization logic for eukaryotes
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
                                    outfile.write(line.replace('>',f'>{sample_id}-----').replace('+', '-----+').replace('*', 'X').replace(' # ', '+').replace(' # ', '-'))
                                else:
                                    outfile.write(line)
        
        # Create final summary for this domain (non-eukaryotes)
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
    parser.add_argument('--restart', type=str, choices=['create-summary'], 
                        help='Restart from specific stage: create-summary')
    
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
