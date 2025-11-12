#!/usr/bin/env python3

import os
import argparse
import logging
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import csv

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Use hmmsearch-g from PATH
HMMSEARCH_BIN = 'hmmsearch-g'

def list_targets(output_dir, domain):
    if domain == 'eukaryotes':
        annot_dir = os.path.join(output_dir, domain, 'annot')
        os.makedirs(annot_dir, exist_ok=True)
        files = [(f.replace('.fas', ''), os.path.join(annot_dir, f)) for f in os.listdir(annot_dir) if f.endswith('.fas') and not f.endswith('.codon.fas')]
        return annot_dir, files
    else:
        manicure_dir = os.path.join(output_dir, domain, 'manicure')
        annot_dir = os.path.join(output_dir, domain, 'annot')
        os.makedirs(annot_dir, exist_ok=True)
        files = [(f.replace('.faa', ''), os.path.join(manicure_dir, f)) for f in os.listdir(manicure_dir) if f.endswith('.faa')]
        return annot_dir, files

def get_input_files(sequence_dir=None, sequence_file=None, extension='faa'):
    if sequence_dir and sequence_file:
        raise ValueError("Cannot specify both --sequence-dir and --sequence-file")
    if not sequence_dir and not sequence_file:
        raise ValueError("Must specify either --sequence-dir or --sequence-file")
    
    if sequence_file:
        seq_path = Path(sequence_file)
        if not seq_path.is_file():
            raise ValueError(f"Sequence file does not exist: {sequence_file}")
        return [seq_path]
    elif sequence_dir:
        seq_dir = Path(sequence_dir)
        if not seq_dir.is_dir():
            raise ValueError(f"Sequence directory does not exist: {sequence_dir}")
        ext = f".{extension}" if not extension.startswith('.') else extension
        files = list(seq_dir.glob(f"*{ext}"))
        if not files:
            raise ValueError(f"No {ext} files found in {sequence_dir}")
        return files

def split_fasta_file(input_file, output_dir, sequences_per_file=100000):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    split_prefix = output_dir / f"{input_file.stem}_split"
    awk_script = f'''BEGIN {{ file_num = 0; seq_count = 0 }}
/^>/ {{
    if (seq_count > 0 && seq_count % {sequences_per_file} == 0) {{
        close(outfile)
        file_num++
    }}
    seq_count++
    outfile = "{split_prefix}_" file_num ".faa"
    print > outfile
    next
}}
{{ print > outfile }}'''
    with open(output_dir / "split_script.awk", 'w') as f:
        f.write(awk_script)
    cmd = ['awk', '-f', str(output_dir / "split_script.awk"), str(input_file)]
    subprocess.run(cmd, check=True)
    (output_dir / "split_script.awk").unlink()
    split_files = sorted(split_prefix.parent.glob(f"{split_prefix.name}_*.faa"))
    return split_files

def run_hmmsearch(sample_id, protein_file, out_dir, hmm_db, threads, suffix=None, Z=None, use_ga=True, nobias=True):
    os.makedirs(out_dir, exist_ok=True)
    tblout = os.path.join(out_dir, f"{sample_id}.hmm.{suffix}.tsv" if suffix else f"{sample_id}.hmm.tsv")
    log_file = os.path.join(out_dir, f"{sample_id}_hmmsearch.log" if not suffix else f"{sample_id}_hmmsearch.{suffix}.log")
    cmd = [HMMSEARCH_BIN, '--tblout', tblout, '--notextw', '--noali', '--cpu', str(threads)]
    if use_ga:
        cmd.append('--cut_ga')
    if nobias:
        cmd.append('--nobias')
    if Z is not None:
        cmd.extend(['-Z', str(Z)])
    cmd.extend([hmm_db, protein_file])
    with open(log_file, 'w') as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)
    return tblout

def parse_tblout_lines(tblout_path, evalue_full_cutoff=None, evalue_dom_cutoff=None):
    rows = {}
    if not os.path.exists(tblout_path):
        return rows
    with open(tblout_path, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) < 18:
                continue
            # Align parsing with call_orfs.parse_hmm_tblout_to_csv
            target_name = parts[0]
            query_name = parts[2] if len(parts) > 2 else ''
            try:
                full_evalue = float(parts[4]) if len(parts) > 4 else 0.0
                full_score = float(parts[5]) if len(parts) > 5 else 0.0
                dom_evalue = float(parts[7]) if len(parts) > 7 else 0.0
                dom_score = float(parts[8]) if len(parts) > 8 else 0.0
            except ValueError:
                continue
            if evalue_full_cutoff is not None and full_evalue > evalue_full_cutoff:
                continue
            if evalue_dom_cutoff is not None and dom_evalue > evalue_dom_cutoff:
                continue
            rows[target_name] = {
                'query_name': query_name,
                'full_evalue': str(full_evalue),
                'full_score': str(full_score),
                'dom_evalue': str(dom_evalue),
                'dom_score': str(dom_score)
            }
    return rows

def parse_tblout_to_scour(tblout_path, scour_path):
    """
    Parse hmmsearch-g tblout format into 5-column scour format:
    gene, match_raw_id, evalue, seq_score, domain_score
    Based on: grep -v '^#' $f | sed 's/ \+/\t/g' | awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\n",$1,$4,$5,$6,$9}'
    After sed converts spaces to tabs, columns are: 1=target, 4=query_accession, 5=evalue, 6=score, 9=domain_score
    """
    os.makedirs(os.path.dirname(scour_path), exist_ok=True)
    with open(tblout_path, 'r') as infile, open(scour_path, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) < 9:
                continue
            target_name = parts[0]
            query_accession = parts[3] if len(parts) > 3 else ''
            full_evalue = parts[4] if len(parts) > 4 else '0.0'
            full_score = parts[5] if len(parts) > 5 else '0.0'
            dom_score = parts[8] if len(parts) > 8 else '0.0'
            outfile.write(f"{target_name}\t{query_accession}\t{full_evalue}\t{full_score}\t{dom_score}\n")

def parse_scour_file(scour_path):
    """
    Parse 5-column TSVs produced by the provided shell script:
    columns: gene, match_raw_id, evalue, seq_score, domain_score
    """
    rows = {}
    if not os.path.exists(scour_path):
        return rows
    with open(scour_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 5:
                continue
            target_name = parts[0]
            query_name = parts[1]
            rows[target_name] = {
                'query_name': query_name,
                'full_evalue': parts[2],
                'full_score': parts[3],
                'dom_score': parts[4]
            }
    return rows

def load_tsv_mapping(tsv_path, key_col=None):
    """
    Load TSV mapping file (database lookup) into a dictionary.
    First column (or specified key_col) is the key (HMM ID), rest are metadata.
    Handles '#' prefix in column names like protparse2.R.
    """
    mapping = {}
    if not os.path.exists(tsv_path):
        return mapping
    with open(tsv_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames:
            return mapping
        
        if key_col:
            key_col_name = key_col
        else:
            key_col_name = reader.fieldnames[0]
            if key_col_name.startswith('#'):
                key_col_name = key_col_name[1:]
        
        for row in reader:
            hmm_id = row.get(key_col_name) or row.get(f"#{key_col_name}") or row[reader.fieldnames[0]]
            if hmm_id:
                mapping[hmm_id] = row
    return mapping

def merge_scours_with_mapping(scour_dir, mapping_tsv, output_file, file_type='pgap'):
    """
    Merge scours files with TSV mapping file and apply cutoffs.
    Based on protparse2.R logic.
    For PGAP: uses first column (NCBI ID) as key
    For Pfam: uses 'source_identifier' or second column as key
    """
    if file_type == 'pgap':
        mapping = load_tsv_mapping(mapping_tsv)
        key_col = None
    else:
        mapping = load_tsv_mapping(mapping_tsv, key_col='source_identifier')
        if not mapping:
            mapping = load_tsv_mapping(mapping_tsv)
        key_col = 'source_identifier'
    
    if not mapping:
        logger.warning(f"Mapping file {mapping_tsv} is empty or not found")
        return
    
    pattern = "*.search" if file_type == 'pgap' else "*.pfam"
    scour_files = list(Path(scour_dir).glob(pattern))
    
    if not scour_files:
        logger.warning(f"No {pattern} files found in {scour_dir}")
        return
    
    with open(output_file, 'w') as outfile:
        header = ['sample_id', 'gene', 'match', 'evalue', 'seq_score', 'domain_score',
                  'label', 'gene_symbol', 'product_name', 'ec_numbers', 'go_terms',
                  'for_AMRFinder', 'comment', 'annotation_source']
        outfile.write('\t'.join(header) + '\n')
        
        for scour_file in scour_files:
            sample_id = scour_file.stem.replace('.search', '').replace('.pfam', '')
            scours = parse_scour_file(str(scour_file))
            
            for gene, hit in scours.items():
                match_id = hit['query_name']
                if match_id not in mapping:
                    continue
                
                meta = mapping[match_id]
                try:
                    seq_score = float(hit['full_score'])
                    domain_score = float(hit['dom_score'])
                except ValueError:
                    continue
                
                # Get cutoffs from mapping file
                seq_cutoff_str = meta.get('sequence_cutoff', '')
                dom_cutoff_str = meta.get('domain_cutoff', '')
                
                # Filter out rows where cutoffs are missing (NA or empty) - protparse2.R line 411
                # merged_dt = merged_dt[!is.na(sequence_cutoff) & !is.na(domain_cutoff) & !is.na(effective_id_for_aggregation)]
                if not seq_cutoff_str or not dom_cutoff_str or seq_cutoff_str.strip() == '' or dom_cutoff_str.strip() == '':
                    continue
                
                try:
                    seq_cutoff = float(seq_cutoff_str)
                    dom_cutoff = float(dom_cutoff_str)
                except (ValueError, TypeError):
                    continue
                
                # Apply thresholds: seq_score >= sequence_cutoff AND domain_score >= domain_cutoff
                # protparse2.R line 415: filtered_hits_dt = merged_dt[seq_score >= sequence_cutoff & domain_score >= domain_cutoff]
                if seq_score < seq_cutoff or domain_score < dom_cutoff:
                    continue
                
                row = [
                    sample_id, gene, match_id,
                    hit['full_evalue'], hit['full_score'], hit['dom_score'],
                    meta.get('label', ''), meta.get('gene_symbol', ''),
                    meta.get('product_name', ''), meta.get('ec_numbers', ''),
                    meta.get('go_terms', ''), meta.get('for_AMRFinder', ''),
                    meta.get('comment', ''), 'PGAP' if file_type == 'pgap' else 'Pfam'
                ]
                outfile.write('\t'.join(str(x) for x in row) + '\n')
    
    logger.info(f"Wrote merged annotations to {output_file}")

def merge_pgap_and_pfam(pgap_file, pfam_file, output_file):
    """
    Merge PGAP and Pfam annotations into a single file with deduplication.
    Based on protparse2.R logic: remove Pfam entries that match PGAP on (gene, label) per sample.
    """
    pgap_keys = set()
    pgap_entries = []
    pfam_entries = []
    
    # Read PGAP annotations and build key set
    if os.path.exists(pgap_file):
        with open(pgap_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                sample_id = row.get('sample_id', '')
                gene = row.get('gene', '')
                label = row.get('label', '')
                key = (sample_id, gene, label)
                pgap_keys.add(key)
                pgap_entries.append(row)
    
    # Read Pfam annotations
    if os.path.exists(pfam_file):
        with open(pfam_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                pfam_entries.append(row)
    
    # Deduplicate Pfam: remove entries that match PGAP on (sample_id, gene, label)
    pfam_dedup = []
    for pfam_row in pfam_entries:
        sample_id = pfam_row.get('sample_id', '')
        gene = pfam_row.get('gene', '')
        label = pfam_row.get('label', '')
        key = (sample_id, gene, label)
        if key not in pgap_keys:
            pfam_dedup.append(pfam_row)
    
    # Combine PGAP + deduplicated Pfam
    all_entries = pgap_entries + pfam_dedup
    
    # Write merged file
    if all_entries:
        header = ['sample_id', 'gene', 'match', 'evalue', 'seq_score', 'domain_score',
                  'label', 'gene_symbol', 'product_name', 'ec_numbers', 'go_terms',
                  'for_AMRFinder', 'comment', 'annotation_source']
        with open(output_file, 'w') as outfile:
            outfile.write('\t'.join(header) + '\n')
            for row in all_entries:
                outfile.write('\t'.join(row.get(col, '') for col in header) + '\n')
        logger.info(f"Wrote merged PGAP+Pfam annotations to {output_file}")
    else:
        logger.warning("No annotations to merge")

def find_scour_for_sample(scour_dir, sample_id, suffix):
    """
    Determine expected scour filename based on suffix:
    - pfam -> {sample_id}.pfam
    - pgap -> {sample_id}.search
    """
    if suffix == 'pfam':
        candidate = os.path.join(scour_dir, f"{sample_id}.pfam")
        return candidate if os.path.exists(candidate) else None
    if suffix == 'pgap':
        candidate = os.path.join(scour_dir, f"{sample_id}.search")
        return candidate if os.path.exists(candidate) else None
    # Fallback: try {sample_id}.{suffix}
    candidate = os.path.join(scour_dir, f"{sample_id}.{suffix}")
    return candidate if os.path.exists(candidate) else None

def write_annotation_summary_bvm_from_tblouts(output_dir, domain, suffix, evalue_full_cutoff=None, evalue_dom_cutoff=None):
    annot_dir = os.path.join(output_dir, domain, 'annot')
    manicure_dir = os.path.join(output_dir, domain, 'manicure')
    summary_file = os.path.join(output_dir, f"{domain}_annotations_{suffix}.tsv" if suffix else f"{domain}_annotations.tsv")
    with open(summary_file, 'w') as f:
        f.write('\t'.join(['sample_id', 'target_name', 'query_name', 'full_evalue', 'full_score', 'dom_evalue', 'dom_score']) + '\n')
    files = [x for x in os.listdir(manicure_dir) if x.endswith('.faa')]
    for file in files:
        sample_id = file.replace('.faa', '')
        hmm_tbl = os.path.join(annot_dir, f"{sample_id}.hmm.{suffix}.tsv" if suffix else f"{sample_id}.hmm.tsv")
        hmm_rows = parse_tblout_lines(hmm_tbl, evalue_full_cutoff=evalue_full_cutoff, evalue_dom_cutoff=evalue_dom_cutoff)
        for target_name, hmm in hmm_rows.items():
            row_data = [
                sample_id, target_name,
                hmm.get('query_name', ''),
                hmm.get('full_evalue', ''), hmm.get('full_score', ''),
                hmm.get('dom_evalue', ''), hmm.get('dom_score', '')
            ]
            with open(summary_file, 'a') as out_f:
                out_f.write('\t'.join(str(x) for x in row_data) + '\n')
    logger.info(f"Wrote annotation summary: {summary_file}")

def write_annotation_summary_euk_from_tblouts(output_dir, suffix, evalue_full_cutoff=None, evalue_dom_cutoff=None):
    annot_dir = os.path.join(output_dir, 'eukaryotes', 'annot')
    summary_file = os.path.join(output_dir, f"eukaryotes_annotations_{suffix}.tsv" if suffix else "eukaryotes_annotations.tsv")
    with open(summary_file, 'w') as f:
        f.write('\t'.join(['sample_id', 'target_name', 'query_name', 'full_evalue', 'full_score', 'dom_evalue', 'dom_score']) + '\n')
    for file in os.listdir(annot_dir):
        if not (file.endswith('.fas') and not file.endswith('.codon.fas')):
            continue
        sample_id = file.replace('.fas', '')
        hmm_tbl = os.path.join(annot_dir, f"{sample_id}.hmm.{suffix}.tsv" if suffix else f"{sample_id}.hmm.tsv")
        hmm_rows = parse_tblout_lines(hmm_tbl, evalue_full_cutoff=evalue_full_cutoff, evalue_dom_cutoff=evalue_dom_cutoff)
        for target_name, hmm in hmm_rows.items():
            row = [sample_id, target_name, hmm.get('query_name', ''), hmm.get('full_evalue', ''), hmm.get('full_score', ''), hmm.get('dom_evalue', ''), hmm.get('dom_score', '')]
            with open(summary_file, 'a') as out_f:
                out_f.write('\t'.join(str(x) for x in row) + '\n')
    logger.info(f"Wrote annotation summary: {summary_file}")

def write_annotation_summary_bvm_from_scours(output_dir, domain, suffix, scour_dir):
    manicure_dir = os.path.join(output_dir, domain, 'manicure')
    summary_file = os.path.join(output_dir, f"{domain}_annotations_{suffix}.tsv")
    with open(summary_file, 'w') as f:
        f.write('\t'.join(['sample_id', 'target_name', 'query_name', 'full_evalue', 'full_score', 'dom_evalue', 'dom_score']) + '\n')
    files = [x for x in os.listdir(manicure_dir) if x.endswith('.faa')]
    for file in files:
        sample_id = file.replace('.faa', '')
        scour_path = find_scour_for_sample(scour_dir, sample_id, suffix)
        if not scour_path:
            continue
        rows = parse_scour_file(scour_path)
        for target_name, hmm in rows.items():
            row_data = [
                sample_id, target_name,
                hmm.get('query_name', ''),
                hmm.get('full_evalue', ''), hmm.get('full_score', ''),
                hmm.get('dom_evalue', ''), hmm.get('dom_score', '')
            ]
            with open(summary_file, 'a') as out_f:
                out_f.write('\t'.join(str(x) for x in row_data) + '\n')
    logger.info(f"Wrote annotation summary: {summary_file}")

def write_annotation_summary_euk_from_scours(output_dir, suffix, scour_dir):
    annot_dir = os.path.join(output_dir, 'eukaryotes', 'annot')
    summary_file = os.path.join(output_dir, f"eukaryotes_annotations_{suffix}.tsv")
    with open(summary_file, 'w') as f:
        f.write('\t'.join(['sample_id', 'target_name', 'query_name', 'full_evalue', 'full_score', 'dom_evalue', 'dom_score']) + '\n')
    for file in os.listdir(annot_dir):
        if not (file.endswith('.fas') and not file.endswith('.codon.fas')):
            continue
        sample_id = file.replace('.fas', '')
        scour_path = find_scour_for_sample(scour_dir, sample_id, suffix)
        if not scour_path:
            continue
        rows = parse_scour_file(scour_path)
        for target_name, hmm in rows.items():
            row = [sample_id, target_name, hmm.get('query_name', ''), hmm.get('full_evalue', ''), hmm.get('full_score', ''), hmm.get('dom_evalue', ''), hmm.get('dom_score', '')]
            with open(summary_file, 'a') as out_f:
                out_f.write('\t'.join(str(x) for x in row) + '\n')
    logger.info(f"Wrote annotation summary: {summary_file}")

def main():
    parser = argparse.ArgumentParser(description='Annotate proteins with hmmsearch-g and merge results.')
    parser.add_argument('--output_directory', type=str, default='magus_output/orf_calling', help='Root ORF output directory produced by call_orfs.')
    parser.add_argument('--faa_dir', type=str, default=None, help='Optional: directory containing .faa files to annotate (overrides discovery).')
    parser.add_argument('--sequence-dir', type=str, default=None, dest='sequence_dir', help='Directory containing FASTA files to annotate.')
    parser.add_argument('--sequence-file', type=str, default=None, dest='sequence_file', help='Single FASTA file to annotate.')
    parser.add_argument('--split-file-size', type=int, default=100000, help='Number of sequences per split file when splitting a single FASTA (default: 100000).')
    parser.add_argument('-x', '--extension', type=str, default='faa', help='File extension when using --sequence-dir (default: faa).')
    parser.add_argument('--domains', type=str, default='bacteria,viruses,metagenomes,eukaryotes', help='Comma-separated domains to process.')
    parser.add_argument('--threads', type=int, default=8, help='Threads per hmmsearch-g job.')
    parser.add_argument('--max_workers', type=int, default=4, help='Parallel samples.')
    
    parser.add_argument('--pfam_tsv', type=str, default=None, help='Path to Pfam mapping TSV file (database lookup).')
    parser.add_argument('--pgap_tsv', type=str, default=None, help='Path to PGAP mapping TSV file (database lookup).')
    parser.add_argument('--pfam_db', type=str, default=None, help='Path to Pfam-A HMM database (.hmm file).')
    parser.add_argument('--pgap_db', type=str, default=None, help='Path to PGAP HMM database (.hmm file).')
    parser.add_argument('--Z_pfam', type=int, default=25545, help='Pfam database size for -Z.')
    parser.add_argument('--Z_pgap', type=int, default=18057, help='PGAP database size for -Z.')
    
    parser.add_argument('--hmmdb', type=str, default=None, help='Path to HMM database (user-specified mode).')
    parser.add_argument('--suffix', type=str, default='custom', help='Suffix tag for outputs.')
    parser.add_argument('--evalue_full', type=float, default=None, help='Full-sequence E-value cutoff.')
    parser.add_argument('--evalue_dom', type=float, default=None, help='Best-domain E-value cutoff.')
    parser.add_argument('--no_cut_ga', action='store_true', help='Do not apply --cut_ga.')
    parser.add_argument('--Z', type=int, default=None, help='Database size for -Z.')

    args = parser.parse_args()

    domains = [d.strip() for d in args.domains.split(',') if d.strip()]
    out_root = args.output_directory

    if args.pfam_db or args.pgap_db or args.pfam_tsv or args.pgap_tsv:
        if not all([args.pfam_db, args.pgap_db, args.pfam_tsv, args.pgap_tsv]):
            parser.error("Default mode requires all four arguments: --pfam_db, --pgap_db, --pfam_tsv, --pgap_tsv")
        
        if args.sequence_dir or args.sequence_file:
            input_files = get_input_files(sequence_dir=args.sequence_dir, sequence_file=args.sequence_file, extension=args.extension)
            
            if args.sequence_file and len(input_files) == 1 and args.split_file_size:
                logger.info(f"Splitting {input_files[0]} into files with {args.split_file_size} sequences each")
                split_dir = Path(out_root) / "split_files"
                input_files = split_fasta_file(input_files[0], split_dir, args.split_file_size)
                logger.info(f"Created {len(input_files)} split files")
            
            annot_dir = Path(out_root) / "annot"
            annot_dir.mkdir(parents=True, exist_ok=True)
            scour_dir = Path(out_root) / "scours"
            scour_dir.mkdir(parents=True, exist_ok=True)
            
            targets = [(f.stem, str(f)) for f in input_files]
            
            def run_pfam_hmmsearch(t):
                sid, p = t
                tblout_path = run_hmmsearch(sid, p, str(annot_dir), args.pfam_db, args.threads, suffix='pfam', Z=args.Z_pfam, use_ga=True, nobias=True)
                scour_path = scour_dir / f"{sid}.pfam"
                parse_tblout_to_scour(tblout_path, str(scour_path))
            
            def run_pgap_hmmsearch(t):
                sid, p = t
                tblout_path = run_hmmsearch(sid, p, str(annot_dir), args.pgap_db, args.threads, suffix='pgap', Z=args.Z_pgap, use_ga=True, nobias=True)
                scour_path = scour_dir / f"{sid}.search"
                parse_tblout_to_scour(tblout_path, str(scour_path))
            
            logger.info(f"Running hmmsearch-g on {len(targets)} files with Pfam database")
            with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
                list(ex.map(run_pfam_hmmsearch, targets))
            
            logger.info(f"Running hmmsearch-g on {len(targets)} files with PGAP database")
            with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
                list(ex.map(run_pgap_hmmsearch, targets))
            
            logger.info("Merging PGAP scours with mapping file")
            pgap_output = Path(out_root) / "pgap_annotations.tsv"
            merge_scours_with_mapping(str(scour_dir), args.pgap_tsv, str(pgap_output), file_type='pgap')
            
            logger.info("Merging Pfam scours with mapping file")
            pfam_output = Path(out_root) / "pfam_annotations.tsv"
            merge_scours_with_mapping(str(scour_dir), args.pfam_tsv, str(pfam_output), file_type='pfam')
            
            logger.info("Merging PGAP and Pfam annotations into single file")
            merged_output = Path(out_root) / "merged_annotations.tsv"
            merge_pgap_and_pfam(str(pgap_output), str(pfam_output), str(merged_output))
        else:
            for domain in domains:
                annot_dir, targets = list_targets(out_root, domain) if not args.faa_dir else (os.path.join(out_root, domain, 'annot'), [(Path(p).stem, str(Path(args.faa_dir) / p)) for p in os.listdir(args.faa_dir) if p.endswith('.faa')])
                os.makedirs(annot_dir, exist_ok=True)
                if not targets:
                    continue
                
                scour_dir = Path(out_root) / domain / "scours"
                scour_dir.mkdir(parents=True, exist_ok=True)
                
                def run_pfam_hmmsearch(t):
                    sid, p = t
                    tblout_path = run_hmmsearch(sid, p, annot_dir, args.pfam_db, args.threads, suffix='pfam', Z=args.Z_pfam, use_ga=True, nobias=True)
                    scour_path = scour_dir / f"{sid}.pfam"
                    parse_tblout_to_scour(tblout_path, str(scour_path))
                
                def run_pgap_hmmsearch(t):
                    sid, p = t
                    tblout_path = run_hmmsearch(sid, p, annot_dir, args.pgap_db, args.threads, suffix='pgap', Z=args.Z_pgap, use_ga=True, nobias=True)
                    scour_path = scour_dir / f"{sid}.search"
                    parse_tblout_to_scour(tblout_path, str(scour_path))
                
                logger.info(f"Running hmmsearch-g on {len(targets)} {domain} samples with Pfam database")
                with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
                    list(ex.map(run_pfam_hmmsearch, targets))
                
                logger.info(f"Running hmmsearch-g on {len(targets)} {domain} samples with PGAP database")
                with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
                    list(ex.map(run_pgap_hmmsearch, targets))
                
                logger.info(f"Merging PGAP scours for {domain}")
                pgap_output = Path(out_root) / domain / "pgap_annotations.tsv"
                merge_scours_with_mapping(str(scour_dir), args.pgap_tsv, str(pgap_output), file_type='pgap')
                
                logger.info(f"Merging Pfam scours for {domain}")
                pfam_output = Path(out_root) / domain / "pfam_annotations.tsv"
                merge_scours_with_mapping(str(scour_dir), args.pfam_tsv, str(pfam_output), file_type='pfam')
                
                logger.info(f"Merging PGAP and Pfam annotations for {domain}")
                merged_output = Path(out_root) / domain / "merged_annotations.tsv"
                merge_pgap_and_pfam(str(pgap_output), str(pfam_output), str(merged_output))

    elif args.hmmdb:
        if not args.evalue_full and not args.evalue_dom:
            parser.error("When using --hmmdb, you must specify either --evalue_full or --evalue_dom")
        
        suffix = args.suffix
        
        if args.sequence_dir or args.sequence_file:
            input_files = get_input_files(sequence_dir=args.sequence_dir, sequence_file=args.sequence_file, extension=args.extension)
            
            if args.sequence_file and len(input_files) == 1 and args.split_file_size:
                logger.info(f"Splitting {input_files[0]} into files with {args.split_file_size} sequences each")
                split_dir = Path(out_root) / "split_files"
                input_files = split_fasta_file(input_files[0], split_dir, args.split_file_size)
                logger.info(f"Created {len(input_files)} split files")
            
            annot_dir = Path(out_root) / "annot"
            annot_dir.mkdir(parents=True, exist_ok=True)
            
            targets = [(f.stem, str(f)) for f in input_files]
            logger.info(f"Annotating {len(targets)} files with {args.hmmdb}")
            
            def task_user(t):
                sid, p = t
                return run_hmmsearch(sid, p, str(annot_dir), args.hmmdb, args.threads, suffix=suffix, Z=args.Z, use_ga=not args.no_cut_ga)
            
            with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
                list(ex.map(task_user, targets))
            
            summary_file = annot_dir.parent / f"annotations_{suffix}.tsv"
            with open(summary_file, 'w') as f:
                f.write('\t'.join(['sample_id', 'target_name', 'query_name', 'full_evalue', 'full_score', 'dom_evalue', 'dom_score']) + '\n')
            
            for sid, p in targets:
                hmm_tbl = annot_dir / f"{sid}.hmm.{suffix}.tsv"
                hmm_rows = parse_tblout_lines(str(hmm_tbl), evalue_full_cutoff=args.evalue_full, evalue_dom_cutoff=args.evalue_dom)
                for target_name, hmm in hmm_rows.items():
                    row_data = [sid, target_name, hmm.get('query_name', ''), hmm.get('full_evalue', ''), hmm.get('full_score', ''), hmm.get('dom_evalue', ''), hmm.get('dom_score', '')]
                    with open(summary_file, 'a') as out_f:
                        out_f.write('\t'.join(str(x) for x in row_data) + '\n')
            logger.info(f"Wrote annotation summary: {summary_file}")
        else:
            for domain in domains:
                annot_dir, targets = list_targets(out_root, domain) if not args.faa_dir else (os.path.join(out_root, domain, 'annot'), [(Path(p).stem, str(Path(args.faa_dir) / p)) for p in os.listdir(args.faa_dir) if p.endswith('.faa')])
                os.makedirs(annot_dir, exist_ok=True)
                if not targets:
                    continue
                logger.info(f"Annotating {len(targets)} {domain} samples with {args.hmmdb}")
                def task_user(t):
                    sid, p = t
                    return run_hmmsearch(sid, p, annot_dir, args.hmmdb, args.threads, suffix=suffix, Z=args.Z, use_ga=not args.no_cut_ga)
                with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
                    list(ex.map(task_user, targets))
                if domain == 'eukaryotes':
                    write_annotation_summary_euk_from_tblouts(out_root, suffix, evalue_full_cutoff=args.evalue_full, evalue_dom_cutoff=args.evalue_dom)
                else:
                    write_annotation_summary_bvm_from_tblouts(out_root, domain, suffix, evalue_full_cutoff=args.evalue_full, evalue_dom_cutoff=args.evalue_dom)
    else:
        parser.error("Must specify either (--pfam_db --pgap_db --pfam_tsv --pgap_tsv) for default mode or --hmmdb for user mode")

if __name__ == '__main__':
    main()


