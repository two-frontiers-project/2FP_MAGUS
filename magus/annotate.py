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

def parse_scour_file(scour_path):
    """
    Parse 5-column TSVs produced by the provided shell script:
    columns: target_name, query_name, query_accession, full_evalue, metric (score/bias)
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
            # We store the commonly used fields; dom_* not available in scours
            rows[target_name] = {
                'query_name': query_name,
                'full_evalue': parts[3],
                'full_score': parts[4],
                'dom_evalue': '',
                'dom_score': ''
            }
    return rows

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
    parser.add_argument('--domains', type=str, default='bacteria,viruses,metagenomes,eukaryotes', help='Comma-separated domains to process.')
    parser.add_argument('--threads', type=int, default=8, help='Threads per hmmsearch-g job.')
    parser.add_argument('--max_workers', type=int, default=4, help='Parallel samples.')

    subparsers = parser.add_subparsers(dest='mode', required=True)

    sp_default = subparsers.add_parser('default', help='Merge from precomputed PFAM/PGAP annotation TSVs, or run hmmsearch-g if DBs are provided.')
    sp_default.add_argument('--pfam_tsv_dir', type=str, default=None, help='Directory of PFAM scours TSVs (e.g., prot/scours/*.pfam).')
    sp_default.add_argument('--pgap_tsv_dir', type=str, default=None, help='Directory of PGAP scours TSVs (e.g., prot/scours/*.search).')
    sp_default.add_argument('--pfam_db', type=str, default=None, help='Optional: Path to Pfam-A HMM database to run hmmsearch-g.')
    sp_default.add_argument('--pgap_db', type=str, default=None, help='Optional: Path to PGAP HMM database to run hmmsearch-g.')
    sp_default.add_argument('--Z_pfam', type=int, default=25545, help='Pfam database size for -Z.')
    sp_default.add_argument('--Z_pgap', type=int, default=18057, help='PGAP database size for -Z.')
    sp_default.add_argument('--no_cut_ga', action='store_true', help='Do not apply --cut_ga.')

    sp_user = subparsers.add_parser('user', help='Run user-specified HMM database and optional e-value cutoffs.')
    sp_user.add_argument('--hmmdb', type=str, required=True, help='Path to HMM database.')
    sp_user.add_argument('--suffix', type=str, default='custom', help='Suffix tag for outputs.')
    sp_user.add_argument('--evalue_full', type=float, default=None, help='Full-sequence E-value cutoff.')
    sp_user.add_argument('--evalue_dom', type=float, default=None, help='Best-domain E-value cutoff.')
    sp_user.add_argument('--no_cut_ga', action='store_true', help='Do not apply --cut_ga.')
    sp_user.add_argument('--Z', type=int, default=None, help='Database size for -Z.')

    args = parser.parse_args()

    domains = [d.strip() for d in args.domains.split(',') if d.strip()]
    out_root = args.output_directory

    if args.mode == 'default':
        # If scour TSVs are provided, merge those; otherwise, run hmmsearch-g only in user mode.
        for domain in domains:
            annot_dir, targets = list_targets(out_root, domain) if not args.faa_dir else (os.path.join(out_root, domain, 'annot'), [(Path(p).stem, str(Path(args.faa_dir) / p)) for p in os.listdir(args.faa_dir) if p.endswith('.faa')])
            os.makedirs(annot_dir, exist_ok=True)
            if not targets:
                continue
            # Merge from provided scours; do not run hmmsearch-g in default mode
            if not args.pgap_tsv_dir and not args.pfam_tsv_dir:
                logger.warning("Default mode requires --pgap_tsv_dir and/or --pfam_tsv_dir. Skipping domain.")
                continue
            if args.pgap_tsv_dir:
                logger.info(f"Merging PGAP scours for {domain} from {args.pgap_tsv_dir}")
                if domain == 'eukaryotes':
                    write_annotation_summary_euk_from_scours(out_root, 'pgap', args.pgap_tsv_dir)
                else:
                    write_annotation_summary_bvm_from_scours(out_root, domain, 'pgap', args.pgap_tsv_dir)
            if args.pfam_tsv_dir:
                logger.info(f"Merging PFAM scours for {domain} from {args.pfam_tsv_dir}")
                if domain == 'eukaryotes':
                    write_annotation_summary_euk_from_scours(out_root, 'pfam', args.pfam_tsv_dir)
                else:
                    write_annotation_summary_bvm_from_scours(out_root, domain, 'pfam', args.pfam_tsv_dir)

    elif args.mode == 'user':
        suffix = args.suffix
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

if __name__ == '__main__':
    main()


