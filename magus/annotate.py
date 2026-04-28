#!/usr/bin/env python3

import argparse
import csv
import logging
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

HMMSEARCH_BIN = 'hmmsearch-g-portable-static-best'


def get_input_files(sequence_dir=None, sequence_file=None, extension='faa'):
    if sequence_dir and sequence_file:
        raise ValueError('Cannot specify both --sequence-dir and --sequence-file')
    if not sequence_dir and not sequence_file:
        raise ValueError('Must specify either --sequence-dir or --sequence-file')

    if sequence_file:
        seq_path = Path(sequence_file)
        if not seq_path.is_file():
            raise ValueError(f'Sequence file does not exist: {sequence_file}')
        return [seq_path]

    seq_dir = Path(sequence_dir)
    if not seq_dir.is_dir():
        raise ValueError(f'Sequence directory does not exist: {sequence_dir}')

    ext = f'.{extension}' if not extension.startswith('.') else extension
    files = sorted(seq_dir.glob(f'*{ext}'))
    if not files:
        raise ValueError(f'No {ext} files found in {sequence_dir}')
    return files


def split_fasta_file(input_file, output_dir, sequences_per_file=100000):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    split_prefix = output_dir / f'{input_file.stem}_split'

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

    script_path = output_dir / 'split_script.awk'
    script_path.write_text(awk_script)
    subprocess.run(['awk', '-f', str(script_path), str(input_file)], check=True)
    script_path.unlink()

    split_files = sorted(split_prefix.parent.glob(f'{split_prefix.name}_*.faa'))
    return split_files


def run_hmmsearch(sample_id, protein_file, out_dir, hmm_db, threads, db_size):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tblout = out_dir / f'{sample_id}.hmm.tsv'
    log_file = out_dir / f'{sample_id}_hmmsearch.log'

    cmd = [
        HMMSEARCH_BIN,
        '-Z',
        str(db_size),
        '-o',
        '/dev/null',
        '--tblout',
        str(tblout),
        '--notextw',
        '--noali',
        '--nobias',
        '--cpu',
        str(threads),
        str(hmm_db),
        str(protein_file),
    ]

    with open(log_file, 'w') as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)

    return tblout


def parse_tblout_to_scour(tblout_path, scour_path):
    scour_path = Path(scour_path)
    scour_path.parent.mkdir(parents=True, exist_ok=True)

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

            outfile.write(f'{target_name}\t{query_accession}\t{full_evalue}\t{full_score}\t{dom_score}\n')


def parse_scour_file(scour_path):
    rows = {}
    if not Path(scour_path).exists():
        return rows

    with open(scour_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 5:
                continue

            rows[parts[0]] = {
                'query_name': parts[1],
                'full_evalue': parts[2],
                'full_score': parts[3],
                'dom_score': parts[4],
            }

    return rows


def load_tsv_mapping(tsv_path, key_col='source_identifier'):
    mapping = {}
    tsv_path = Path(tsv_path)
    if not tsv_path.exists():
        return mapping

    with open(tsv_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames:
            return mapping

        key_col_name = key_col
        if key_col_name not in reader.fieldnames and f'#{key_col_name}' in reader.fieldnames:
            key_col_name = f'#{key_col_name}'

        for row in reader:
            key = row.get(key_col_name, '')
            if key:
                mapping[key] = row

    return mapping


def merge_scours_with_mapping(scour_dir, mapping_tsv, output_file):
    mapping = load_tsv_mapping(mapping_tsv, key_col='source_identifier')
    if not mapping:
        logger.warning(f'Mapping file {mapping_tsv} is empty or not found')
        return

    scour_files = sorted(Path(scour_dir).glob('*.search'))
    if not scour_files:
        logger.warning(f'No *.search files found in {scour_dir}')
        return

    header = [
        'sample_id',
        'gene',
        'match',
        'evalue',
        'seq_score',
        'domain_score',
        'label',
        'gene_symbol',
        'product_name',
        'ec_numbers',
        'go_terms',
        'for_AMRFinder',
        'comment',
        'annotation_source',
    ]

    with open(output_file, 'w') as outfile:
        outfile.write('\t'.join(header) + '\n')

        for scour_file in scour_files:
            sample_id = scour_file.stem.replace('.search', '')
            scours = parse_scour_file(str(scour_file))

            for gene, hit in scours.items():
                match_id = hit['query_name']
                meta = mapping.get(match_id)
                if not meta:
                    continue

                try:
                    seq_score = float(hit['full_score'])
                    domain_score = float(hit['dom_score'])
                except ValueError:
                    continue

                seq_cutoff_str = meta.get('sequence_cutoff', '')
                dom_cutoff_str = meta.get('domain_cutoff', '')
                if not seq_cutoff_str or not dom_cutoff_str:
                    continue

                try:
                    seq_cutoff = float(seq_cutoff_str)
                    dom_cutoff = float(dom_cutoff_str)
                except (ValueError, TypeError):
                    continue

                if seq_score < seq_cutoff or domain_score < dom_cutoff:
                    continue

                row = [
                    sample_id,
                    gene,
                    match_id,
                    hit['full_evalue'],
                    hit['full_score'],
                    hit['dom_score'],
                    meta.get('label', ''),
                    meta.get('gene_symbol', ''),
                    meta.get('product_name', ''),
                    meta.get('ec_numbers', ''),
                    meta.get('go_terms', ''),
                    meta.get('for_AMRFinder', ''),
                    meta.get('comment', ''),
                    meta.get('annotation_source', ''),
                ]
                outfile.write('\t'.join(str(x) for x in row) + '\n')

    logger.info(f'Wrote merged annotations to {output_file}')


def write_annotation_summary(targets, annot_dir, output_file):
    with open(output_file, 'w') as f:
        f.write('\t'.join(['sample_id', 'target_name', 'query_name', 'full_evalue', 'full_score', 'dom_evalue', 'dom_score']) + '\n')

    for sample_id, _ in targets:
        hmm_tbl = Path(annot_dir) / f'{sample_id}.hmm.tsv'
        if not hmm_tbl.exists():
            continue
        with open(hmm_tbl, 'r') as f:
            for line in f:
                if not line or line.startswith('#'):
                    continue
                parts = line.strip().split()
                if len(parts) < 9:
                    continue

                row = [
                    sample_id,
                    parts[0],
                    parts[2] if len(parts) > 2 else '',
                    parts[4] if len(parts) > 4 else '',
                    parts[5] if len(parts) > 5 else '',
                    parts[7] if len(parts) > 7 else '',
                    parts[8] if len(parts) > 8 else '',
                ]
                with open(output_file, 'a') as out_f:
                    out_f.write('\t'.join(row) + '\n')


def main():
    parser = argparse.ArgumentParser(description='Annotate proteins with merged PGAP+Pfam database.')
    parser.add_argument('--sequence-dir', type=str, default=None, dest='sequence_dir', help='Directory containing FASTA files to annotate.')
    parser.add_argument('--sequence-file', type=str, default=None, dest='sequence_file', help='Single FASTA file to annotate.')
    parser.add_argument('--split-file-size', type=int, default=100000, help='Number of sequences per split file when splitting --sequence-file.')
    parser.add_argument('-x', '--extension', type=str, default='faa', help='File extension when using --sequence-dir (default: faa).')
    parser.add_argument('--output-directory', '--output_directory', type=str, default='magus_output/orf_calling', dest='output_directory', help='Directory for output files.')
    parser.add_argument('--threads', type=int, default=8, help='Threads per hmmsearch job.')
    parser.add_argument('--maxworkers', '--max-workers', dest='max_workers', type=int, default=4, help='Number of files to process in parallel.')
    parser.add_argument('--annotation-tsv', '--annotations', '--annntoation', dest='annotation_tsv', required=True, help='Merged PGAP/Pfam annotation TSV.')
    parser.add_argument('--hmmdb', required=True, help='Path to merged PGAP/Pfam HMM database.')
    parser.add_argument('--db-size', type=int, default=48647, help='Value passed to hmmsearch -Z (default: 48647).')

    args = parser.parse_args()

    input_files = get_input_files(sequence_dir=args.sequence_dir, sequence_file=args.sequence_file, extension=args.extension)

    if args.sequence_file and args.split_file_size:
        logger.info(f'Splitting {input_files[0]} into files with {args.split_file_size} sequences each')
        split_dir = Path(args.output_directory) / 'split_files'
        input_files = split_fasta_file(input_files[0], split_dir, args.split_file_size)
        logger.info(f'Created {len(input_files)} split files')

    annot_dir = Path(args.output_directory) / 'annot'
    annot_dir.mkdir(parents=True, exist_ok=True)

    scour_dir = Path(args.output_directory) / 'scours'
    scour_dir.mkdir(parents=True, exist_ok=True)

    targets = [(f.stem, str(f)) for f in input_files]

    def task(target):
        sample_id, protein_path = target
        tblout_path = run_hmmsearch(sample_id, protein_path, annot_dir, args.hmmdb, args.threads, db_size=args.db_size)
        scour_path = scour_dir / f'{sample_id}.search'
        parse_tblout_to_scour(tblout_path, scour_path)

    total_targets = len(targets)
    logger.info(f'Running {HMMSEARCH_BIN} on {total_targets} sequence files')
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        future_to_target = {executor.submit(task, target): target for target in targets}
        for idx, future in enumerate(as_completed(future_to_target), start=1):
            sample_id, _ = future_to_target[future]
            future.result()
            pct_complete = (idx / total_targets) * 100
            logger.info(f'Annotation progress: {idx}/{total_targets} files ({pct_complete:.1f}%) complete - {sample_id}')

    merged_output = Path(args.output_directory) / 'merged_annotations.tsv'
    merge_scours_with_mapping(scour_dir, args.annotation_tsv, merged_output)

    summary_output = Path(args.output_directory) / 'annotations.tsv'
    write_annotation_summary(targets, annot_dir, summary_output)
    logger.info(f'Wrote annotation summary: {summary_output}')


if __name__ == '__main__':
    main()
