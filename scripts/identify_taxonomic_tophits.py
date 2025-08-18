#!/usr/bin/env python3
import glob
import sys
import os
import gzip

def parse_file(path):
    entries = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#') or s.startswith('==>'):
                continue
            try:
                name, count_str = s.rsplit(None, 1)
                count = int(count_str)
                entries.append((name, count))
            except ValueError:
                continue
    return entries

def count_genome_bases(fasta_path):
    total = 0
    opener = gzip.open if fasta_path.endswith(".gz") else open
    try:
        with opener(fasta_path, 'rt', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('>'):
                    continue
                total += len(line.strip())
        return total
    except FileNotFoundError:
        return None

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: python findtophits.py /path/to/fa_dir/\n")
        sys.exit(1)

    fa_dir = sys.argv[1]
    ref_files = sorted(glob.glob("*ref"))

    # Header
    cols = (
        ["file", "genome_size"] +
        [f"count_{i}" for i in range(1, 6)] +
        [f"pct_{i}" for i in range(1, 6)] +
        [f"name_{i}" for i in range(1, 6)]
    )
    print("\t".join(cols))

    for ref_path in ref_files:
        entries = parse_file(ref_path)
        total_reads = sum(c for _, c in entries)

        # Top 5 hits
        top = sorted(entries, key=lambda x: x[1], reverse=True)[:5]
        while len(top) < 5:
            top.append(("", 0))

        counts = [c for _, c in top]
        pcts = [(c / total_reads * 100.0) if total_reads > 0 else 0.0 for c in counts]
        names = [n for n, _ in top]

        # Genome size: look for .fa then .fa.gz
        fa_base = os.path.basename(ref_path)[:-4]  # remove ".ref"
        fa_file = os.path.join(fa_dir, fa_base + ".fa")
        fa_gz_file = os.path.join(fa_dir, fa_base + ".fa.gz")
        genome_size = None
        if os.path.isfile(fa_file):
            genome_size = count_genome_bases(fa_file)
        elif os.path.isfile(fa_gz_file):
            genome_size = count_genome_bases(fa_gz_file)

        genome_size_str = str(genome_size) if genome_size is not None else "NA"

        row = (
            [os.path.basename(ref_path), genome_size_str] +
            [str(c) for c in counts] +
            [f"{p:.4f}" for p in pcts] +
            names
        )
        print("\t".join(row))

if __name__ == "__main__":
    main()

