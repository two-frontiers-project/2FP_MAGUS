import os
import subprocess
import pandas as pd
import argparse
from concurrent.futures import ThreadPoolExecutor

class Sample:
    def __init__(self, name, pe1, pe2=None):
        self.name = name
        self.pe1 = pe1
        self.pe2 = pe2

    def is_paired(self):
        return bool(self.pe2)

class XTreeRunner:
    def __init__(self, sample, output_dir, db_path, threads):
        self.sample = sample
        self.output_dir_base = output_dir
        self.output_dir = self.output_dir_base + '/raw_alignments'
        self.db_path = db_path
        self.threads = threads

    def run(self):
        if not os.path.exists(self.sample.pe1):
            raise FileNotFoundError(f"File not found: {self.sample.pe1}")
        if self.sample.is_paired() and not os.path.exists(self.sample.pe2):
            raise FileNotFoundError(f"File not found: {self.sample.pe2}")

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        cov_out = os.path.join(self.output_dir, f"{self.sample.name}.cov")
        perq_out = os.path.join(self.output_dir, f"{self.sample.name}.perq")

        # Determine command based on file compression
        def decompress_cmd(filepath):
            return f"zcat {filepath}" if filepath.endswith(".gz") else f"cat {filepath}"

        if self.sample.is_paired():
            cmd = f"{decompress_cmd(self.sample.pe1)} && {decompress_cmd(self.sample.pe2)}"
            cmd = f"cat <({decompress_cmd(self.sample.pe1)}) <({decompress_cmd(self.sample.pe2)})"
        else:
            cmd = decompress_cmd(self.sample.pe1)

        cmd += f" | xtree ALIGN --seqs - --cov-out {cov_out} --db {self.db_path} --perq-out {perq_out} --threads {self.threads} --redistribute --doforage"

        print(f"Running xtree for sample {self.sample.name}:")
        print(cmd)
        subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)

def merge_abundance_matrix(output_dir, tax_file=None, coverage_cutoff=0.05):
    import glob

    alignment_dir = os.path.join(output_dir, "raw_alignments")
    cov_files = glob.glob(os.path.join(alignment_dir, "*.cov"))

    if not cov_files:
        print("No .cov files found for merging.")
        return

    dfs = []
    for file in cov_files:
        df = pd.read_csv(file, sep='\t')
        if df.empty:
            continue
        sample_name = os.path.splitext(os.path.basename(file))[0]
        df['sample'] = sample_name
        df['min_coverage'] = df[[
            'Adamantium_covered',
            'Unique_proportion_covered',
            'Proportion_covered']].min(axis=1)
        df['RA'] = df['Coverage_est'] / df['Coverage_est'].sum()
        dfs.append(df)

    if not dfs:
        print("All .cov files were empty after filtering.")
        return

    merged_df = pd.concat(dfs, ignore_index=True)
    merged_df = merged_df[merged_df['min_coverage'] >= coverage_cutoff]

    if tax_file:
        taxonomy_df = pd.read_csv(tax_file, sep='\t', header=None, names=['Reference', 'taxonomy'])
        taxonomy_df['Reference'] = taxonomy_df['Reference'].str.replace('RS_', '', regex=False).str.replace('GB_', '', regex=False)
        merged_df = pd.merge(merged_df, taxonomy_df, on='Reference', how='inner')

    abundance_matrix = merged_df.pivot_table(index='Reference', columns='sample', values='RA', fill_value=0)
    abundance_matrix.to_csv(os.path.join(output_dir, "merged_xtree.csv"))
    print("Merged abundance matrix written to merged_xtree.csv")

def load_samples(config_path):
    df = pd.read_csv(config_path, sep='\t')
    samples = []
    for _, row in df.iterrows():
        filename = str(row['filename'])
        pe1 = str(row['pe1'])
        pe2 = str(row['pe2']) if pd.notna(row['pe2']) and row['pe2'].strip() else None
        samples.append(Sample(filename, pe1, pe2))
    return samples

def parse_args():
    parser = argparse.ArgumentParser(description="Run XTree to get taxonomic classification on SE/PE reads using parallel processing.")
    parser.add_argument("--config", required=True, help="Path to TSV config file with columns: filename, pe1, pe2 (pe2 can just be empty for se files)")
    parser.add_argument("--output", default='magus/xtree_output', help="Directory to store output files (default magus/xtree_output)")
    parser.add_argument("--db", required=True, help="Path to relevant xtree database")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use for each xtree run (default 4)")
    parser.add_argument("--max_workers", type=int, default=1, help="Max number of samples to process in parallel (default 1)")
    parser.add_argument("--taxmap", help="Path to GTDB taxonomy file; if not provided, taxonomy names will not be merged")
    parser.add_argument("--coverage-cutoff", type=float, default=0.05, help="Coverage cutoff for filtering alignments (default 0.05)")
    return parser.parse_args()

def main():
    args = parse_args()

    samples = load_samples(args.config)

    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = []
        for sample in samples:
            runner = XTreeRunner(sample, args.output, args.db, args.threads)
            futures.append(executor.submit(runner.run))

        for future in futures:
            try:
                future.result()
            except Exception as e:
                print(f"Error processing sample: {e}")

    merge_abundance_matrix(args.output, tax_file=args.taxmap, coverage_cutoff=args.coverage_cutoff)

if __name__ == "__main__":
    main()
