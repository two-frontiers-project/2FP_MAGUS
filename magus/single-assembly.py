import subprocess
import os
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

class Assembly:
    def __init__(self, config, outdir="asm", threads=14):
        self.config = self.load_config(config)
        self.outdir = outdir
        self.threads = threads
        os.makedirs(self.outdir, exist_ok=True)

    def load_config(self, config_path):
        config_df = pd.read_csv(config_path, sep='\t')
        return config_df.to_dict(orient='records')

    def run_megahit(self, sample):
        sample_name = sample['filename']
        r1 = sample['pe1'] #f"{sample_name}.R1.fa.gz"
        r2 = sample['pe2'] #f"{sample_name}.R2.fa.gz"
        out_sample_dir = f"{self.outdir}/{sample_name}"

        # Run megahit assembly
        cmd = (f"megahit --k-min 75 --k-max 333 --k-step 6 --cleaning-rounds 1 "
               f"--merge-level 100,.999 --min-count 1 --min-contig-len 1000 "
               f"--continue -t {self.threads} -1 {r1} -2 {r2} -o {out_sample_dir}")

        print(f"Assembling {sample_name}")
        subprocess.run(cmd, shell=True)

        # Post-process the contigs
        self.post_process_contigs(sample_name)

    def post_process_contigs(self, sample_name):
        out_sample_dir = f"{self.outdir}/{sample_name}"
        contig_file = f"{out_sample_dir}/final.contigs.fa"
        filtered_contig_file = f"{out_sample_dir}/temp0.fa"

        # Filter contigs based on length and multiplicity
        cmd = (f"grep -A1 --no-group-separator 'multi=1\\.[5-9]\\|multi=[2-9]\\|multi=[0-9][0-9]' {contig_file} "
               f"| awk '{{if (0==(NR % 2) && length >= 1000) {{print x; print $0}}; x=$0}}' > {filtered_contig_file}")

        subprocess.run(cmd, shell=True)
        print(f"Filtered contigs for {sample_name}")

    def run(self, max_workers=None):
        # Run assemblies in parallel using ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            executor.map(self.run_megahit, self.config)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Run assembly with megahit on genomic data.")
    parser.add_argument('--config', type=str, required=True, help='Path to the configuration TSV file')
    parser.add_argument('--max_workers', type=int, default=1, help='Number of parallel workers (default: 1)')
    parser.add_argument('--threads', type=int, default=14, help='Number of threads for megahit (default: 14)')

    # Parse arguments
    args = parser.parse_args()

    # Initialize and run the Assembly instance
    assembly = Assembly(config=args.config, threads=args.threads)
    assembly.run(max_workers=args.max_workers)

if __name__ == '__main__':
    main()
