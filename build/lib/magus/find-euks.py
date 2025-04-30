import os
import glob
import argparse
import subprocess
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

class EukRepRunner:
    def __init__(self, bin_dirs, wildcards, size_threshold, euk_binning_outputdir, dblocs, max_workers=4, threads=8, skip_eukrep=False, skip_eukcc=False, eukrepenv=None):
        self.bin_dirs = bin_dirs
        self.wildcards = wildcards
        self.size_threshold = size_threshold
        self.euk_binning_outputdir = euk_binning_outputdir
        self.eukcc_db = self.get_db_location(dblocs, 'eukccdb')
        self.max_workers = max_workers
        self.threads = threads
        self.skip_eukrep = skip_eukrep
        self.skip_eukcc = skip_eukcc
        self.eukrepenv = eukrepenv

        self.input_bins_dir = os.path.join(self.euk_binning_outputdir, "input_bins")
        os.makedirs(self.input_bins_dir, exist_ok=True)
        self.bin_sizes = {}
        self.bin_contigs = {}

    def get_db_location(self, dblocs, db_name):
        db_df = pd.read_csv(dblocs, header=None, index_col=0)
        return db_df.loc[db_name, 1]

    def find_bins(self):
        self.bin_paths = []
        for directory in self.bin_dirs:
            for pattern in self.wildcards:
                pattern="*" + pattern +"*"
                self.bin_paths.extend(glob.glob(os.path.join(directory, pattern)))
        
        for bin_path in self.bin_paths:
            unique_name = os.path.relpath(bin_path).rsplit('.', 1)[0].replace('/', '-')
            is_large, bin_size, contig_count = self.is_bin_large(bin_path)
            self.bin_sizes[unique_name] = bin_size
            self.bin_contigs[unique_name] = contig_count

            symlink_dest = os.path.join(self.input_bins_dir, unique_name.replace('/','-') + ".fa")
            if not os.path.exists(symlink_dest):
                os.symlink(os.path.abspath(bin_path), symlink_dest)
                print(f"Symlinked {bin_path} to {symlink_dest}")

    def is_bin_large(self, bin_path):
        bin_size = 0
        contig_count = 0
        with open(bin_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    contig_count += 1
                else:
                    bin_size += len(line.strip())
        return bin_size >= self.size_threshold, bin_size, contig_count

    def run_eukrep(self):
        futures = []
        print("Running EukRep...")
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            for bin_file in glob.glob(os.path.join(self.input_bins_dir, "*.fa")):
                bin_name = os.path.basename(bin_file).replace('.fa', '')
                output_dir = os.path.join(self.euk_binning_outputdir, bin_name)
                os.makedirs(output_dir, exist_ok=True)
                output_file = os.path.join(output_dir, f"EUKREP_{bin_name}_eukrepcontigs.fa")
                cmd = f"bash -c 'source activate {self.eukrepenv} && EukRep -i {bin_file} -o {output_file}'"
                futures.append(executor.submit(subprocess.run, cmd, shell=True))

            for future in as_completed(futures):
                future.result()
        print("EukRep completed.")

    def run_eukcc(self):
        futures = []
        print("Running EukCC...")
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            for bin_file in glob.glob(os.path.join(self.input_bins_dir, "*.fa")):
                bin_name = os.path.basename(bin_file).replace('.fa', '')
                output_dir = os.path.join(self.euk_binning_outputdir, bin_name, 'eukcc')
                os.makedirs(output_dir, exist_ok=True)
                cmd = f"eukcc single --out {output_dir} --threads {self.threads} --db {self.eukcc_db} {bin_file}"
                futures.append(executor.submit(subprocess.run, cmd, shell=True))

            for future in as_completed(futures):
                future.result()
        print("EukCC completed.")

    def process_output(self):
        summary_data = []
        for bin_file in glob.glob(os.path.join(self.input_bins_dir, "*.fa")):
            bin_name = os.path.basename(bin_file).replace('.fa', '')
            output_dir = os.path.join(self.euk_binning_outputdir, bin_name)
            eukcc_file = os.path.join(output_dir, 'eukcc/eukcc.csv')
            contig_file = os.path.join(output_dir, f"EUKREP_{bin_name}_eukrepcontigs.fa")

            completeness = contamination = float('nan')
            if os.path.exists(eukcc_file) and os.path.getsize(eukcc_file) > 0:
                eukcc_data = pd.read_csv(eukcc_file, sep='\t')
                if 'completeness' in eukcc_data.columns:
                    completeness = eukcc_data['completeness'].iloc[0]
                if 'contamination' in eukcc_data.columns:
                    contamination = eukcc_data['contamination'].iloc[0]

            eukrep_contig_count = 0
            if os.path.exists(contig_file):
                with open(contig_file, 'r') as f:
                    eukrep_contig_count = sum(1 for line in f if line.startswith('>'))

            row = {
                'bin_name': bin_name,
                'bin_size': self.bin_sizes.get(bin_name, float('nan')),
                'bin_contig_count': self.bin_contigs.get(bin_name, float('nan')),
                'completeness': completeness,
                'contamination': contamination,
                'eukrep_contig_count': eukrep_contig_count
            }
            summary_data.append(row)

        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(os.path.join(self.euk_binning_outputdir, 'eukaryotic_summary_table.csv'), index=False)
        print("Summary table saved.")

    def run(self):
        self.find_bins()
        #if not self.skip_eukrep:
        #    self.run_eukrep()
        #if not self.skip_eukcc:
        #    self.run_eukcc()
        self.process_output()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin_dirs", type=str, required=True, help="Pipe-separated list of directories containing bins, quoted")
    parser.add_argument("--wildcards", type=str, default = "", required=False, help="Pipe-separated list of glob patterns for bin files, quoted")
    parser.add_argument("--size_threshold", type=int, default=10000000)
    parser.add_argument("--euk_binning_outputdir", type=str, default="magus_output/magus_euks")
    parser.add_argument("--dblocs", type=str, required=True)
    parser.add_argument("--max_workers", type=int, default=1)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--skip_eukrep", type=lambda x: x.lower() == 'true', default=False)
    parser.add_argument("--skip_eukcc", type=lambda x: x.lower() == 'true', default=False)
    parser.add_argument("--eukrep_env", type=str, default=None)
    args = parser.parse_args()

    bin_dirs = args.bin_dirs.split('|')
    wildcards = args.wildcards.split('|')

    runner = EukRepRunner(
        bin_dirs=bin_dirs,
        wildcards=wildcards,
        size_threshold=args.size_threshold,
        euk_binning_outputdir=args.euk_binning_outputdir,
        dblocs=args.dblocs,
        max_workers=args.max_workers,
        threads=args.threads,
        skip_eukrep=args.skip_eukrep,
        skip_eukcc=args.skip_eukcc,
        eukrepenv=args.eukrep_env
    )
    runner.run()
