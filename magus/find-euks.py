import os
import glob
import argparse
import subprocess
import pandas as pd
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class EukRepRunner:
    def __init__(self, bin_dirs, wildcards, size_threshold, euk_binning_outputdir, dblocs, max_workers=4, threads=8, skip_eukrep=False, skip_eukcc=False, eukrepenv=None, checkm2_file=None):
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
        self.checkm2_file = checkm2_file

        self.input_bins_dir = os.path.join(self.euk_binning_outputdir, "input_bins")
        os.makedirs(self.input_bins_dir, exist_ok=True)
        self.bin_sizes = {}
        self.bin_contigs = {}
        self.bin_name_mapping = {}  # Maps short names to full names
        self.checkm2_data = None

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
                logging.info(f"Symlinked {bin_path} to {symlink_dest}")
            else:
                logging.info(f"Symlink already exists for {bin_path}")

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
                
                # Check if EukRep output already exists and is not empty
                if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                    logging.info(f"EukRep output already exists for {bin_name}, skipping")
                    continue
                
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
                eukcc_file = os.path.join(output_dir, 'eukcc.csv')
                
                # Check if EukCC output already exists and is not empty
                if os.path.exists(eukcc_file) and os.path.getsize(eukcc_file) > 0:
                    logging.info(f"EukCC output already exists for {bin_name}, skipping")
                    continue
                
                cmd = f"eukcc single --out {output_dir} --threads {self.threads} --db {self.eukcc_db} {bin_file}"
                futures.append(executor.submit(subprocess.run, cmd, shell=True))

            for future in as_completed(futures):
                future.result()
        print("EukCC completed.")

    def find_checkm2_outputs(self):
        """Find CheckM2 quality report files in bin directories."""
        checkm2_files = []
        for directory in self.bin_dirs:
            for root, dirs, files in os.walk(directory):
                if 'quality_report.tsv' in files:
                    checkm2_files.append(os.path.join(root, 'quality_report.tsv'))
        
        if checkm2_files:
            logging.info(f"Found {len(checkm2_files)} CheckM2 quality report files")
            return checkm2_files
        else:
            logging.info("No CheckM2 quality report files found")
            return []

    def create_bin_name_mapping(self):
        """Create mapping between short bin names and full bin names."""
        for bin_path in self.bin_paths:
            full_name = os.path.basename(bin_path).replace('.fa', '')
            # Extract the last part of the name (e.g., "30.7" from the long name)
            short_name = full_name.split('-')[-1]
            self.bin_name_mapping[short_name] = full_name
        logging.info(f"Created mapping for {len(self.bin_name_mapping)} bins")

    def load_checkm2_data(self):
        """Load and process CheckM2 data."""
        if self.checkm2_file:
            logging.info(f"Loading user-provided CheckM2 file: {self.checkm2_file}")
            self.checkm2_data = pd.read_csv(self.checkm2_file, sep='\t')
            return

        checkm2_files = self.find_checkm2_outputs()
        if not checkm2_files:
            return

        # Combine all CheckM2 files
        dfs = []
        for file in checkm2_files:
            df = pd.read_csv(file, sep='\t')
            dfs.append(df)
        
        if dfs:
            self.checkm2_data = pd.concat(dfs, ignore_index=True)
            logging.info(f"Loaded CheckM2 data for {len(self.checkm2_data)} bins")
        else:
            logging.info("No CheckM2 data found to load")

    def process_output(self):
        summary_data = []
        for bin_file in glob.glob(os.path.join(self.input_bins_dir, "*.fa")):
            bin_name = os.path.basename(bin_file).replace('.fa', '')
            output_dir = os.path.join(self.euk_binning_outputdir, bin_name)
            eukcc_file = os.path.join(output_dir, 'eukcc/eukcc.csv')
            contig_file = os.path.join(output_dir, f"EUKREP_{bin_name}_eukrepcontigs.fa")

            completeness = contamination = float('nan')
            if os.path.exists(eukcc_file) and os.path.getsize(eukcc_file) > 0:
                try:
                    eukcc_data = pd.read_csv(eukcc_file, sep='\t')
                    if 'completeness' in eukcc_data.columns:
                        completeness = eukcc_data['completeness'].iloc[0]
                    if 'contamination' in eukcc_data.columns:
                        contamination = eukcc_data['contamination'].iloc[0]
                except Exception as e:
                    logging.warning(f"Error reading EukCC file for {bin_name}: {str(e)}")

            eukrep_contig_count = 0
            if os.path.exists(contig_file) and os.path.getsize(contig_file) > 0:
                try:
                    with open(contig_file, 'r') as f:
                        eukrep_contig_count = sum(1 for line in f if line.startswith('>'))
                except Exception as e:
                    logging.warning(f"Error reading EukRep file for {bin_name}: {str(e)}")

            row = {
                'bin_name': bin_name,
                'bin_size': self.bin_sizes.get(bin_name, float('nan')),
                'bin_contig_count': self.bin_contigs.get(bin_name, float('nan')),
                'completeness': completeness,
                'contamination': contamination,
                'eukrep_contig_count': eukrep_contig_count
            }

            # Add CheckM2 data if available
            if self.checkm2_data is not None:
                short_name = bin_name.split('-')[-1]
                checkm2_row = self.checkm2_data[self.checkm2_data.iloc[:, 0] == short_name]
                if not checkm2_row.empty:
                    for col in checkm2_row.columns[1:]:  # Skip the bin name column
                        row[f'checkm2_{col}'] = checkm2_row[col].iloc[0]

            summary_data.append(row)

        summary_df = pd.DataFrame(summary_data)
        output_file = os.path.join(self.euk_binning_outputdir, 'eukaryotic_summary_table.csv')
        summary_df.to_csv(output_file, index=False)
        logging.info(f"Summary table saved to {output_file} with CheckM2 data incorporated")

    def run(self):
        self.find_bins()
        self.create_bin_name_mapping()
        self.load_checkm2_data()
        if not self.skip_eukrep:
            self.run_eukrep()
        if not self.skip_eukcc:
            self.run_eukcc()
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
    parser.add_argument("--skip_eukrep", action='store_true', help="Skip EukRep step")
    parser.add_argument("--skip_eukcc", action='store_true', help="Skip EukCC step")
    parser.add_argument("--eukrep_env", type=str, default=None)
    parser.add_argument("--checkm2_file", type=str, help="Optional: Direct path to CheckM2 quality report file")
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
        eukrepenv=args.eukrep_env,
        checkm2_file=args.checkm2_file
    )
    runner.run()
