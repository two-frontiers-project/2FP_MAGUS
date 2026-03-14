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
    def __init__(self, bin_dirs, asmdir, wildcards, size_threshold, euk_binning_outputdir, dblocs=None, max_workers=4, threads=8, run_eukrep=False, run_eukcc=False, eukrepenv="eukrep-env", checkm2_file=None):
        self.bin_dirs = bin_dirs
        self.asmdir = asmdir
        self.wildcards = wildcards
        self.size_threshold = size_threshold
        self.euk_binning_outputdir = euk_binning_outputdir
        self.dblocs = dblocs
        self.eukcc_db = None
        self.max_workers = max_workers
        self.threads = threads
        self.run_eukrep_flag = run_eukrep
        self.run_eukcc_flag = run_eukcc
        self.eukrepenv = eukrepenv
        self.checkm2_file = checkm2_file

        self.input_bins_dir = os.path.join(self.euk_binning_outputdir, "tmp_input_bins")
        os.makedirs(self.input_bins_dir, exist_ok=True)
        self.bin_sizes = {}
        self.bin_contigs = {}
        self.bin_name_mapping = {}  # Maps short names to full names
        self.checkm2_data = None

    def get_db_location(self, dblocs, db_name):
        db_df = pd.read_csv(dblocs, header=None, index_col=0)
        return db_df.loc[db_name, 1]

    def discover_default_bins_from_asmdir(self):
        """Discover bins from single-assembly output directory structure: <asmdir>/*/bins/*.fa."""
        if not self.asmdir:
            return []

        pattern = os.path.join(self.asmdir, "*", "bins", "*.fa")
        discovered_bins = glob.glob(pattern)
        logging.info(f"Discovered {len(discovered_bins)} bins from asmdir pattern: {pattern}")
        return discovered_bins

    def find_bins(self):
        self.bin_paths = []

        # Default workflow for MAGUS users after single-assembly binning.
        if self.bin_dirs is None:
            self.bin_paths = self.discover_default_bins_from_asmdir()
            logging.info(f"Using asmdir-based bin discovery from: {self.asmdir}")
        else:
            self.bin_paths = []
        # Explicit bin_dirs workflow (legacy mode)
        if self.bin_dirs is not None:
            # First expand any glob patterns in bin_dirs to get actual directory paths
            expanded_dirs = []
            for directory in self.bin_dirs:
                if '*' in directory:
                    # Expand glob pattern to get actual directories
                    expanded_dirs.extend(glob.glob(directory))
                else:
                    expanded_dirs.append(directory)

            logging.info(f"Expanded directories: {expanded_dirs}")

            # Process each expanded directory
            for directory in expanded_dirs:
                logging.info(f"Searching in directory: {directory}")
                if not os.path.exists(directory):
                    logging.warning(f"Directory does not exist: {directory}")
                    continue

                for pattern in self.wildcards:
                    original_pattern = pattern
                    logging.info(f"Processing wildcard: {original_pattern}")

                    # Check if pattern looks like a directory name (no file extension)
                    if '.' not in pattern and pattern in ['bins', 'bin']:
                        # Special case: look for files within directories named 'bins'
                        bins_dirs = glob.glob(os.path.join(directory, "*", pattern), recursive=False)
                        logging.info(f"Found {len(bins_dirs)} 'bins' directories under {directory}")

                        for bins_dir in bins_dirs:
                            logging.info(f"Searching for files in bins directory: {bins_dir}")
                            # Look for common genome file extensions in bins directory
                            for ext in ['*.fa', '*.fasta', '*.fna', '*.fas']:
                                files = glob.glob(os.path.join(bins_dir, ext))
                                logging.info(f"  Found {len(files)} files matching {ext}")
                                for f in files:
                                    logging.info(f"    File: {f}")
                                self.bin_paths.extend(files)
                    else:
                        # Regular file pattern matching
                        pattern = "*" + pattern + "*"
                        logging.info(f"Looking for file pattern: {pattern}")

                        # First try direct pattern match
                        direct_pattern = os.path.join(directory, pattern)
                        logging.info(f"Direct search pattern: {direct_pattern}")
                        direct_matches = glob.glob(direct_pattern)
                        logging.info(f"Direct matches found: {len(direct_matches)}")
                        for match in direct_matches:
                            logging.info(f"  Direct match: {match}")
                        self.bin_paths.extend(direct_matches)

                        # Also try recursive search for nested directories
                        recursive_pattern = os.path.join(directory, "**", pattern)
                        logging.info(f"Recursive search pattern: {recursive_pattern}")
                        recursive_matches = glob.glob(recursive_pattern, recursive=True)
                        logging.info(f"Recursive matches found: {len(recursive_matches)}")
                        for match in recursive_matches:
                            logging.info(f"  Recursive match: {match}")
                        self.bin_paths.extend(recursive_matches)
        
        logging.info(f"Total bin paths found: {len(self.bin_paths)}")

        # Remove duplicates
        self.bin_paths = list(set(self.bin_paths))

        for bin_path in self.bin_paths:
            unique_name = self.create_unique_bin_name(bin_path)
            is_large, bin_size, contig_count = self.is_bin_large(bin_path)
            
            # Only process bins that meet the size threshold
            if is_large:
                self.bin_sizes[unique_name] = bin_size
                self.bin_contigs[unique_name] = contig_count

                symlink_dest = os.path.join(self.input_bins_dir, unique_name.replace('/','-') + ".fa")
                if not os.path.exists(symlink_dest):
                    os.symlink(os.path.abspath(bin_path), symlink_dest)
                    logging.info(f"Symlinked {bin_path} to {symlink_dest} (size: {bin_size:,} bp)")
                else:
                    logging.info(f"Symlink already exists for {bin_path}")
            else:
                logging.info(f"Skipping {bin_path} - size {bin_size:,} bp below threshold {self.size_threshold:,} bp")

    def create_unique_bin_name(self, bin_path):
        """Create a stable, non-hidden bin name for symlinks and downstream merges."""
        if self.bin_dirs is None and self.asmdir:
            asmdir_base = os.path.basename(os.path.normpath(self.asmdir)) or "asm"
            rel = os.path.relpath(bin_path, start=self.asmdir).rsplit('.', 1)[0]
            unique_name = f"{asmdir_base}-{rel.replace('/', '-')}"
        else:
            unique_name = os.path.relpath(bin_path).rsplit('.', 1)[0].replace('/', '-')

        return unique_name.lstrip('.-')

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
        # If user provides a direct file, we should not auto-discover.
        if self.checkm2_file:
            return []

        if self.bin_dirs is None:
            expanded_dirs = [self.asmdir]
        else:
            # Use the same expanded directories logic as find_bins
            expanded_dirs = []
            for directory in self.bin_dirs:
                if '*' in directory:
                    expanded_dirs.extend(glob.glob(directory))
                else:
                    expanded_dirs.append(directory)
        
        for directory in expanded_dirs:
            # Run find command in each bin directory
            cmd = f"find {directory} -name 'quality_report.tsv'"
            try:
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                if result.stdout:
                    files = result.stdout.strip().split('\n')
                    checkm2_files.extend(files)
                    for file in files:
                        logging.info(f"Found CheckM2 file at {file}")
            except Exception as e:
                logging.error(f"Error running find command in {directory}: {str(e)}")
        
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
            try:
                self.checkm2_data = pd.read_csv(self.checkm2_file, sep='\t')
                logging.info(f"Successfully loaded CheckM2 data for {len(self.checkm2_data)} bins")
            except Exception as e:
                logging.error(f"Error loading CheckM2 file: {str(e)}")
                self.checkm2_data = None
            return

        checkm2_files = self.find_checkm2_outputs()
        if not checkm2_files:
            self.checkm2_data = None
            return

        # Combine all CheckM2 files
        dfs = []
        for file in checkm2_files:
            try:
                # Get the root directory (up to /bins)
                root_dir = os.path.dirname(os.path.dirname(file))
                logging.info(f"\nProcessing CheckM2 file: {file}")
                logging.info(f"Root directory: {root_dir}")
                
                # Load the file
                df = pd.read_csv(file, sep='\t')
                
                # Format the directory name to match symlink naming scheme
                if self.bin_dirs is None and self.asmdir:
                    asmdir_base = os.path.basename(os.path.normpath(self.asmdir)) or "asm"
                    rel_root = os.path.relpath(root_dir, start=self.asmdir)
                    dir_name = f"{asmdir_base}-{rel_root.replace('/', '-')}"
                else:
                    dir_name = root_dir.replace('/', '-')
                logging.info(f"Formatted directory name: {dir_name}")
                
                # Show original bin names
                logging.info("Original bin names in CheckM2 file:")
                for bin_name in df.iloc[:, 0]:
                    logging.info(f"  {bin_name}")
                
                # Append directory to bin names
                df.iloc[:, 0] = dir_name + '-' + df.iloc[:, 0].astype(str)
                
                # Show new bin names
                logging.info("New bin names after transformation:")
                for bin_name in df.iloc[:, 0]:
                    logging.info(f"  {bin_name}")
                
                dfs.append(df)
                
                # Write the transformed data to a file for inspection
                output_file = os.path.join(self.euk_binning_outputdir, f"checkm2_transformed_{os.path.basename(root_dir)}.tsv")
                df.to_csv(output_file, sep='\t', index=False)
                logging.info(f"Wrote transformed CheckM2 data to {output_file}")
                
            except Exception as e:
                logging.error(f"Error reading CheckM2 file {file}: {str(e)}")
        
        if dfs:
            self.checkm2_data = pd.concat(dfs, ignore_index=True)
            # Write the combined data to a file
            combined_file = os.path.join(self.euk_binning_outputdir, "checkm2_combined.tsv")
            self.checkm2_data.to_csv(combined_file, sep='\t', index=False)
            logging.info(f"Wrote combined CheckM2 data to {combined_file}")
            logging.info(f"Successfully loaded CheckM2 data for {len(self.checkm2_data)} bins")
        else:
            self.checkm2_data = None
            logging.info("No valid CheckM2 data found to load")

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
            summary_data.append(row)

        # Create summary DataFrame
        summary_df = pd.DataFrame(summary_data)
        logging.info("Created summary DataFrame")
        
        # Merge with CheckM2 data if available
        if self.checkm2_data is not None:
            logging.info("Merging CheckM2 data...")
            # Rename first column to bin_name and add checkm_ prefix to rest
            checkm2_cols = {col: f'checkm_{col}' for col in self.checkm2_data.columns if col != 'Name'}
            checkm2_df = self.checkm2_data.rename(columns={'Name': 'bin_name', **checkm2_cols})
            
            # Merge on bin name
            summary_df = pd.merge(summary_df, checkm2_df, on='bin_name', how='left')
            logging.info("CheckM2 data merged")
        
        # Save the final summary
        output_file = os.path.join(self.euk_binning_outputdir, 'eukaryotic_summary_table.csv')
        summary_df.to_csv(output_file, index=False)
        logging.info(f"Summary table saved to {output_file}")

    def run(self):
        if self.run_eukcc_flag:
            if not self.dblocs:
                raise ValueError("--dblocs is required when --run-eukcc is used")
            self.eukcc_db = self.get_db_location(self.dblocs, 'eukccdb')

        self.find_bins()
        self.create_bin_name_mapping()
        self.load_checkm2_data()
        if self.run_eukrep_flag:
            self.run_eukrep()
        if self.run_eukcc_flag:
            self.run_eukcc()
        self.process_output()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin_dirs", type=str, required=False, default=None, help="Pipe-separated list of directories containing bins, quoted. If omitted, MAGUS discovers bins from --asmdir.")
    parser.add_argument("--asmdir", type=str, default="asm", help="Assembly directory for default bin discovery (expects <asmdir>/*/bins/*.fa).")
    parser.add_argument("--wildcards", type=str, default = "", required=False, help="Pipe-separated list of patterns for bin files, quoted. Use 'bins' to search in subdirectories named 'bins', or file patterns like '.fa|.fasta'")
    parser.add_argument("--size_threshold", type=int, default=10000000)
    parser.add_argument("--euk-binning-outputdir", "--euk_binning_outputdir", dest="euk_binning_outputdir", type=str, default="magus_output/magus_euks")
    parser.add_argument("--dblocs", type=str, required=False, default=None, help="Path to database locations file (required when --run-eukcc is used)")
    parser.add_argument("--max-workers", dest="max_workers", type=int, default=1)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--run-eukrep", dest="run_eukrep", action='store_true', help="Run EukRep step (disabled by default)")
    parser.add_argument("--run-eukcc", dest="run_eukcc", action='store_true', help="Run EukCC step (disabled by default)")
    parser.add_argument("--eukrep_env", type=str, default="eukrep-env")
    parser.add_argument("--checkm2_file", type=str, help="Optional: Direct path to CheckM2 quality report file")
    args = parser.parse_args()

    bin_dirs = args.bin_dirs.split('|') if args.bin_dirs else None
    wildcards = args.wildcards.split('|')

    runner = EukRepRunner(
        bin_dirs=bin_dirs,
        asmdir=args.asmdir,
        wildcards=wildcards,
        size_threshold=args.size_threshold,
        euk_binning_outputdir=args.euk_binning_outputdir,
        dblocs=args.dblocs,
        max_workers=args.max_workers,
        threads=args.threads,
        run_eukrep=args.run_eukrep,
        run_eukcc=args.run_eukcc,
        eukrepenv=args.eukrep_env,
        checkm2_file=args.checkm2_file
    )
    runner.run()
