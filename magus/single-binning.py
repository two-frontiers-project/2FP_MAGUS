import subprocess
import os
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil

class Binning:
    def __init__(self, config, tmpdir, asmdir, magdir="mags", threads=14, checkm_db=None,
                 test_mode=False, max_workers=4, completeness=None,
                 contamination=None, quality=None, restart=None):
        self.config = self.load_config(config)
        self.asmdir = asmdir
        self.magdir = asmdir+"/"+magdir
        self.tmpdir = tmpdir
        self.threads = threads
        self.checkm_db = checkm_db  # Custom CheckM database path
        self.test_mode = test_mode  # Flag for test mode
        self.max_workers = max_workers  # Maximum parallel workers
        self.completeness = completeness
        self.contamination = contamination
        self.quality = quality
        # restart can be 'binning', 'checkm', or 'filtering' to skip completed stages
        self.restart = restart

        # Create the output directories if they don't exist
        os.makedirs(f"{self.magdir}", exist_ok=True)
        os.makedirs(self.tmpdir, exist_ok=True)

    def load_config(self, config_path):
        config_df = pd.read_csv(config_path, sep='\t')
        if 'pe2' in config_df.columns:
            return config_df.to_dict(orient='records')
        else:
            return config_df.assign(pe2=None).to_dict(orient='records')

    def recreate_tmp_structure_for_restart(self, sample_name):
        """
        Recreate the tmp directory structure for restart by symlinking to preserved data 
        in the assembly directory.
        """
        sample_tmp_dir = f"{self.tmpdir}/{sample_name}"
        sample_asm_dir = f"{self.asmdir}/{sample_name}"
        
        # Check if assembly directory exists (where data was moved to)
        if not os.path.exists(sample_asm_dir):
            print(f"Warning: Assembly directory {sample_asm_dir} does not exist. Cannot recreate tmp structure for restart.")
            return False
        
        # Create the sample tmp directory if it doesn't exist
        os.makedirs(sample_tmp_dir, exist_ok=True)
        
        def safe_symlink(source, target, description):
            """Helper function to safely create symlinks, removing broken ones first."""
            if os.path.islink(target) and not os.path.exists(target):
                # Remove broken symlink
                os.unlink(target)
                print(f"Removed broken symlink: {target}")
            
            if os.path.exists(source) and not os.path.exists(target):
                os.symlink(os.path.abspath(source), target)
                print(f"Created symlink for {description}: {target} -> {source}")
                return True
            return False
        
        # Determine what needs to be recreated based on restart mode
        if self.restart == "binning":
            # Need covs.txt for metabat
            covs_source = f"{sample_asm_dir}/covs.txt"
            covs_target = f"{sample_tmp_dir}/covs.txt"
            safe_symlink(covs_source, covs_target, "covs.txt")
                
        elif self.restart == "checkm":
            # Need bins directory for checkm
            bins_source = f"{sample_asm_dir}/bins"
            bins_target = f"{sample_tmp_dir}/bins"
            safe_symlink(bins_source, bins_target, "bins directory")
                
        elif self.restart == "filtering":
            # Need bins directory with checkm results for filtering
            bins_source = f"{sample_asm_dir}/bins"
            bins_target = f"{sample_tmp_dir}/bins"
            safe_symlink(bins_source, bins_target, "bins directory")
            
            # For filtering restarts, clean up downstream files that depend on filtering thresholds
            # This ensures a fresh restart when using different contamination/completeness cutoffs
            downstream_files = [
                "worthwhile.tsv", "worthwhile.stat", "good.fasta", "good.dm", 
                "L", "L.txt", "bestmags.txt"
            ]
            
            for file_name in downstream_files:
                target_file = f"{sample_tmp_dir}/{file_name}"
                if os.path.exists(target_file) or os.path.islink(target_file):
                    os.unlink(target_file)
                    print(f"Removed existing {file_name} for clean filtering restart")
            
            # Also remove the good directory if it exists (will be recreated with new thresholds)
            good_target = f"{sample_tmp_dir}/good"
            if os.path.exists(good_target) or os.path.islink(good_target):
                if os.path.islink(good_target):
                    os.unlink(good_target)
                else:
                    shutil.rmtree(good_target)
                print(f"Removed existing good directory for clean filtering restart")
        
        return True

    def run_sorenson(self, sample_name):
        out_sample_dir = f"{self.tmpdir}/{sample_name}"
        os.makedirs(out_sample_dir, exist_ok=True)
        temp0_file = f"{self.asmdir}/{sample_name}/temp0.fa"
        sample = next(item for item in self.config if item["filename"] == sample_name)
        r1 = sample["pe1"]
        r2 = sample.get("pe2", None)
        cov_file = f"{out_sample_dir}/covs.txt"
        
        if pd.isna(r2) or r2 is None:
            cmd = f"sorenson-g -db {temp0_file} -qc -r1 {r1} -t {self.threads} -e 0.01 -o {cov_file}"
            print(f"Running sorenson-g in single-end mode on {sample_name}")
        else:
            cmd = f"sorenson-g -db {temp0_file} -qc -r1 {r1} -r2 {r2} -t {self.threads} -e 0.01 -o {cov_file}"
            print(f"Running sorenson-g in paired-end mode on {sample_name}")
        print(cmd) 
        subprocess.run(cmd, shell=True)

    def run_metabat(self, sample_name):
        out_sample_dir = f"{self.tmpdir}/{sample_name}"
        temp0_file = f"{self.asmdir}/{sample_name}/temp0.fa"
        cov_file = f"{out_sample_dir}/covs.txt"
        bins_dir = f"{out_sample_dir}/bins"
        
        # Ensure bins directory exists
        os.makedirs(bins_dir, exist_ok=True)
        
        # Run MetaBAT2 with different seed values
        for k in [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 27, 30]:
            bin_output = f"{bins_dir}/{k}"
            cmd = f"metabat2 --seed {k} -i {temp0_file} -a {cov_file} -t {self.threads} -o {bin_output} -m 1500"
            print(cmd)
            print(f"Running MetaBAT2 with seed {k} for {sample_name}")
            subprocess.run(cmd, shell=True)
        
        # Check if any bins were created
        bins_found = any(os.path.isfile(os.path.join(bins_dir, file)) for file in os.listdir(bins_dir))
        
        # If no bins were found and test_mode is enabled, create a fallback bin
        if self.test_mode and not bins_found:
            fallback_bin = os.path.join(bins_dir, "bin0.fa")
            final_contigs_path = f"{self.asmdir}/{sample_name}/final.contigs.fa"
            print(f"IN TEST MODE, GENERATING SYNTHETIC BIN. Copying {final_contigs_path} to {fallback_bin} as fallback bin.")
            shutil.copy(final_contigs_path, fallback_bin)

    def run_checkm(self, sample_name):
        out_sample_dir = f"{self.tmpdir}/{sample_name}"
        bins_dir = f"{out_sample_dir}/bins"
        checkm_output = f"{bins_dir}/checkm"
        
        if self.checkm_db:
            db_flag = f"--database_path {self.checkm_db}"
        else:
            db_flag = ""  # Default database path

        cmd = f"checkm2 predict -x fa -i {bins_dir} -o {checkm_output} --force --remove_intermediates -t {self.threads} --tmpdir {self.tmpdir} {db_flag}"
        print(f"Running CheckM for {sample_name} with database: {self.checkm_db if self.checkm_db else 'default'}")
        subprocess.run(cmd, shell=True)

    def filter_good_bins(self, sample_name):
        out_sample_dir = f"{self.tmpdir}/{sample_name}"
        checkm_report = f"{out_sample_dir}/bins/checkm/quality_report.tsv"
        worthwhile_bins = f"{out_sample_dir}/worthwhile.tsv"
        
        # Determine filtering thresholds
        if self.quality == 'low':
            completeness_cutoff = 0
            contamination_cutoff = 100
        elif self.quality == 'medium':
            completeness_cutoff = 50
            contamination_cutoff = 10
        elif self.quality == 'high':
            completeness_cutoff = 90
            contamination_cutoff = 5
        else:
            completeness_cutoff = 50
            contamination_cutoff = 5

        if self.completeness is not None:
            completeness_cutoff = self.completeness
        if self.contamination is not None:
            contamination_cutoff = self.contamination

        if self.test_mode:
            completeness_cutoff = 0
            contamination_cutoff = 100

        cmd = f"awk -F'\\t' '$2 >= {completeness_cutoff} && $3 <= {contamination_cutoff}' {checkm_report} > {worthwhile_bins}"
        subprocess.run(cmd, shell=True)
        print(f"Filtered good bins for {sample_name} (completeness >= {completeness_cutoff}%, contamination <= {contamination_cutoff}%)")

    def copy_good_bins(self, sample_name):
        out_sample_dir = f"{self.tmpdir}/{sample_name}"
        worthwhile_bins = f"{out_sample_dir}/worthwhile.tsv"
        bins_dir = f"{out_sample_dir}/bins"
        good_dir = f"{out_sample_dir}/good"
        
        os.makedirs(good_dir, exist_ok=True)

        cmd = f"cut -f1 {worthwhile_bins} | sed 's/$/.fa/' | sed 's:^:{bins_dir}/:'"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        bin_files = result.stdout.strip().split('\n')

        if not bin_files or bin_files == ['']:
            print(f"No bins to copy for {sample_name}.")
            return

        print(f"Files to copy for {sample_name}: {bin_files}")
        
        for bin_file in bin_files:
            cmd_copy = f"cp {bin_file} {good_dir}/"
            subprocess.run(cmd_copy, shell=True)
            print(f"Copied {bin_file} to {good_dir}")

    def run_lingenome(self, sample_name):
        out_sample_dir = f"{self.tmpdir}/{sample_name}/good"
        parent_dir = f"{self.tmpdir}/{sample_name}"
        good_fasta = f"{parent_dir}/good.fasta"
        cmd = f"lingenome {out_sample_dir} {good_fasta} FILENAME"
        print(f"Running lingenome for {sample_name}")
        subprocess.run(cmd, shell=True)

    def run_akmer102(self, sample_name):
        parent_dir = f"{self.tmpdir}/{sample_name}"
        good_fasta = f"{parent_dir}/good.fasta"
        good_dm = f"{parent_dir}/good.dm"
        cmd = f"OMP_NUM_THREADS={self.threads} akmer102 {good_fasta} {good_dm} 13 ANI CHANCE GC LOCAL RC"
        print(f"Running akmer102 for {sample_name}")
        subprocess.run(cmd, shell=True)

    def run_spamw(self, sample_name):
        parent_dir = f"{self.tmpdir}/{sample_name}"
        good_dm = f"{parent_dir}/good.dm"
        L_output = f"{parent_dir}/L"
        cmd = f"spamw2 {good_dm} {L_output} 0 {self.threads} D2"
        print(f"Running spamw2 for {sample_name}")
        subprocess.run(cmd, shell=True)
        worthwhile_stat = f"{parent_dir}/worthwhile.stat"
        worthwhile_tsv = f"{parent_dir}/worthwhile.tsv"
        cmd_process_bins = (
            f"awk -F'\\t' '{{printf \"%s\\t%s\\t%s\\t%s\\t%.3f\\n\", $1, $4, $2, $3, $12/10}}' {worthwhile_tsv} > {worthwhile_stat}"
        )
        print(f"Processing worthwhile bins for {sample_name}")
        subprocess.run(cmd_process_bins, shell=True)

    def run_bestmag(self, sample_name):
        parent_dir = f"{self.tmpdir}/{sample_name}"
        worthwhile_stat = f"{parent_dir}/worthwhile.stat"
        L_output = f"{parent_dir}/L.txt"
        bestmags_txt = f"{parent_dir}/bestmags.txt"
        cmd = f"bestmag2 {parent_dir}/good.dm {L_output} {worthwhile_stat} {bestmags_txt} SELF S1"
        print(f"Running bestmag for {sample_name}")
        subprocess.run(cmd, shell=True)

    def run_final_copy(self, sample_name):
        parent_dir = f"{self.tmpdir}/{sample_name}"
        mags_dir = self.magdir
        bestmags_txt = f"{parent_dir}/bestmags.txt"
        for z in subprocess.run(f"tail -n+2 {bestmags_txt} | cut -f1", shell=True, capture_output=True, text=True).stdout.splitlines():
            print(z)
            cmd = f"cp {parent_dir}/bins/{z}.fa {mags_dir}/{sample_name}_{z}.fa"
            subprocess.run(cmd, shell=True)
            print(f"Copied final bin {z} for {sample_name} to {mags_dir}")
        checks_file = f"{mags_dir}/checks_single_assembly.txt"
        cmd_append = f"tail -n+2 {bestmags_txt} | sed 's/^/{sample_name}_/' >> {checks_file}; mv {self.tmpdir}/{sample_name}/* {self.asmdir}/{sample_name}"
        subprocess.run(cmd_append, shell=True)

    def run_sample(self, sample_name):
        """Run the pipeline for one sample starting from the chosen step."""
        try:
            # If restarting, check and recreate tmp structure if needed
            if self.restart is not None:
                success = self.recreate_tmp_structure_for_restart(sample_name)
                if not success:
                    print(f"Failed to recreate tmp structure for {sample_name}. Cannot restart.")
                    return
            
            if self.restart is None:
                self.run_sorenson(sample_name)
                self.run_metabat(sample_name)
                self.run_checkm(sample_name)
            elif self.restart == "binning":
                self.run_metabat(sample_name)
                self.run_checkm(sample_name)
            elif self.restart == "checkm":
                self.run_checkm(sample_name)
            # restart == "filtering" skips directly to filtering steps
            self.filter_good_bins(sample_name)
            self.copy_good_bins(sample_name)
            self.run_lingenome(sample_name)
            self.run_akmer102(sample_name)
            self.run_spamw(sample_name)
            self.run_bestmag(sample_name)
            self.run_final_copy(sample_name)
            print(f"Completed processing for {sample_name}")
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while processing {sample_name}: {e}")

    def run(self):
        """Runs the binning pipeline for all samples in parallel."""
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {executor.submit(self.run_sample, sample['filename']): sample for sample in self.config}
            for future in as_completed(futures):
                sample = futures[future]
                sample_name = sample['filename']
                try:
                    future.result()  # Retrieving result or exception
                except Exception as exc:
                    print(f"Sample {sample_name} generated an exception: {exc}")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Run binning pipeline on genomic data.")
    parser.add_argument('--config', type=str, required=True, help='Path to the configuration TSV file')
    parser.add_argument('--threads', type=int, default=14, help='Number of threads for tools (default: 14)')
    parser.add_argument('--checkm_db', type=str, default=None, help='Path to custom CheckM database')
    parser.add_argument('--asmdir', type=str, default="asm", help='Output directory for assembly (default: asm)')
    parser.add_argument('--max-workers', '--max_workers', dest='max_workers', type=int, default=4, help='Maximum number of parallel workers (default: 4)')
    parser.add_argument('--tmpdir', type=str, default='tmp/binning', help='Temp directory. Default tmp/binning.')
    parser.add_argument('--test_mode', action='store_true', help='Enable test mode with relaxed filtering criteria')
    parser.add_argument('--completeness', type=float, default=None, help='Completeness threshold for filtering bins')
    parser.add_argument('--contamination', type=float, default=None, help='Contamination threshold for filtering bins')
    parser.add_argument('--restart', choices=['binning','checkm','filtering'], default=None,
                        help='Resume pipeline from the specified step')
    quality_group = parser.add_mutually_exclusive_group()
    quality_group.add_argument('--low-quality', dest='low_quality', action='store_true', help='Allow low quality bins (not recommended)')
    quality_group.add_argument('--medium-quality', dest='medium_quality', action='store_true', help='Use medium quality thresholds (completeness >=50, contamination <=10)')
    quality_group.add_argument('--high-quality', dest='high_quality', action='store_true', help='Use high quality thresholds (completeness >=90, contamination <=5)')

    # Parse arguments
    args = parser.parse_args()

    # Initialize and run the Binning instance
    quality = None
    if args.low_quality:
        quality = 'low'
    elif args.medium_quality:
        quality = 'medium'
    elif args.high_quality:
        quality = 'high'

    binning = Binning(
        config=args.config,
        threads=args.threads,
        checkm_db=args.checkm_db,
        asmdir=args.asmdir,
        tmpdir=args.tmpdir,
        max_workers=args.max_workers,
        test_mode=args.test_mode,
        completeness=args.completeness,
        contamination=args.contamination,
        quality=quality,
        restart=args.restart
    )
    binning.run()

if __name__ == '__main__':
    main()
