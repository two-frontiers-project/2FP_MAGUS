import subprocess
import os
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil

class Binning:
    def __init__(self, config,tmp_dir, asmdir, magdir="mags", threads=14, checkm_db=None, test_mode=False, max_workers=4,skipcm = False):
        self.config = self.load_config(config)
        self.asmdir = asmdir
        self.magdir = asmdir+"/"+magdir
        self.tmpdir = tmp_dir
        self.threads = threads
        self.checkm_db = checkm_db  # Custom CheckM database path
        self.test_mode = test_mode  # Flag for test mode
        self.max_workers = max_workers  # Maximum parallel workers
        self.skip_checkm = skipcm

        # Create the output directories if they don't exist
        os.makedirs(f"{self.magdir}", exist_ok=True)
        os.makedirs(self.tmpdir, exist_ok=True)

    def load_config(self, config_path):
        config_df = pd.read_csv(config_path, sep='\t')
        if 'pe2' in config_df.columns:
            return config_df.to_dict(orient='records')
        else:
            return config_df.assign(pe2=None).to_dict(orient='records')

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
        
        if self.test_mode:
            completeness_cutoff = 0
            contamination_cutoff = 100
        else:
            completeness_cutoff = 50
            contamination_cutoff = 5

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

    def run_akmer100b(self, sample_name):
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
        """Runs the entire binning pipeline for a single sample."""
        try:
            #self.run_sorenson(sample_name)
            #self.run_metabat(sample_name)
            if self.skip_checkm:
                return 'Stopping after binning.'
            self.run_checkm(sample_name)
            self.filter_good_bins(sample_name)
            self.copy_good_bins(sample_name)
            self.run_lingenome(sample_name)
            self.run_akmer100b(sample_name)
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
    parser.add_argument('--max_workers', type=int, default=4, help='Maximum number of parallel workers (default: 4)')
    parser.add_argument('--tmp_dir', type=str, default='tmp/single-binning', help='Temp directory. Default tmp/single-binning.')
    parser.add_argument('--test_mode', action='store_true', help='Enable test mode with relaxed filtering criteria')
    parser.add_argument('--skipcheckm',action = 'store_true', help='Skips checkm and downstream steps (so stops immediately after binning)')

    # Parse arguments
    args = parser.parse_args()

    # Initialize and run the Binning instance
    binning = Binning(
        config=args.config,
        threads=args.threads,
        checkm_db=args.checkm_db,
        asmdir=args.asmdir,
        tmp_dir=args.tmp_dir,
        max_workers=args.max_workers,
        test_mode=args.test_mode,
        skipcm=args.skipcheckm
    )
    binning.run()

if __name__ == '__main__':
    main()
