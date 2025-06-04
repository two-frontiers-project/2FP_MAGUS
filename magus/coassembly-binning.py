import subprocess
import os
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

class CoAssemblyBinning:
    def __init__(self, config, tmpdir, outdir, threads=28, checkm_db=None, max_workers=4, test_mode=False):
        self.config = self.load_config(config)
        self.outdir = outdir
        self.magdir = f"{self.outdir}/mags"
        self.tmpdir = tmpdir
        self.threads = threads
        self.checkm_db = checkm_db
        self.max_workers = max_workers
        self.test_mode = test_mode

        os.makedirs(self.outdir, exist_ok=True)
        os.makedirs(self.magdir, exist_ok=True)
        os.makedirs(self.tmpdir, exist_ok=True)

    def load_config(self, config_path):
        config_df = pd.read_csv(config_path, sep='\t')
        return config_df.to_dict(orient='records')


    def run_sorenson_alignment(self, BS):
        OUTF = f"{self.outdir}/{BS}"
        temp0_file = f"{OUTF}/temp0.fa"
        temp1_file = f"{OUTF}/temp1.fa"
        
        for sample in self.config:
            sample_name = sample['filename']
            r1 = sample['pe1']
            r2 = sample['pe2']
            print(f"Aligning {sample_name} to temp0.fa and temp1.fa for {BS}")
            
            subprocess.run(f"sorenson-g -db {temp0_file} -r1 {r1} -r2 {r2} -t {self.threads} -e .01 -o {OUTF}/CX-{sample_name}.cov", shell=True)
            subprocess.run(f"sorenson-g -db {temp1_file} -r1 {r1} -r2 {r2} -t {self.threads} -e .01 -o {OUTF}/CXb-{sample_name}.cov", shell=True)

    def merge_coverage_data(self, BS):
        OUTF = f"{self.outdir}/{BS}"

        print(f"Merging coverage data for {BS} temp0.fa")
        subprocess.run(f"cut -f1-3 {OUTF}/covs.txt > {OUTF}/combo.cov", shell=True)
        for cov_file in os.listdir(OUTF):
            if cov_file.startswith("CX-") and cov_file.endswith(".cov"):
                SN = cov_file.split('-')[1].replace(".cov", "")
                subprocess.run(f"cut -f4,5 {OUTF}/{cov_file} | sed 's/\tsample/\t{SN}/' | sed 's/^sample/{SN}/' | paste {OUTF}/combo.cov - > {OUTF}/temp.cov", shell=True)
                subprocess.run(f"mv {OUTF}/temp.cov {OUTF}/combo.cov", shell=True)
        subprocess.run(f"mv {OUTF}/combo.cov {OUTF}/combo1.cov", shell=True)

        print(f"Merging coverage data for {BS} temp1.fa")
        subprocess.run(f"cut -f1-3 {OUTF}/covs1.txt > {OUTF}/combo.cov", shell=True)
        for cov_file in os.listdir(OUTF):
            if cov_file.startswith("CXb-") and cov_file.endswith(".cov"):
                SN = cov_file.split('-')[1].replace(".cov", "")
                subprocess.run(f"cut -f4,5 {OUTF}/{cov_file} | sed 's/\tsample/\t{SN}/' | sed 's/^sample/{SN}/' | paste {OUTF}/combo.cov - > {OUTF}/temp.cov", shell=True)
                subprocess.run(f"mv {OUTF}/temp.cov {OUTF}/combo.cov", shell=True)
        subprocess.run(f"mv {OUTF}/combo.cov {OUTF}/combo1b.cov", shell=True)

    def run_metabat(self, BS):
        OUTF = f"{self.outdir}/{BS}"

        print(f"Running metabat2 on {BS} temp0.fa with combo1.cov")
        for k in range(15, 31):
            subprocess.run(f"metabat2 --seed {k} -i {OUTF}/temp0.fa -a {OUTF}/combo1.cov -t {self.threads} -o {OUTF}/bins/{k}CX1 -m {k*100}", shell=True)

        print(f"Running metabat2 on {BS} temp1.fa with combo1b.cov")
        for k in range(15, 31):
            subprocess.run(f"metabat2 --seed {k} -i {OUTF}/temp1.fa -a {OUTF}/combo1b.cov -t {self.threads} -o {OUTF}/bins/{k}CX1b -m {k*100}", shell=True)

    def run_checkm_binning(self, BS):
        OUTF = f"{self.outdir}/{BS}"
        checkm_output = f"{OUTF}/bins/temp"
        db_flag = f"--database_path {self.checkm_db}" if self.checkm_db else ""

        print(f"Running CheckM2 for {BS}")
        subprocess.run(f"checkm2 predict -x fa -i {OUTF}/bins/ -o {checkm_output} --force --remove_intermediates -t {self.threads} --tmpdir {self.tmpdir} {db_flag}", shell=True)

        # Set filtering criteria based on test_mode
        completeness_cutoff = 0 if self.test_mode else 50
        contamination_cutoff = 100 if self.test_mode else 5

        print(f"Filtering CheckM2 results for {BS} (completeness >= {completeness_cutoff}%, contamination <= {contamination_cutoff}%)")
        subprocess.run(f"awk -F'\\t' '$2 >= {completeness_cutoff} && $3 <= {contamination_cutoff}' {checkm_output}/quality_report.tsv > {OUTF}/worthwhile.tsv", shell=True)
        subprocess.run(f"mkdir -p {OUTF}/good && rm -f {OUTF}/good/*", shell=True)
        subprocess.run(f"cp -l $(cut -f1 {OUTF}/worthwhile.tsv | sed 's/$/.fa/' | sed 's:^:{OUTF}/bins/:') {OUTF}/good/", shell=True)

        print(f"Running lingenome for {BS}")
        subprocess.run(f"lingenome {OUTF}/good {OUTF}/good.fasta FILENAME", shell=True)

        print(f"Running akmer102 for {BS}")
        subprocess.run(f"OMP_NUM_THREADS={self.threads} akmer102 {OUTF}/good.fasta {OUTF}/good.dm 13 ANI CHANCE GC LOCAL RC", shell=True)

        print(f"Running spamw for {BS}")
        subprocess.run(f"spamw2 {OUTF}/good.dm {OUTF}/L 0 {self.threads} D2", shell=True)

        # Ensure worthwhile.stat is created correctly
        print(f"Generating worthwhile statistics for {BS}")
        subprocess.run(
            f"awk -F'\\t' '{{printf \"%s\\t%s\\t%s\\t%s\\t%.3f\\n\", $1, $4, $2, $3, $12/10}}' {OUTF}/worthwhile.tsv > {OUTF}/worthwhile.stat", 
            shell=True
        )

        # Check if worthwhile.tsv is empty
        if self.test_mode and os.path.getsize(f"{OUTF}/worthwhile.tsv") == 0:
            print(f"GENERATING SYNTHETIC MAGS FOR TEST MODE")
            temp0_path = os.path.abspath(f"{OUTF}/temp0.fa")
            bin_0_path = os.path.abspath(f"{self.magdir}/bin_0_coasm.fa")
            os.symlink(temp0_path, bin_0_path)

    def run_bestmag(self, BS):
        OUTF = f"{self.outdir}/{BS}"
        bestmags_txt = f"{OUTF}/bestmags.txt"

        print(f"Running bestmag for {BS}")
        subprocess.run(f"bestmag2 {OUTF}/good.dm {OUTF}/L.txt {OUTF}/worthwhile.stat {OUTF}/bestmags.txt SELF S1", shell=True)

        print(f"Moving best bins for {BS}")
        for z in subprocess.run(f"tail -n+2 {bestmags_txt} | cut -f1", shell=True, capture_output=True, text=True).stdout.splitlines():
            subprocess.run(f"cp -l {OUTF}/good/{z}.fa {self.magdir}/{BS}_{z}.fa", shell=True)
            subprocess.run(f"echo '{self.magdir}/{BS}_{z}.fa' >> {self.magdir}/checks_{os.getenv('NODE')}.txt", shell=True)

    def run(self):
        """ Run all steps in sequence with parallel execution. """
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {executor.submit(self.process_sample, BS): BS for BS in os.listdir(self.outdir) if os.path.isdir(f"{self.outdir}/{BS}") and BS != 'mags'}
            for future in as_completed(futures):
                BS = futures[future]
                try:
                    future.result()
                    print(f"Processing complete for {BS}")
                except Exception as exc:
                    print(f"Error processing {BS}: {exc}")

    def process_sample(self, BS):
        self.run_sorenson_alignment(BS)
        self.merge_coverage_data(BS)
        self.run_metabat(BS)
        self.run_checkm_binning(BS)
        self.run_bestmag(BS)

def main():
    parser = argparse.ArgumentParser(description="Run co-assembly and binning pipeline.")
    parser.add_argument('--config', type=str, required=True, help='Path to the configuration TSV file')
    parser.add_argument('--coasm_outdir', type=str, default="coasm", help='Output directory from previous step (default: coasm)')
    parser.add_argument('--tmp_dir', default = 'tmp/coassembly-binning',type=str, help='Temporary directory (default: tmp/coassembly-binning)')
    parser.add_argument('--threads', type=int, default=28, help='Number of threads (default: 48)')
    parser.add_argument('--checkm_db', type=str, help='Path to a custom CheckM database')
    parser.add_argument('--max_workers', type=int, default=4, help='Maximum number of parallel workers (default: 4)')
    parser.add_argument('--test_mode', action='store_true', help='Enable test mode to allow all bins regardless of quality (and generate spoof bins in the case of none being found)')

    args = parser.parse_args()

    binning = CoAssemblyBinning(
        config=args.config,
        outdir=args.coasm_outdir,
        tmpdir=args.tmp_dir,
        threads=args.threads,
        checkm_db=args.checkm_db,
        max_workers=args.max_workers,
        test_mode=args.test_mode
    )
    binning.run()

if __name__ == '__main__':
    main()