import subprocess
import os
import argparse
import pandas as pd
import shutil

class CoAssembly:
    def __init__(self, config, tmpdir, coasm_todo,outdir="coasm", threads=48, test_mode=False):
        self.config = self.load_config(config)
        self.outdir = outdir
        self.magdir = f"{self.outdir}/mags"
        self.tmpdir = tmpdir
        self.threads = threads
        self.test_mode = test_mode  # Added test_mode flag
        self.coasm_todo = coasm_todo

        # Create necessary directories
        os.makedirs(self.outdir, exist_ok=True)
        os.makedirs(self.magdir, exist_ok=True)
        os.makedirs(self.tmpdir, exist_ok=True)

    def load_config(self, config_path):
        config_df = pd.read_csv(config_path, sep='\t')
        return config_df.to_dict(orient='records')

    def run_coassembly(self):
        coasm_todo = self.coasm_todo
        
        with open(coasm_todo) as f:
            for line in f:
                x = line.strip().replace('\t', ';')
                BS = x.split(';')[0]
                OUTF = f"{self.outdir}/{BS}"
                R1 = f"{OUTF}/{BS}.merged.R1.gz"
                R2 = f"{OUTF}/{BS}.merged.R2.gz"

                os.makedirs(OUTF, exist_ok=True)
                
                if os.path.exists(R1):
                    os.remove(R1)
                if os.path.exists(R2):
                    os.remove(R2)

                for s in x.split(';'):
                    sample_config = next((item for item in self.config if item['filename'] == s), None)
                    if sample_config:
                        r1_cat = f"cat {sample_config['pe1']} >> {R1}"
                        subprocess.run(r1_cat, shell=True)
                        print(f"Concatenating {sample_config['pe1']} to {R1}")
                        
                        r2 = sample_config.get('pe2', None)
                        if pd.isna(r2) or r2 is None:
                            R2 = None
                        else:
                            r2_cat = f"cat {r2} >> {R2}"
                            subprocess.run(r2_cat, shell=True)
                            print(f"Concatenating {r2} to {R2}")

                if R2 is None:
                    megahit_cmd = (
                        f"megahit-g -f --k-min 75 --k-max 333 --k-step 6 "
                        f"--cleaning-rounds 1 --merge-level 100,.999 "
                        f"--min-count 1 --min-contig-len 1000 --continue "
                        f"-t {self.threads} -r {R1} -o {OUTF}"
                    )
                    print(f"Running megahit in single-end mode for {BS}")
                else:
                    megahit_cmd = (
                        f"megahit-g -f --k-min 75 --k-max 333 --k-step 6 "
                        f"--cleaning-rounds 1 --merge-level 100,.999 "
                        f"--min-count 1 --min-contig-len 1000 --continue "
                        f"-t {self.threads} -1 {R1} -2 {R2} -o {OUTF}"
                    )
                    print(f"Running megahit in paired-end mode for {BS}")
                
                subprocess.run(megahit_cmd, shell=True)

                # Rename headers in final.contigs.fa with 'coassembly_all'
                final_contigs_file = f"{OUTF}/final.contigs.fa"
                renamed_final_contigs_file = f"{OUTF}/final_renamed.contigs.fa"
                rename_final_cmd = f"sed 's/^>/>{BS}_coassembly_all_/' {final_contigs_file} > {renamed_final_contigs_file}"
                subprocess.run(rename_final_cmd, shell=True, check=True)
                subprocess.run(f"mv {renamed_final_contigs_file} {final_contigs_file}", shell=True)

                # Move and rename k147.contigs.fa with 'coassembly_k147'
                mv_cmd = f"mv {OUTF}/intermediate_contigs/k147.contigs.fa {OUTF}/k147.contigs.fa"
                print(f"Moving k147 contigs and cleaning up for {BS}")
                subprocess.run(mv_cmd, shell=True)

                renamed_k147_contigs_file = f"{OUTF}/k147_renamed.contigs.fa"
                rename_k147_cmd = f"sed 's/^>/>{BS}_coassembly_k147_/' {OUTF}/k147.contigs.fa > {renamed_k147_contigs_file}"
                subprocess.run(rename_k147_cmd, shell=True, check=True)
                subprocess.run(f"mv {renamed_k147_contigs_file} {OUTF}/k147.contigs.fa", shell=True)

                # Cleanup intermediate directory
                rm_cmd = f"rm -r {OUTF}/intermediate_contigs/"
                subprocess.run(rm_cmd, shell=True)

    def run_filtering_binning(self):
        coasm_dirs = [d for d in os.listdir(self.outdir) if os.path.isdir(f"{self.outdir}/{d}")]
        coasm_dirs = [d for d in coasm_dirs if d!="mags"]

        for BS in coasm_dirs:
            OUTF = f"{self.outdir}/{BS}"
            R1 = f"{OUTF}/{BS}.merged.R1.gz"
            R2 = f"{OUTF}/{BS}.merged.R2.gz"
            bins_dir = f"{OUTF}/bins"

            if not os.path.exists(R1):
                print(f"Merged reads for {BS} not found. Skipping...")
                continue
            
            os.makedirs(bins_dir, exist_ok=True)

            print(f"Running stronger filtering for {BS}")
            subprocess.run(f"rm -rf {bins_dir}/*", shell=True)
            subprocess.run(f"grep -A1 --no-group-separator 'multi=[2-9]\\|multi=[0-9][0-9]' {OUTF}/final.contigs.fa | "
                           f"awk '{{if (0==(NR % 2) && length >= 1000) {{print x; print $0}}; x=$0}}' > {OUTF}/temp0.fa", shell=True)
            
            if os.path.exists(R2):
                sorenson_cmd = f"sorenson-g -db {OUTF}/temp0.fa -qc -r1 {R1} -r2 {R2} -t {self.threads} -e 0.01 -o {OUTF}/covs.txt"
            else:
                sorenson_cmd = f"sorenson-g -db {OUTF}/temp0.fa -qc -r1 {R1} -t {self.threads} -e 0.01 -o {OUTF}/covs.txt"
            
            subprocess.run(sorenson_cmd, shell=True)

            print(f"Running metabat2 for binning {BS}")
            for k in range(15, 31):
                subprocess.run(f"metabat2 --seed {k} -i {OUTF}/temp0.fa -a {OUTF}/covs.txt -t {self.threads} -o {bins_dir}/{k}CX0 -m {k*100}", shell=True)

            bins_found = any(os.path.isfile(os.path.join(bins_dir, file)) for file in os.listdir(bins_dir))
            
            if self.test_mode and not bins_found:
                fallback_bin = os.path.join(bins_dir, "bin0.fa")
                final_contigs_path = f"{OUTF}/final.contigs.fa"
                print(f"No bins found for {BS} in test mode. Copying {final_contigs_path} to {fallback_bin} as fallback bin.")
                shutil.copy(final_contigs_path, fallback_bin)

            print(f"Running lighter filtering for shorter contigs for {BS}")
            subprocess.run(f"grep -A1 --no-group-separator 'multi=[2-9]\\|multi=[0-9][0-9]' {OUTF}/k147.contigs.fa | "
                           f"awk '{{if (0==(NR % 2) && length >= 500) {{print x; print $0}}; x=$0}}' > {OUTF}/temp1.fa", shell=True)
            
            if os.path.exists(R2):
                sorenson_cmd2 = f"sorenson-g -db {OUTF}/temp1.fa -qc -r1 {R1} -r2 {R2} -t {self.threads} -e 0.01 -o {OUTF}/covs1.txt"
            else:
                sorenson_cmd2 = f"sorenson-g -db {OUTF}/temp1.fa -qc -r1 {R1} -t {self.threads} -e 0.01 -o {OUTF}/covs1.txt"
            
            subprocess.run(sorenson_cmd2, shell=True)

            print(f"Running second metabat2 pass for {BS}")
            for k in range(15, 31):
                subprocess.run(f"metabat2 --seed {k+100} -i {OUTF}/temp1.fa -a {OUTF}/covs1.txt -t {self.threads} -o {bins_dir}/{k}CX0b -m {k*100}", shell=True)


    def run(self):
        self.run_coassembly()
        self.run_filtering_binning()

def main():
    parser = argparse.ArgumentParser(description="Run co-assembly pipeline on genomic data.")
    parser.add_argument('--config', type=str, required=True, help='Path to the configuration TSV file')
    parser.add_argument('--coasm_todo', type=str, required=True, help='Path to output of cluster_contigs, the coasm_todo file')
    parser.add_argument('--outdir', type=str, default="coasm", help='Output directory for co-assembly (default: coasm)')
    parser.add_argument(
        '--tmpdir',
        type=str,
        default='tmp/coasm',
        help='Temporary directory (default: tmp/coasm)'
    )
    parser.add_argument('--threads', type=int, default=48, help='Number of threads for tools (default: 48)')
    parser.add_argument('--test_mode', action='store_true', help='Enable test mode with relaxed filtering criteria')

    args = parser.parse_args()

    coassembly = CoAssembly(
        config=args.config,
        outdir=args.outdir,
        tmpdir=args.tmpdir,
        coasm_todo=args.coasm_todo,
        threads=args.threads,
        test_mode=args.test_mode
    )
    coassembly.run()

if __name__ == '__main__':
    main()
