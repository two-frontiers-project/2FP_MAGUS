import subprocess
import os

class CoAssembly:
    def __init__(self, config, outdir="coasm", magdir=None, tmpdir=None, threads=48):
        self.config = config
        self.outdir = outdir
        self.magdir = magdir if magdir else f"{self.outdir}/mags"
        self.tmpdir = tmpdir if tmpdir else os.path.expanduser("~/partmp")
        self.threads = threads
        
        # Create necessary directories
        os.makedirs(self.outdir, exist_ok=True)
        os.makedirs(self.magdir, exist_ok=True)
        os.makedirs(self.tmpdir, exist_ok=True)

    def run_coassembly(self):
        coasm_todo = "coasm.todo"
        
        with open(coasm_todo) as f:
            for line in f:
                # Extract the base sample name and file list
                x = line.strip().replace('\t', ';')
                BS = x.split(';')[0]
                OUTF = f"{self.outdir}/{BS}"
                R1 = f"{OUTF}/{BS}.merged.R1.gz"
                R2 = f"{OUTF}/{BS}.merged.R2.gz"

                # Create output directory for the base sample
                os.makedirs(OUTF, exist_ok=True)
                
                # Remove any previous R1/R2 files
                if os.path.exists(R1):
                    os.remove(R1)
                if os.path.exists(R2):
                    os.remove(R2)

                # Concatenate all R1 and R2 files from the qc directory into a single R1 and R2 file for co-assembly
                for s in x.split(';'):
                    # Find the sample in the config
                    sample_config = next((item for item in self.config if item['filename'] == s), None)
                    if sample_config:
                        # Concatenate the appropriate files from the qc directory into one R1 and one R2
                        r1_cat = f"cat {sample_config['pe1']} >> {R1}"
                        r2_cat = f"cat {sample_config['pe2']} >> {R2}"
                        print(f"Concatenating {sample_config['pe1']} to {R1}")
                        subprocess.run(r1_cat, shell=True)
                        print(f"Concatenating {sample_config['pe2']} to {R2}")
                        subprocess.run(r2_cat, shell=True)

                # Run megahit for co-assembly
                megahit_cmd = (
                    f"megahit -f --k-min 75 --k-max 333 --k-step 6 "
                    f"--cleaning-rounds 1 --merge-level 100,.999 "
                    f"--min-count 1 --min-contig-len 1000 --continue "
                    f"-t {self.threads} -1 {R1} -2 {R2} -o {OUTF}"
                )
                print(f"Running megahit for {BS}")
                subprocess.run(megahit_cmd, shell=True)

                # Move intermediate contigs and clean up
                mv_cmd = f"mv {OUTF}/intermediate_contigs/k147.contigs.fa {OUTF}/"
                rm_cmd = f"rm -r {OUTF}/intermediate_contigs/"
                print(f"Moving contigs and cleaning up for {BS}")
                subprocess.run(mv_cmd, shell=True)
                subprocess.run(rm_cmd, shell=True)

    def run_filtering_binning(self):
        """ Perform stronger filtering by depth, and lighter short-contig filtering, followed by binning. """
        coasm_dirs = [d for d in os.listdir(self.outdir) if os.path.isdir(f"{self.outdir}/{d}")]
        
        for BS in coasm_dirs:  # Loop through coassembly directories only
            OUTF = f"{self.outdir}/{BS}"
            R1 = f"{OUTF}/{BS}.merged.R1.gz"  # Use the merged reads from coassembly
            R2 = f"{OUTF}/{BS}.merged.R2.gz"

            # Ensure merged reads exist
            if not (os.path.exists(R1) and os.path.exists(R2)):
                print(f"Merged reads for {BS} not found. Skipping...")
                continue

            # Step 1: Stronger filtering by depth
            print(f"Running stronger filtering for {BS}")
            subprocess.run(f"rm -rf {OUTF}/bins/*", shell=True)  # Ensure bins directory is clean
            subprocess.run(f"grep -A1 --no-group-separator 'multi=[2-9]\\|multi=[0-9][0-9]' {OUTF}/final.contigs.fa | "
                           f"awk '{{if (0==(NR % 2) && length >= 1000) {{print x; print $0}}; x=$0}}' > {OUTF}/temp0.fa", shell=True)
            subprocess.run(f"sorenson-g -db {OUTF}/temp0.fa -qc -r1 {R1} -r2 {R2} -t {self.threads} -e 0.01 -o {OUTF}/covs.txt", shell=True)

            # Step 2: Run metabat2 for binning
            print(f"Running metabat2 for binning {BS}")
            for k in range(15, 31):
                subprocess.run(f"metabat2 --seed {k} -i {OUTF}/temp0.fa -a {OUTF}/covs.txt -t {self.threads} -o {OUTF}/bins/{k}CX0 -m {k*100}", shell=True)

            # Step 3: Lighter filtering for shorter contigs
            print(f"Running lighter filtering for shorter contigs for {BS}")
            subprocess.run(f"grep -A1 --no-group-separator 'multi=[2-9]\\|multi=[0-9][0-9]' {OUTF}/k147.contigs.fa | "
                           f"awk '{{if (0==(NR % 2) && length >= 500) {{print x; print $0}}; x=$0}}' > {OUTF}/temp1.fa", shell=True)
            subprocess.run(f"sorenson-g -db {OUTF}/temp1.fa -qc -r1 {R1} -r2 {R2} -t {self.threads} -e 0.01 -o {OUTF}/covs1.txt", shell=True)

            # Step 4: Run metabat2 again for shorter contigs
            print(f"Running second metabat2 pass for {BS}")
            for k in range(15, 31):
                subprocess.run(f"metabat2 --seed {k+100} -i {OUTF}/temp1.fa -a {OUTF}/covs1.txt -t {self.threads} -o {OUTF}/bins/{k}CX0b -m {k*100}", shell=True)

    def run(self):
        #self.run_coassembly()
        self.run_filtering_binning()


# Example usage with config:
config = [
    {
        'filename': 'SRR15373729', 
        'pe1': './qc/SRR15373729.R1.fa.gz', 
        'pe2': './qc/SRR15373729.R2.fa.gz'
    },
    {
        'filename': 'SRR15373733', 
        'pe1': './qc/SRR15373733.R1.fa.gz', 
        'pe2': './qc/SRR15373733.R2.fa.gz'
    }
]

coassembly = CoAssembly(config)
coassembly.run()
