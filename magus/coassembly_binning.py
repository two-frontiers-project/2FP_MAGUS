import subprocess
import os

class CoAssemblyBinning:
    def __init__(self, config, outdir="coasm", magdir=None, tmpdir=None, threads=48, checkm_db=None):
        self.config = config
        self.outdir = outdir
        self.magdir = magdir if magdir else f"{self.outdir}/mags"
        self.tmpdir = tmpdir if tmpdir else os.path.expanduser("~/partmp")
        self.threads = threads
        self.checkm_db = checkm_db  # Optional custom CheckM database path
        
        # Create necessary directories
        os.makedirs(self.outdir, exist_ok=True)
        os.makedirs(self.magdir, exist_ok=True)
        os.makedirs(self.tmpdir, exist_ok=True)

    def run_sorenson_alignment(self, BS):
        """ Align constituent samples to filtered contigs without QC. """
        OUTF = f"{self.outdir}/{BS}"
        temp0_file = f"{OUTF}/temp0.fa"
        temp1_file = f"{OUTF}/temp1.fa"
        
        # Align each constituent sample (from the original config) to temp0.fa and temp1.fa
        for sample in self.config:
            sample_name = sample['filename']
            r1 = sample['pe1']
            r2 = sample['pe2']
            print(f"Aligning {sample_name} to temp0.fa and temp1.fa for {BS}")
            
            # Align temp0
            subprocess.run(f"sorenson-g -db {temp0_file} -r1 {r1} -r2 {r2} -t {self.threads} -e .01 -o {OUTF}/CX-{sample_name}.cov", shell=True)
            # Align temp1
            subprocess.run(f"sorenson-g -db {temp1_file} -r1 {r1} -r2 {r2} -t {self.threads} -e .01 -o {OUTF}/CXb-{sample_name}.cov", shell=True)

    def merge_coverage_data(self, BS):
        """ Merge coverage data from all aligned samples. """
        OUTF = f"{self.outdir}/{BS}"

        # Step 1: Merge coverage data for temp0.fa
        print(f"Merging coverage data for {BS} temp0.fa")
        subprocess.run(f"cut -f1-3 {OUTF}/covs.txt > {OUTF}/combo.cov", shell=True)
        for cov_file in os.listdir(OUTF):
            if cov_file.startswith("CX-") and cov_file.endswith(".cov"):
                SN = cov_file.split('-')[1].replace(".cov", "")
                subprocess.run(f"cut -f4,5 {OUTF}/{cov_file} | sed 's/\tsample/\t{SN}/' | sed 's/^sample/{SN}/' | paste {OUTF}/combo.cov - > {OUTF}/temp.cov", shell=True)
                subprocess.run(f"mv {OUTF}/temp.cov {OUTF}/combo.cov", shell=True)
        subprocess.run(f"mv {OUTF}/combo.cov {OUTF}/combo1.cov", shell=True)

        # Step 2: Merge coverage data for temp1.fa
        print(f"Merging coverage data for {BS} temp1.fa")
        subprocess.run(f"cut -f1-3 {OUTF}/covs1.txt > {OUTF}/combo.cov", shell=True)
        for cov_file in os.listdir(OUTF):
            if cov_file.startswith("CXb-") and cov_file.endswith(".cov"):
                SN = cov_file.split('-')[1].replace(".cov", "")
                subprocess.run(f"cut -f4,5 {OUTF}/{cov_file} | sed 's/\tsample/\t{SN}/' | sed 's/^sample/{SN}/' | paste {OUTF}/combo.cov - > {OUTF}/temp.cov", shell=True)
                subprocess.run(f"mv {OUTF}/temp.cov {OUTF}/combo.cov", shell=True)
        subprocess.run(f"mv {OUTF}/combo.cov {OUTF}/combo1b.cov", shell=True)

    def run_metabat(self, BS):
        """ Run MetaBAT2 on merged coverage files. """
        OUTF = f"{self.outdir}/{BS}"
        
        # MetaBAT2 on temp0.fa with combo1.cov
        print(f"Running metabat2 on {BS} temp0.fa with combo1.cov")
        for k in range(15, 31):
            subprocess.run(f"metabat2 --seed {k} -i {OUTF}/temp0.fa -a {OUTF}/combo1.cov -t {self.threads} -o {OUTF}/bins/{k}CX1 -m {k*100}", shell=True)

        # MetaBAT2 on temp1.fa with combo1b.cov
        print(f"Running metabat2 on {BS} temp1.fa with combo1b.cov")
        for k in range(15, 31):
            subprocess.run(f"metabat2 --seed {k} -i {OUTF}/temp1.fa -a {OUTF}/combo1b.cov -t {self.threads} -o {OUTF}/bins/{k}CX1b -m {k*100}", shell=True)

    def run_checkm_binning(self, BS):
        """ Run CheckM2 and binning process. """
        OUTF = f"{self.outdir}/{BS}"
        checkm_output = f"{OUTF}/bins/temp"

        # Run CheckM2
        db_flag = f"--database_path {self.checkm_db}" if self.checkm_db else ""
        print(f"Running CheckM2 for {BS}")
        subprocess.run(f"checkm2 predict -x fa -i {OUTF}/bins/ -o {checkm_output} --force --remove_intermediates -t {self.threads} --tmpdir {self.tmpdir} {db_flag}", shell=True)

        # Filter bins based on CheckM2 quality report
        print(f"Filtering CheckM2 results for {BS}")
        subprocess.run(f"awk -F'\\t' '$2 >= 50 && $3 <= 5' {checkm_output}/quality_report.tsv > {OUTF}/worthwhile.tsv", shell=True)
        subprocess.run(f"mkdir {OUTF}/good; rm {OUTF}/good/*", shell=True)
        subprocess.run(f"cp -l $(cut -f1 {OUTF}/worthwhile.tsv | sed 's/$/.fa/' | sed 's:^:{OUTF}/bins/:') {OUTF}/good/", shell=True)

        # Run lingenome
        print(f"Running lingenome for {BS}")
        subprocess.run(f"lingenome {OUTF}/good {OUTF}/good.fasta FILENAME", shell=True)
        
        # Run akmer100b
        print(f"Running akmer100b for {BS}")
        subprocess.run(f"OMP_NUM_THREADS={self.threads} akmer100b {OUTF}/good.fasta {OUTF}/good.dm 13 ANI CHANCE GC LOCAL RC", shell=True)

        # Run spamw
        print(f"Running spamw for {BS}")
        subprocess.run(f"spamw2 {OUTF}/good.dm {OUTF}/L 0 {self.threads} D2", shell=True)

        # Generate worthwhile stats
        print(f"Generating worthwhile statistics for {BS}")
        cmd_process_bins = (
            f"awk -F'\\t' '{{printf \"%s\\t%s\\t%s\\t%s\\t%.3f\\n\", $1, $4, $2, $3, $12/10}}' {OUTF}/worthwhile.tsv > {OUTF}/worthwhile.stat"
        )
        subprocess.run(cmd_process_bins, shell=True)

    def run_bestmag(self, BS):
        """ Run bestmag and move the best bins to the final MAG directory. """
        OUTF = f"{self.outdir}/{BS}"
        bestmags_txt = f"{OUTF}/bestmags.txt"
        
        print(f"Running bestmag for {BS}")
        subprocess.run(f"bestmag2 {OUTF}/good.dm {OUTF}/L.txt {OUTF}/worthwhile.stat {OUTF}/bestmags.txt SELF S1", shell=True)
        # Move best bins to MAG directory
        print(f"Moving best bins for {BS}")
        for z in subprocess.run(f"tail -n+2 {bestmags_txt} | cut -f1", shell=True, capture_output=True, text=True).stdout.splitlines():
            subprocess.run(f"cp -l {OUTF}/good/{z}.fa {self.magdir}/{BS}_{z}.fa", shell=True)
            subprocess.run(f"echo '{self.magdir}/{BS}_{z}.fa' >> {self.magdir}/checks_{os.getenv('NODE')}.txt", shell=True)

    def run(self):
        """ Run all steps in sequence. """
        for BS in os.listdir(self.outdir):
            if os.path.isdir(f"{self.outdir}/{BS}") and BS!='mags':
                print(BS)
                self.run_sorenson_alignment(BS)
                self.merge_coverage_data(BS)
                self.run_metabat(BS)
                self.run_checkm_binning(BS)
                self.run_bestmag(BS)


# Example usage:
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

binning = CoAssemblyBinning(config)
binning.run()
