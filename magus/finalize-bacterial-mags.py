import subprocess
import os
import argparse
import pandas as pd

class FinalMAGMerge:
    def __init__(self, tmp_dir,singleassembly_mag_dir, coasm_mag_dir, outdir, threads=28):
        self.original_mag_dir = singleassembly_mag_dir
        self.new_coasm_dir = coasm_mag_dir
        self.outdir = outdir
        self.merged_output = outdir + '/magus_consolidated_bac_arc_mags.fasta'
        self.threads = threads
        self.tmpdir = tmp_dir
        os.makedirs(self.tmpdir, exist_ok=True)
        os.makedirs(self.outdir, exist_ok=True)

    def load_config(self, config_path):
        config_df = pd.read_csv(config_path, sep='\t')
        return config_df.to_dict(orient='records')

    def clean_old_mags(self):
        print(f"Cleaning single assembly MAG directory: {self.original_mag_dir}")
        for file in os.listdir(self.original_mag_dir):
            if not file.endswith(".fa") and not file.endswith(".fasta"):
                os.remove(os.path.join(self.original_mag_dir, file))
        print(f"Cleaning co-assembled MAG directory: {self.new_coasm_dir}")
        for file in os.listdir(self.new_coasm_dir):
            if not file.endswith(".fa") and not file.endswith(".fasta"):
                os.remove(os.path.join(self.new_coasm_dir, file))

    def run_lingenome(self, input_dir, output_fasta, label):
        print(f"Running lingenome for {label} on {input_dir}")
        cmd = f"lingenome {input_dir} {output_fasta} FILENAME"
        subprocess.run(cmd, shell=True)

    def filter_canolax(self, original_fasta, new_fasta, filtered_output):
        print(f"Filtering co-assembled genomes using canolax5")
        cmd = (f"canolax4b -db {original_fasta} -q {new_fasta} -o {filtered_output} "
               f"-k 14 -t {self.threads} -local -filter .02")
        subprocess.run(cmd, shell=True)

    def merge_fasta_files(self, new_fasta, original_fasta):
        print(f"Merging {new_fasta} into {original_fasta}")
        with open(original_fasta, 'a') as orig_fasta, open(new_fasta, 'r') as new_f:
            orig_fasta.write(new_f.read())

    def annotate_bins(self, output_file):
        """Annotate the final bins with their source (single assembly or co-assembled)."""
        mags_in_merged = set()
        
        # Extract MAG names from the merged_output fasta file
        with open(self.merged_output, 'r') as merged_fasta:
            for line in merged_fasta:
                if line.startswith(">"):
                    mag_name = line[1:].strip().split()[0]
                    mags_in_merged.add(mag_name)

        # Determine source and annotate the file
        with open(output_file, 'w') as summary_out:
            summary_out.write("MAG\tSource\n")
            for mag in mags_in_merged:
                if f"{mag}.fa" in os.listdir(self.original_mag_dir):
                    source = "Single Assembly"
                elif f"{mag}.fa" in os.listdir(self.new_coasm_dir):
                    source = "Co-assembled"
                else:
                    source = "Unknown"
                summary_out.write(f"{mag}\t{source}\n")
        print(f"Annotation complete. Results written to {output_file}")

    def run(self):
        self.clean_old_mags()
        
        original_fasta = f"{self.tmpdir}/AllOrig.fasta"
        self.run_lingenome(self.original_mag_dir, original_fasta, label="original")
        
        new_fasta = f"{self.tmpdir}/NewCoAsm.fasta"
        self.run_lingenome(self.new_coasm_dir, new_fasta, label="co-assembly")
        
        filtered_fasta = f"{self.tmpdir}/CoAsm_to_keep.fasta"
        self.filter_canolax(original_fasta, new_fasta, filtered_fasta)
        
        self.merge_fasta_files(filtered_fasta, original_fasta)
        
        subprocess.run(f"mv {original_fasta} {self.merged_output}", shell=True)
        print(f"Final merged genome set saved as {self.merged_output}")
        
        self.annotate_bins(f"{self.outdir}/bacterial_mag_sources.tsv")

def main():
    parser = argparse.ArgumentParser(description="Merge and annotate final MAG files.")
    parser.add_argument('--singleassembly_mag_dir', type=str, default="asm/mags", help='Original MAGs directory (default: mags)')
    parser.add_argument('--coasm_mag_dir', type=str, default="coasm/mags", help='New co-assembled MAGs directory (default: coasm/mags)')
    parser.add_argument('--outdir', type=str, default="magus_output/magus_bacteria_archaea", help='Merged output file (default: magus_output/magus_bacteria_archaea/)')
    parser.add_argument('--threads', type=int, default=28, help='Number of threads (default: 28)')
    parser.add_argument('--tmp_dir', default = 'tmp/finalize-bacterial-mags',type=str, help='Temporary directory (default: tmp/finalize-bacterial-mags)')

    args = parser.parse_args()

    final_merge = FinalMAGMerge(
        singleassembly_mag_dir=args.singleassembly_mag_dir,
        coasm_mag_dir=args.coasm_mag_dir,
        outdir=args.outdir,
        tmp_dir=args.tmp_dir,
        threads=args.threads
    )
    final_merge.run()

if __name__ == '__main__':
    main()