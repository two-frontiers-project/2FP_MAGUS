import subprocess
import os

class FinalMAGMerge:
    def __init__(self, original_mag_dir="mags", new_coasm_dir="coasm/mags", merged_output="Complete_MAG_masterpiece.fasta", threads=28):
        self.original_mag_dir = original_mag_dir
        self.new_coasm_dir = new_coasm_dir
        self.merged_output = merged_output
        self.threads = threads
        self.checkm_db = None

    def clean_old_mags(self):
        """Clean up the original MAGs directory to remove non-FASTA files."""
        print(f"Cleaning single assembly MAG directory: {self.original_mag_dir}")
        for file in os.listdir(self.original_mag_dir):
            if not file.endswith(".fa") and not file.endswith(".fasta"):
                os.remove(os.path.join(self.original_mag_dir, file))
        print(f"Cleaning co-assembled MAG directory: {self.new_coasm_dir}")
        for file in os.listdir(self.new_coasm_dir):
            if not file.endswith(".fa") and not file.endswith(".fasta"):
                os.remove(os.path.join(self.new_coasm_dir, file))


    def run_lingenome(self, input_dir, output_fasta, label):
        """Run lingenome on the input directory and output to the specified FASTA."""
        print(f"Running lingenome for {label} on {input_dir}")
        cmd = f"lingenome {input_dir} {output_fasta} FILENAME"
        subprocess.run(cmd, shell=True)

    def filter_canolax(self, original_fasta, new_fasta, filtered_output):
        """Use canolax5 to filter out genomes that don't map to the old set."""
        print(f"Filtering co-assembled genomes using canolax5")
        cmd = (f"canolax4b -db {original_fasta} -q {new_fasta} -o {filtered_output} "
               f"-k 14 -t {self.threads} -local -filter .02")
        subprocess.run(cmd, shell=True)

    def merge_fasta_files(self, new_fasta, original_fasta):
        """Merge the filtered co-assembled genomes with the original set."""
        print(f"Merging {new_fasta} into {original_fasta}")
        with open(original_fasta, 'a') as orig_fasta, open(new_fasta, 'r') as new_f:
            orig_fasta.write(new_f.read())

    def annotate_bins(self, quality_report, output_file):
        """Annotate the final bins with source (co-assembly/single assembly)."""
        print(f"Annotating bins in {output_file} with source (co-assembly/single-assembly)")
        with open(quality_report, 'r') as quality_in, open(output_file, 'w') as summary_out:
            summary_out.write("Bin\tCompleteness\tContamination\tSource\n")
            for line in quality_in:
                if line.startswith("Bin"):
                    continue
                bin_name, completeness, contamination, *_ = line.strip().split('\t')
                source = "co-assembly" if "CX" in bin_name else "single-assembly"
                summary_out.write(f"{bin_name}\t{completeness}\t{contamination}\t{source}\n")

    def run(self):
        # Step 1: Clean the old MAGs directory (remove non-FASTA files)
        self.clean_old_mags()

        # Step 2: Run lingenome for single-sample bins
        original_fasta = "AllOrig.fasta"
        self.run_lingenome(self.original_mag_dir, original_fasta, label="original")

        # Step 3: Run lingenome for new co-assembled bins
        new_fasta = "NewCoAsm.fasta"
        self.run_lingenome(f"{self.new_coasm_dir}", new_fasta, label="co-assembly")

        # Step 4: Filter co-assembled genomes using canolax5
        filtered_fasta = "CoAsm_to_keep.fasta"
        self.filter_canolax(original_fasta, new_fasta, filtered_fasta)

        # Step 5: Merge filtered co-assembled genomes with original single-sample genomes
        self.merge_fasta_files(filtered_fasta, original_fasta)

        # Step 6: Move final merged file to the final output location
        subprocess.run(f"mv {original_fasta} {self.merged_output}", shell=True)
        print(f"Final merged genome set saved as {self.merged_output}")

        # Step 7: Annotate final selected bins with source (co-assembly or single-sample)
       # quality_report = f"{self.new_coasm_dir}/bins/temp/quality_report.tsv"
       # self.annotate_bins(quality_report, "quality_summary.tsv")
       # print("Annotation and quality summary complete.")

# Example usage:
final_merge = FinalMAGMerge()
final_merge.run()
