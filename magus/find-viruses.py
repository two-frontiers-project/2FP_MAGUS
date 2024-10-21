import os
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import glob
import pandas as pd

class CheckVRunner:
    def __init__(self, asm_dir, coasm_dir,dblocs,combined_contig_file, filtered_contig_file, min_length, max_length, threads, quality, tmp_dir="tmp/run_checkv"):
        self.asm_dir = asm_dir
        self.coasm_dir = coasm_dir
        self.combined_contig_file = os.path.join(tmp_dir, combined_contig_file)
        self.filtered_contig_file = os.path.join(tmp_dir, filtered_contig_file)
        self.min_length = min_length
        self.max_length = max_length
        self.threads = threads
        self.tmp_dir = tmp_dir
        self.quality = set(quality)
        self.virus_dir = os.path.join("magus_viruses")
        self.checkv_db = self.get_db_location(dblocs, 'checkv')
        os.makedirs(self.virus_dir, exist_ok=True)
        os.makedirs(tmp_dir, exist_ok=True)

    def get_db_location(self, dblocs, db_name):
        """Retrieve the path of the specified database from the dblocs configuration file."""
        db_df = pd.read_csv(dblocs, header=None, index_col=0)
        return db_df.loc[db_name, 1]

    def merge_contig_files(self):
        """Concatenate all contig files from both coassembly and single assembly into a single merged file using `cat`."""
        coasm_files = glob.glob(os.path.join(self.coasm_dir, "*/final.contigs.fa"))
        single_asm_files = glob.glob(os.path.join(self.asm_dir, "*/final.contigs.fa"))
        
        if not coasm_files and not single_asm_files:
            print("No final.contigs.fa files found in either coassembly or single assembly directories.")
            return
        
        all_files = coasm_files + single_asm_files
        cat_command = f"cat {' '.join(all_files)} > {self.combined_contig_file}"
        subprocess.run(cat_command, shell=True, check=True)
        print(f"All contigs from coassembly and single assembly merged into {self.combined_contig_file}")

    def filter_contigs(self, filtered_file):
        """Filter merged contigs by length and save to a new file using `awk`."""
        awk_command = (
            f"awk '/^>/ {{if (seqlen >= {self.min_length} && seqlen <= {self.max_length}) "
            f"print header\"\\n\"seq; header=$0; seq=\"\"; seqlen=0; next}} "
            f"{{seq=seq$0; seqlen+=length($0)}} "
            f"END {{if (seqlen >= {self.min_length} && seqlen <= {self.max_length}) print header\"\\n\"seq}}' "
            f"{self.combined_contig_file} > {filtered_file}"
        )
        subprocess.run(awk_command, shell=True, check=True)
        print(f"Filtered contigs saved to {filtered_file}")

    def run_checkv_single(self, filtered_file):
        """Run CheckV on a filtered contig file."""
        output_dir = os.path.join(self.virus_dir, 'checkv_output')
        os.makedirs(output_dir, exist_ok=True)
        cmd = f"checkv end_to_end {filtered_file} {output_dir} -d {self.checkv_db} -t {self.threads}"
        subprocess.run(cmd, shell=True)
        print(f"CheckV analysis completed for {filtered_file}")
        return output_dir

    def process_results(self):
        quality_levels = {
            'C': "Complete",
            'H': "High-quality",
            'M': "Medium-quality",
            'L': "Low-quality"
        }
        
        selected_qualities = {quality_levels[q] for q in self.quality if q in quality_levels}
        print(f"Filtering CheckV output for quality levels: {selected_qualities}")
        
        binids = []
        for result_dir in os.listdir(self.virus_dir):
            result_path = os.path.join(self.virus_dir, result_dir)
            summary_file = os.path.join(result_path, "quality_summary.tsv")
            if os.path.exists(summary_file):
                good_file = os.path.join(result_path, "checkv_quality_summary_medium-high-complete.tsv")
                with open(good_file, "w") as good_out:
                    with open(summary_file) as summary_in:
                        header = summary_in.readline().strip()
                        good_out.write(f"{header}\tRepresentative\n")  # Add new column for Representative
                        for line in summary_in:
                            columns = line.strip().split("\t")
                            contig_id, quality = columns[0], columns[7]
                            if quality in selected_qualities:
                                binids.append(contig_id)
                                good_out.write(f"{line.strip()}\tNO\n")  # Initialize as 'NO'
        print(f"Processed CheckV results and saved filtered summary with selected qualities.")
        return binids

    def dereplicate_viruses(self, binids):
        """Subset viral contigs, run canolax5 to dereplicate, and update good.tsv with Representative sequences."""
        viruses_fna = os.path.join(self.virus_dir, "filtered_all_contigs/viruses.fna")
        tobe_dereplicated_viruses = os.path.join(self.virus_dir, "tobe_dereplicated_viruses.fasta")
        good_summary_output = os.path.join(self.virus_dir, "good_dereplicated_viruses.tsv")
        
        with open(viruses_fna, 'r') as infile, open(tobe_dereplicated_viruses, 'w') as outfile:
            write_flag = False
            for line in infile:
                if line.startswith('>'):
                    contig_id = line[1:].strip().split()[0]
                    write_flag = contig_id in binids
                if write_flag:
                    outfile.write(line)

        # Run canolax5 for dereplication
        canola_cmd = (
            f"canolax5 -db {tobe_dereplicated_viruses} -q {tobe_dereplicated_viruses} -o {self.virus_dir}/dereplicated_viruses.fasta -k 14 -t {self.threads} -fitC 0.02; rm {tobe_dereplicated_viruses}"
        )
        subprocess.run(canola_cmd, shell=True, check=True)

        # Collect dereplicated IDs
        dereplicated_ids = set()
        with open(os.path.join(self.virus_dir, "dereplicated_viruses.fasta"), 'r') as derep_in:
            for line in derep_in:
                if line.startswith('>'):
                    dereplicated_ids.add(line[1:].strip().split()[0])

        # Update good.tsv to mark representative sequences
        with open(os.path.join(self.virus_dir, "filtered_all_contigs/checkv_quality_summary_medium-high-complete.tsv"), 'r') as good_in, \
             open(good_summary_output, 'w') as good_out:
            for line in good_in:
                if line.startswith("contig_id"):  # Write header
                    good_out.write(line)
                else:
                    contig_id = line.split("\t")[0]
                    representative = "YES" if contig_id in dereplicated_ids else "NO"
                    good_out.write(f"{line.strip()}\t{representative}\n")
        print(f"Updated {good_summary_output} with representative information.")

    def run(self):
        self.merge_contig_files()
        self.filter_contigs(self.filtered_contig_file)
        self.run_checkv_single(self.filtered_contig_file)
        #binids = self.process_results()
        #self.dereplicate_viruses(binids)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CheckV on merged and filtered contigs")
    parser.add_argument("--asm_dir", type=str, default='asm', help="Directory containing singe assemblies")
    parser.add_argument("--coasm_dir", type=str, default='coasm', help="Directory containing coassemblies")
    parser.add_argument("--combined_contig_file", type=str, default='all_contigs.fasta', help="Filename for merged contigs file")
    parser.add_argument("--filtered_contig_file", type=str, default='filtered_all_contigs.fasta', help="Filename for length-filtered merged contigs file")
    parser.add_argument("--min_length", type=int, default=500, help="Minimum length of contigs to include")
    parser.add_argument("--max_length", type=int, default=1000000, help="Maximum length of contigs to include")
    parser.add_argument("--threads", type=int, default=28, help="Number of threads for CheckV")
    parser.add_argument("--quality", type=str, default="CHM", help="Viral contig levels to include (C [Complete], H [High], M [Medium], L [Low])")
    parser.add_argument("--tmp_dir", type=str, default="tmp/run_checkv", help="Temporary directory for storing intermediate files")
    parser.add_argument("--dblocs", type=str, required=True,default = 'configs/db_locs', help="Path to the dblocs configuration file")
    args = parser.parse_args()

    runner = CheckVRunner(
        asm_dir=args.asm_dir,
        coasm_dir=args.coasm_dir,
        combined_contig_file=args.combined_contig_file,
        filtered_contig_file=args.filtered_contig_file,
        min_length=args.min_length,
        dblocs=args.dblocs,
        max_length=args.max_length,
        threads=args.threads,
        quality=args.quality,
        tmp_dir=args.tmp_dir
    )
    runner.run()
