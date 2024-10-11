import os
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import glob

class CheckVRunner:
	def __init__(self, coasm_dir, combined_contig_file,filtered_contig_file, min_length, max_length, threads, parallel_jobs):
		self.coasm_dir = coasm_dir
		self.combined_contig_file = combined_contig_file
		self.filtered_contig_file = filtered_contig_file
		self.min_length = min_length
		self.max_length = max_length
		self.threads = threads
		self.parallel_jobs = parallel_jobs
		self.virus_dir = "checkv_output"
		self.processed_dir = "checkv_output"
		os.makedirs(self.virus_dir, exist_ok=True)
		os.makedirs(self.processed_dir, exist_ok=True)

	def merge_contig_files(self):
		"""Concatenate all contig files from both coassembly and single assembly into a single merged file using `cat`."""
		coasm_files = glob.glob(os.path.join(self.coasm_dir, "*/final.contigs.fa"))
		single_asm_files = glob.glob(os.path.join("asm", "*/final.contigs.fa"))
		
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
		output_dir = os.path.join(self.virus_dir, os.path.basename(filtered_file).replace(".fasta", ""))
		os.makedirs(output_dir, exist_ok=True)
		cmd = f"checkv end_to_end {filtered_file} {output_dir} -t {self.threads}"
		subprocess.run(cmd, shell=True)
		print(f"CheckV analysis completed for {filtered_file}")
		return output_dir

	def process_results(self):
		"""Process CheckV output and extract high-quality sequences."""
		binids = []
		for result_dir in os.listdir(self.virus_dir):
			result_path = os.path.join(self.virus_dir, result_dir)
			summary_file = os.path.join(result_path, "quality_summary.tsv")
			if os.path.exists(summary_file):
				good_file = os.path.join(result_path, "checkv_quality_summary_medium-high-complete.tsv")
				with open(good_file, "w") as good_out:
					with open(summary_file) as summary_in:
						for line in summary_in:
							if line.split("\t")[7] in ["Low-quality","Medium-quality", "High-quality", "Complete"]:
								good_out.write(line)
								binids.append(line[0])
				print(f"Processed CheckV results for {result_dir}")
		return binids

	def dereplicate_viruses(self, good_file):
		# subset viruses.fna ( in checkv directory with everything else) based on binids from previous function and write sequence 
		# run canola command to dereplicate viruses (replace filter with -fitC 0.02 from previous canola5x example)
		# save dereplicated viruses in root directory
		# subset good.tsv to contig ids from dereplicated viruses and put in root directory

	def run(self):
		self.merge_contig_files()
		self.filter_contigs(self.filtered_contig_file)
		self.run_checkv_single(self.filtered_contig_file)
		binids = self.process_results()
		self.dereplicate_viruses(binids)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Run CheckV on merged and filtered contigs")
	parser.add_argument("--coasm_dir", type=str, default = 'coasm', required=False, help="Directory containing coassembly contig files")
	parser.add_argument("--combined_contig_file", default = 'all_contigs.fasta',type=str, required=False, help="Output path for and name of the merged contigs file")
	parser.add_argument("--filtered_contig_file", default = 'filtered_all_contigs.fasta',type=str, required=False, help="Output path for and name of the length filtered, merged contigs file")
	parser.add_argument("--min_length", type=int, default = 500, required=False, help="Minimum length of contigs to include")
	parser.add_argument("--max_length", type=int, default = 1000000, required=False, help="Maximum length of contigs to include")
	parser.add_argument("--threads", type=int, default=8, help="Number of threads for CheckV")
	parser.add_argument("--parallel_jobs", type=int, default=1, help="Number of parallel jobs")
	args = parser.parse_args()

	runner = CheckVRunner(
		coasm_dir=args.coasm_dir,
		combined_contig_file=args.combined_contig_file,
		filtered_contig_file=args.filtered_contig_file,
		min_length=args.min_length,
		max_length=args.max_length,
		threads=args.threads,
		parallel_jobs=args.parallel_jobs
	)
	runner.run()
