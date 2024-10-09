import subprocess
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

class Binning:
	def __init__(self, config, outdir="asm", magdir="mags", tmpdir="tmp", threads=14, checkm_db=None, test_mode=False, max_workers=4):
		self.config = config
		self.outdir = outdir
		self.magdir = magdir
		self.tmpdir = tmpdir
		self.threads = threads
		self.checkm_db = checkm_db  # Custom CheckM database path
		self.test_mode = test_mode  # Flag for test mode
		self.max_workers = max_workers  # Maximum parallel workers
		
		# Create the output directories if they don't exist
		os.makedirs(self.outdir, exist_ok=True)
		os.makedirs(self.magdir, exist_ok=True)
		os.makedirs(self.tmpdir, exist_ok=True)

	def run_sorenson(self, sample_name):
		out_sample_dir = f"{self.outdir}/{sample_name}"
		temp0_file = f"{out_sample_dir}/temp0.fa"
		r1 = f"qc/{sample_name}.R1.fa.gz"
		r2 = f"qc/{sample_name}.R2.fa.gz"
		cov_file = f"{out_sample_dir}/covs.txt"
		
		# Run sorenson-g to generate coverage information
		cmd = f"sorenson-g -db {temp0_file} -qc -r1 {r1} -r2 {r2} -t {self.threads} -e 0.01 -o {cov_file}"
		print(f"Running sorenson-g for {sample_name}")
		subprocess.run(cmd, shell=True)

	def run_metabat(self, sample_name):
		out_sample_dir = f"{self.outdir}/{sample_name}"
		temp0_file = f"{out_sample_dir}/temp0.fa"
		cov_file = f"{out_sample_dir}/covs.txt"
		
		# Run MetaBAT2 with different seed values
		for k in [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 27, 30]:
			bin_output = f"{out_sample_dir}/bins/{k}"
			cmd = f"metabat2 --seed {k} -i {temp0_file} -a {cov_file} -t {self.threads} -o {bin_output} -m {k*100}"
			print(f"Running MetaBAT2 with seed {k} for {sample_name}")
			subprocess.run(cmd, shell=True)

	def run_checkm(self, sample_name):
		out_sample_dir = f"{self.outdir}/{sample_name}"
		bins_dir = f"{out_sample_dir}/bins"
		checkm_output = f"{bins_dir}/temp"
		
		# Use the provided database if available, otherwise default to the standard CheckM database
		if self.checkm_db:
			db_flag = f"--database_path {self.checkm_db}"
		else:
			db_flag = ""  # Default database path

		# Run CheckM to validate and filter the bins
		cmd = f"checkm2 predict -x fa -i {bins_dir} -o {checkm_output} --force --remove_intermediates -t {self.threads} --tmpdir {self.tmpdir} {db_flag}"
		print(f"Running CheckM for {sample_name} with database: {self.checkm_db if self.checkm_db else 'default'}")
		subprocess.run(cmd, shell=True)

	def filter_good_bins(self, sample_name):
		out_sample_dir = f"{self.outdir}/{sample_name}"
		checkm_report = f"{out_sample_dir}/bins/temp/quality_report.tsv"
		worthwhile_bins = f"{out_sample_dir}/worthwhile.tsv"
		
		# Test mode: Loosen completeness and contamination cutoffs
		if self.test_mode:
			completeness_cutoff = 0
			contamination_cutoff = 100
		else:
			completeness_cutoff = 50
			contamination_cutoff = 5

		# Filter the CheckM results based on cutoffs
		cmd = f"awk -F'\\t' '$2 >= {completeness_cutoff} && $3 <= {contamination_cutoff}' {checkm_report} > {worthwhile_bins}"
		subprocess.run(cmd, shell=True)
		print(f"Filtered good bins for {sample_name} (completeness >= {completeness_cutoff}%, contamination <= {contamination_cutoff}%)")

	def copy_good_bins(self, sample_name):
		out_sample_dir = f"{self.outdir}/{sample_name}"
		worthwhile_bins = f"{out_sample_dir}/worthwhile.tsv"
		bins_dir = f"{out_sample_dir}/bins"
		good_dir = f"{out_sample_dir}/good"
		
		# Create the 'good' directory if it doesn't exist
		os.makedirs(good_dir, exist_ok=True)

		# Extract the bin files to copy
		cmd = f"cut -f1 {worthwhile_bins} | sed 's/$/.fa/' | sed 's:^:{bins_dir}/:'"
		result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
		bin_files = result.stdout.strip().split('\n')
		
		# Check if bin_files is empty
		if not bin_files or bin_files == ['']:
			print(f"No bins to copy for {sample_name}.")
			return

		# Print the files to be copied
		print(f"Files to copy for {sample_name}: {bin_files}")
		
		# Copy the bins to the 'good' directory
		for bin_file in bin_files:
			cmd_copy = f"cp -l {bin_file} {good_dir}/"
			subprocess.run(cmd_copy, shell=True)
			print(f"Copied {bin_file} to {good_dir}")

	# Updated lingenome to output one directory up
	def run_lingenome(self, sample_name):
		out_sample_dir = f"{self.outdir}/{sample_name}/good"
		parent_dir = f"{self.outdir}/{sample_name}"  # Directory one level up
		good_fasta = f"{parent_dir}/good.fasta"  # Place good.fasta one level up
		cmd = f"lingenome {out_sample_dir} {good_fasta} FILENAME"
		print(f"Running lingenome for {sample_name}")
		subprocess.run(cmd, shell=True)

	# Update following paths to reflect new location of good.fasta
	def run_akmer100b(self, sample_name):
		parent_dir = f"{self.outdir}/{sample_name}"
		good_fasta = f"{parent_dir}/good.fasta"
		good_dm = f"{parent_dir}/good.dm"
		cmd = f"OMP_NUM_THREADS={self.threads} akmer100b {good_fasta} {good_dm} 13 ANI CHANCE GC LOCAL RC"
		print(f"Running akmer100b for {sample_name}")
		subprocess.run(cmd, shell=True)

	def run_spamw(self, sample_name):
		parent_dir = f"{self.outdir}/{sample_name}"
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
		parent_dir = f"{self.outdir}/{sample_name}"
		worthwhile_stat = f"{parent_dir}/worthwhile.stat"
		L_output = f"{parent_dir}/L.txt"
		bestmags_txt = f"{parent_dir}/bestmags.txt"
		good_dir = f"{parent_dir}/good"
		# Run bestmag
		cmd = f"bestmag2 {parent_dir}/good.dm {L_output} {worthwhile_stat} {bestmags_txt} SELF S1"
		print(f"Running bestmag for {sample_name}")
		print(cmd)
		subprocess.run(cmd, shell=True)

	def run_final_copy(self, sample_name):
		parent_dir = f"{self.outdir}/{sample_name}"
		mags_dir = self.magdir
		bestmags_txt = f"{parent_dir}/bestmags.txt"
		# Copy final good bins to MAG directory
		for z in subprocess.run(f"tail -n+2 {bestmags_txt} | cut -f1", shell=True, capture_output=True, text=True).stdout.splitlines():
			cmd = f"cp -l {parent_dir}/bins/{z}.fa {mags_dir}/{sample_name}_{z}.fa"
			subprocess.run(cmd, shell=True)
			print(f"Copied final bin {z} for {sample_name} to {mags_dir}")
		# Update the MAGDIR checks file with new bins
		checks_file = f"{mags_dir}/checks_{os.uname().nodename}.txt"
		cmd_append = f"tail -n+2 {bestmags_txt} | sed 's/^/{sample_name}_/' >> {checks_file}"
		subprocess.run(cmd_append, shell=True)
		print(f"Updated checks file: {checks_file}")
		# Cleanup: Remove intermediate files from bins and good directories
		cmd_cleanup = f"rm -r {parent_dir}/bins/* {parent_dir}/good/* {parent_dir}/temp*.fa {parent_dir}/good.dm"
		subprocess.run(cmd_cleanup, shell=True)
		print(f"Cleaned up intermediate files for {sample_name}")

	def run_sample(self, sample_name):
		"""Runs the entire binning pipeline for a single sample."""
		try:
			self.run_sorenson(sample_name)
			self.run_metabat(sample_name)
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

# Example usage:
config = [
	{'filename': 'sample_1', 'pe1': '/mnt/b/2FP_MAGUS/dev/qc/sample_1.R1.fa.gz', 'pe2': '/mnt/b/2FP_MAGUS/dev/qc/sample_1.R2.fa.gz'},
	{'filename': 'sample_2', 'pe1': '/mnt/b/2FP_MAGUS/dev/qc/sample_2.R1.fa.gz', 'pe2': '/mnt/b/2FP_MAGUS/dev/qc/sample_2.R2.fa.gz'},
	{'filename': 'sample_3', 'pe1': '/mnt/b/2FP_MAGUS/dev/qc/sample_3.R1.fa.gz', 'pe2': '/mnt/b/2FP_MAGUS/dev/qc/sample_3.R2.fa.gz'},
	{'filename': 'sample_4', 'pe1': '/mnt/b/2FP_MAGUS/dev/qc/sample_4.R1.fa.gz', 'pe2': '/mnt/b/2FP_MAGUS/dev/qc/sample_4.R2.fa.gz'},
	{'filename': 'sample_5', 'pe1': '/mnt/b/2FP_MAGUS/dev/qc/sample_5.R1.fa.gz', 'pe2': '/mnt/b/2FP_MAGUS/dev/qc/sample_5.R2.fa.gz'},
	{'filename': 'sample_6', 'pe1': '/mnt/b/2FP_MAGUS/dev/qc/sample_6.R1.fa.gz', 'pe2': '/mnt/b/2FP_MAGUS/dev/qc/sample_6.R2.fa.gz'},
	{'filename': 'sample_7', 'pe1': '/mnt/b/2FP_MAGUS/dev/qc/sample_7.R1.fa.gz', 'pe2': '/mnt/b/2FP_MAGUS/dev/qc/sample_7.R2.fa.gz'},
	{'filename': 'sample_8', 'pe1': '/mnt/b/2FP_MAGUS/dev/qc/sample_8.R1.fa.gz', 'pe2': '/mnt/b/2FP_MAGUS/dev/qc/sample_8.R2.fa.gz'},
	{'filename': 'sample_9', 'pe1': '/mnt/b/2FP_MAGUS/dev/qc/sample_9.R1.fa.gz', 'pe2': '/mnt/b/2FP_MAGUS/dev/qc/sample_9.R2.fa.gz'},
	{'filename': 'sample_10', 'pe1': '/mnt/b/2FP_MAGUS/dev/qc/sample_10.R1.fa.gz', 'pe2': '/mnt/b/2FP_MAGUS/dev/qc/sample_10.R2.fa.gz'}
]

binning = Binning(config, test_mode=False, max_workers=7) 
binning.run()
