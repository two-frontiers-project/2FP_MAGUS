import subprocess
import os
import argparse
import pandas as pd
import math
from concurrent.futures import ThreadPoolExecutor

class Assembly:
	def __init__(self, config, outdir="asm", threads=14, max_workers=1, mode="local", slurm_config=None):
		self.config_path = config
		self.config = self.load_config(config)
		self.outdir = outdir
		self.threads = threads
		self.max_workers = max_workers
		self.mode = mode
		self.slurm_config = self.load_slurm_config(slurm_config) if slurm_config else None
		os.makedirs(self.outdir, exist_ok=True)

	def load_config(self, config_path):
		config_df = pd.read_csv(config_path, sep='\t')
		return config_df.to_dict(orient='records')

	def load_slurm_config(self, slurm_config_path):
		slurm_config_df = pd.read_csv(slurm_config_path, sep=',', header=None)
		slurm_config_dict = dict(zip(slurm_config_df[0], slurm_config_df[1]))
		return 	slurm_config_dict

	def split_config(self):
		# Split the original config into smaller batchconfig files based on max_workers
		config_df = pd.read_csv(self.config_path, sep='\t')
		batch_size = math.ceil(len(config_df) / self.max_workers)
		batch_files = []
		
		for i in range(self.max_workers):
			batch_df = config_df.iloc[i * batch_size: (i + 1) * batch_size]
			batch_file = os.path.join(os.path.dirname(self.config_path), f"batchconfig_{i+1}.tsv")
			batch_df.to_csv(batch_file, sep='\t', index=False)
			batch_files.append(batch_file)
			print(f"Created batch configuration file: {batch_file}")

		return batch_files

	def write_deploy_script(self, batch_files):
		deploy_script_path = "deploy_assembly.sh"
		
		with open(deploy_script_path, 'w') as f:
			f.write("#!/bin/bash\n\n")
			
			for batch_file in batch_files:
				# Adjusted sbatch command to call magus_assembly_slurm_helper.sh with arguments
				sbatch_cmd = (
					f"sbatch -p {self.slurm_config['queue']} "
					f"-t {self.slurm_config['time']} "
					f"--cpus-per-task={self.threads} "
					f"magus_assembly_slurm_helper.sh {batch_file} {self.threads}\n"
				)
				f.write(sbatch_cmd)

		os.chmod(deploy_script_path, 0o755)  # Make the deploy script executable
		print(f"Deploy script created: {deploy_script_path}")


	def run_megahit(self, sample):
		sample_name = sample['filename']
		r1 = sample['pe1']
		r2 = sample['pe2']
		out_sample_dir = f"{self.outdir}/{sample_name}"

		# Run megahit assembly
		cmd = (f"megahit --k-min 75 --k-max 333 --k-step 6 --cleaning-rounds 1 "
			   f"--merge-level 100,.999 --min-count 1 --min-contig-len 1000 "
			   f"--continue -t {self.threads} -1 {r1} -2 {r2} -o {out_sample_dir}")

		print(f"Assembling {sample_name}")
		subprocess.run(cmd, shell=True)

		# Post-process the contigs
		self.post_process_contigs(sample_name)

	def post_process_contigs(self, sample_name):
		out_sample_dir = f"{self.outdir}/{sample_name}"
		contig_file = f"{out_sample_dir}/final.contigs.fa"
		renamed_contig_file = f"{out_sample_dir}/renamed_final.contigs.fa"
		filtered_contig_file = f"{out_sample_dir}/temp0.fa"

		rename_cmd = f"sed 's/^>/>{sample_name}_singleassembly_/' {contig_file} > {renamed_contig_file}"
		subprocess.run(rename_cmd, shell=True, check=True)

		mv_cmd = f"mv {renamed_contig_file} {contig_file}"
		subprocess.run(mv_cmd, shell=True, check=True)

		filter_cmd = (
			f"grep -A1 --no-group-separator 'multi=1\\.[5-9]\\|multi=[2-9]\\|multi=[0-9][0-9]' {contig_file} "
			f"| awk '{{if (0==(NR % 2) && length >= 1000) {{print x; print $0}}; x=$0}}' > {filtered_contig_file}"
		)
		subprocess.run(filter_cmd, shell=True, check=True)
		print(f"Filtered contigs for {sample_name}")

	def run_local_mode(self):
		with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
			executor.map(self.run_megahit, self.config)

	def run_slurm_mode(self):
		batch_files = self.split_config()
		self.write_deploy_script(batch_files)
		print("Setup complete. Run 'deploy_assembly.sh' to submit jobs to Slurm.")

	def run(self):
		if self.mode == "slurm":
			self.run_slurm_mode()
		else:
			self.run_local_mode()

def main():
	parser = argparse.ArgumentParser(description="Run assembly with megahit on genomic data.")
	parser.add_argument('--config', type=str, required=True, help='Path to the configuration TSV file')
	parser.add_argument('--slurm_config', type=str, help='Path to the Slurm configuration TSV file')
	parser.add_argument('--max_workers', type=int, default=1, help='Number of parallel workers (default: 1)')
	parser.add_argument('--threads', type=int, default=14, help='Number of threads for megahit (default: 14)')
	parser.add_argument('--mode', type=str, default="local", help="Execution mode: local or slurm (default: local)")

	args = parser.parse_args()

	assembly = Assembly(
		config=args.config,
		slurm_config=args.slurm_config,
		max_workers=args.max_workers,
		threads=args.threads,
		mode=args.mode
	)
	assembly.run()

if __name__ == '__main__':
	main()
