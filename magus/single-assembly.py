import subprocess
import os
from concurrent.futures import ThreadPoolExecutor

class Assembly:
	def __init__(self, config, outdir="asm", magdir="mags", tmpdir="tmp", threads=14):
		self.config = config
		self.outdir = outdir
		self.magdir = magdir
		self.tmpdir = tmpdir
		self.threads = threads
		
		# Create the output directories if they don't exist
		os.makedirs(self.outdir, exist_ok=True)
		os.makedirs(self.magdir, exist_ok=True)
		os.makedirs(self.tmpdir, exist_ok=True)

	def run_megahit(self, sample):
		sample_name = sample['filename']
		r1 = f"qc/{sample_name}.R1.fa.gz"
		r2 = f"qc/{sample_name}.R2.fa.gz"
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
		filtered_contig_file = f"{out_sample_dir}/temp0.fa"
		
		# Filter contigs based on length and multiplicity
		cmd = (f"grep -A1 --no-group-separator 'multi=1\\.[5-9]\\|multi=[2-9]\\|multi=[0-9][0-9]' {contig_file} "
			   f"| awk '{{if (0==(NR % 2) && length >= 1000) {{print x; print $0}}; x=$0}}' > {filtered_contig_file}")
		
		subprocess.run(cmd, shell=True)
		print(f"Filtered contigs for {sample_name}")

	def run(self, max_workers=None):
		# Run assemblies in parallel using ThreadPoolExecutor or ProcessPoolExecutor
		with ThreadPoolExecutor(max_workers=max_workers) as executor:
			executor.map(self.run_megahit, self.config)

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

assembly = Assembly(config)
assembly.run(max_workers=7)  # Adjust max_workers based on desired parallelism
