import subprocess
import os

class Assembly:
	def __init__(self, config, outdir="asm", magdir="mags", tmpdir="tmp", threads=28):
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
		print(f"Sample received: {sample}")
		sample_name = sample['filename']
		r1 = f"qc/{sample_name}.R1.fa.gz"  # Adjusted input path
		r2 = f"qc/{sample_name}.R2.fa.gz"  # Adjusted input path
		out_sample_dir = f"{self.outdir}/{sample_name}"
		
		# Run megahit assembly
		cmd = (f"megahit --k-min 75 --k-max 333 --k-step 6 --cleaning-rounds 1 "
			   f"--merge-level 100,.999 --min-count 1 --min-contig-len 1000 "
			   f"--continue -t {self.threads} -1 {r1} -2 {r2} -o {out_sample_dir}")
		
		print(f"Assembling {sample_name}")
		subprocess.run(cmd, shell=True)

	def post_process_contigs(self, sample_name):
		out_sample_dir = f"{self.outdir}/{sample_name}"
		contig_file = f"{out_sample_dir}/final.contigs.fa"
		filtered_contig_file = f"{out_sample_dir}/temp0.fa"
		
		# Filter contigs based on length and multiplicity
		cmd = (f"grep -A1 --no-group-separator 'multi=1\\.[5-9]\\|multi=[2-9]\\|multi=[0-9][0-9]' {contig_file} "
			   f"| awk '{{if (0==(NR % 2) && length >= 1000) {{print x; print $0}}; x=$0}}' > {filtered_contig_file}")
		
		subprocess.run(cmd, shell=True)
		print(f"Filtered contigs for {sample_name}")

	def run(self):
		for sample in self.config:
			#self.run_megahit(sample)
			self.post_process_contigs(sample['filename'])

# Example usage:
#config = [
#	{'filename': 'SRR30713783', 'pe1': 'qc/SRR30713783.R1.fa.gz', 'pe2': 'qc/SRR30713783.R2.fa.gz'}]
#config = [{'filename': 'testsample', 'pe1': 'qc/testsample.R1.fa.gz', 'pe2': 'qc/testsample.R2.fa.gz'},{'filename': 'testsample2', 'pe1': 'qc/testsample2.R1.fa.gz', 'pe2': 'qc/testsample2.R2.fa.gz'}]
config = [
    {
        'filename': 'SRR15373729', 
        'pe1': './qc/SRR15373729_1.fa.gz', 
        'pe2': './qc/SRR15373729_2.fa.gz'
    },
    {
        'filename': 'SRR15373733', 
        'pe1': './qc/SRR15373733_1.fa.gz', 
        'pe2': './qc/SRR15373733_2.fa.gz'
    }
]
assembly = Assembly(config)
assembly.run()
