#!/usr/bin/env python3

import os
import argparse
import shutil
import subprocess
from pathlib import Path
import glob

class Dereplicator:
	def __init__(self, mag_glob, tmp, threads, extension, wildcard, output, kmer_size, individual_reps, contig_level):
		self.input_paths = [Path(p).resolve() for p in glob.glob(mag_glob, recursive=True)]
		self.tmp = Path(tmp)
		self.threads = threads
		self.extension = extension
		self.wildcard = wildcard
		self.derep_out = output
		self.kmer_size = kmer_size
		self.individual_reps = individual_reps
		self.contig_level = contig_level
		self.derep_tmp = self.tmp / "dereplicate_tmp"
		self.tmp_input_bins = self.derep_tmp / "input_bins"
		self.linearized_fa = self.derep_tmp / "lin.fa"
		self.output_reps = f"{self.derep_out}/reps.fa"
		self.output_log = f"{self.derep_out}/clus.tsv"
		self.individual_reps_dir = f"{self.derep_out}/genome_representatives"

		# Store original file mapping for individual reps
		self.original_file_map = {}

		self._init_tmp_dir()

	def _init_tmp_dir(self):
		if self.derep_tmp.exists():
			shutil.rmtree(self.derep_tmp)
		os.makedirs(self.tmp_input_bins, exist_ok=True)
		os.makedirs(self.derep_out, exist_ok=True)
		if self.individual_reps:
			os.makedirs(self.individual_reps_dir, exist_ok=True)

	def symlink_bins(self):
		all_matches = []
		for base_path in self.input_paths:
			if base_path.is_dir():
				suffix = f".{self.extension}"
				for path in base_path.rglob(f"*{suffix}"):
					if self.wildcard in str(path):
						all_matches.append(path.resolve())

		if not all_matches:
			raise RuntimeError("No matching MAG files found.")

		common_prefix = os.path.commonpath([str(p) for p in all_matches])

		for path in all_matches:
			rel_path = path.relative_to(common_prefix)
			new_name = str(rel_path).replace("/", "-")
			target = self.tmp_input_bins / new_name
			os.symlink(path, target)
			
			# Store mapping for individual reps feature
			if self.individual_reps:
				self.original_file_map[new_name] = path

	def run_lingenome(self):
		cmd = f'lingenome {self.tmp_input_bins} {self.linearized_fa} HEADFIX FILENAME; sleep 5'
		#cmd = [
	#		"lingenome",
#			str(self.tmp_input_bins),
#			str(self.linearized_fa),
#			"HEADFIX",
#			"FILENAME; sleep 5",
#		]
		#subprocess.run(cmd, check=True)
		os.system(cmd)

	def run_linfasta_contig_level(self):
		"""Run linfasta on each genome individually, then concatenate results."""
		linfasta_dir = self.derep_tmp / "linfasta_output"
		os.makedirs(linfasta_dir, exist_ok=True)
		
		# Process each genome with linfasta
		for genome_file in self.tmp_input_bins.iterdir():
			if genome_file.is_symlink() and genome_file.suffix == f".{self.extension}":
				output_file = linfasta_dir / f"{genome_file.stem}_lin.fa"
				cmd = f'linfasta {genome_file} {output_file}'
				print(f"Running linfasta on {genome_file.name}")
				os.system(cmd)
		
		# Concatenate all linearized files
		with open(self.linearized_fa, 'w') as outfile:
			for lin_file in linfasta_dir.glob("*_lin.fa"):
				if lin_file.exists():
					with open(lin_file, 'r') as infile:
						outfile.write(infile.read())
		
		print(f"Concatenated all linearized contigs into {self.linearized_fa}")

	def run_canolax5(self):
		cmd = [
			"canolax5",
			"-k", str(self.kmer_size),
			"-db", str(self.linearized_fa),
			"-o", str(self.output_reps),
			"-local",
			"-fitC", "0.05",
			"-log", str(self.output_log)
		]
		env = os.environ.copy()
		env["OMP_NUM_THREADS"] = str(self.threads)
		subprocess.run(cmd, check=True, env=env)

	def create_individual_representatives(self):
		"""Parse canolax5 output and extract individual representative genome files from reps.fa."""
		if not self.individual_reps:
			return
			
		if not os.path.exists(self.output_reps):
			print(f"Warning: Representative file {self.output_reps} not found. Skipping individual representatives.")
			return
		
		# Parse the representative genomes and group by original genome
		genome_contigs = {}
		
		with open(self.output_reps, 'r') as f:
			current_header = None
			current_sequence = []
			
			for line in f:
				line = line.strip()
				if line.startswith('>'):
					# Save previous genome if exists
					if current_header and current_sequence:
						# Extract genome name from header (everything before first semicolon)
						genome_name = current_header.split(';')[0].lstrip('>')
						if genome_name not in genome_contigs:
							genome_contigs[genome_name] = []
						genome_contigs[genome_name].append((current_header, ''.join(current_sequence)))
					
					# Start new genome
					current_header = line
					current_sequence = []
				else:
					current_sequence.append(line)
			
			# Don't forget the last genome
			if current_header and current_sequence:
				genome_name = current_header.split(';')[0].lstrip('>')
				if genome_name not in genome_contigs:
					genome_contigs[genome_name] = []
				genome_contigs[genome_name].append((current_header, ''.join(current_sequence)))
		
		# Write individual genome files
		created_count = 0
		for genome_name, contigs in genome_contigs.items():
			output_file = os.path.join(self.individual_reps_dir, f"{genome_name}.fa")
			try:
				with open(output_file, 'w') as f:
					for header, sequence in contigs:
						f.write(f"{header}\n")
						f.write(f"{sequence}\n")
				created_count += 1
				print(f"Created representative genome: {genome_name} with {len(contigs)} contigs")
			except Exception as e:
				print(f"Error creating {genome_name}: {e}")
		
		print(f"Created {created_count} individual representative genome files in {self.individual_reps_dir}")

	def create_contig_level_representatives(self):
		"""For contig-level mode: separate representative contigs back into individual genome files."""
		if not self.individual_reps or not self.contig_level:
			return
			
		if not os.path.exists(self.output_reps):
			print(f"Warning: Representative file {self.output_reps} not found. Skipping contig-level representatives.")
			return
		
		# Parse the representative contigs and group by original genome
		genome_contigs = {}
		
		with open(self.output_reps, 'r') as f:
			current_header = None
			current_sequence = []
			
			for line in f:
				line = line.strip()
				if line.startswith('>'):
					# Save previous contig if exists
					if current_header and current_sequence:
						# Extract genome name from header (everything before first semicolon)
						genome_name = current_header.split(';')[0].lstrip('>')
						if genome_name not in genome_contigs:
							genome_contigs[genome_name] = []
						genome_contigs[genome_name].append((current_header, ''.join(current_sequence)))
					
					# Start new contig
					current_header = line
					current_sequence = []
				else:
					current_sequence.append(line)
			
			# Don't forget the last contig
			if current_header and current_sequence:
				genome_name = current_header.split(';')[0].lstrip('>')
				if genome_name not in genome_contigs:
					genome_contigs[genome_name] = []
				genome_contigs[genome_name].append((current_header, ''.join(current_sequence)))
		
		# Write individual genome files
		created_count = 0
		for genome_name, contigs in genome_contigs.items():
			output_file = os.path.join(self.individual_reps_dir, f"{genome_name}.fa")
			try:
				with open(output_file, 'w') as f:
					for header, sequence in contigs:
						f.write(f"{header}\n")
						f.write(f"{sequence}\n")
				created_count += 1
				print(f"Created representative genome: {genome_name} with {len(contigs)} contigs")
			except Exception as e:
				print(f"Error creating {genome_name}: {e}")
		
		print(f"Created {created_count} individual representative genome files in {self.individual_reps_dir}")

def main():
	parser = argparse.ArgumentParser(description="Dereplicate MAGs using lingenome and canolax5.")
	parser.add_argument("-m", "--mag_dir", type=str, required=True, help="Path or glob to MAGs (e.g. asm/*/bins).")
	parser.add_argument("--tmp", type=str, default="tmp", help="Temporary working directory.")
	parser.add_argument("--threads", type=int, default=4, help="Number of threads for canolax5.")
	parser.add_argument("--extension", type=str, default="fa", help="File extension of MAGs (default: fa).")
	parser.add_argument("-w", "--wildcard", type=str, default="", help="Pattern to match anywhere in MAG path.")
	parser.add_argument("-o", "--output", type=str, default="dereplicated_genomes", help="Output directory (default: dereplicated_genomes).")
	parser.add_argument("-k", "--kmer_size", type=int, default=16, help="K-mer size for canolax5 (default: 16).")
	parser.add_argument("--individual_reps", action="store_true", help="Create individual files for each representative genome with original contigs in genome_representatives/ folder.")
	parser.add_argument("--contig-level", action="store_true", help="Use contig-level deduplication with linfasta (for large eukaryotic contigs).")

	args = parser.parse_args()

	runner = Dereplicator(
		args.mag_dir,
		args.tmp,
		args.threads,
		args.extension,
		args.wildcard,
		args.output,
		args.kmer_size,
		args.individual_reps,
		args.contig_level
	)
	runner.symlink_bins()
	
	if args.contig_level:
		runner.run_linfasta_contig_level()
	else:
		runner.run_lingenome()
	
	runner.run_canolax5()
	
	if args.contig_level:
		runner.create_contig_level_representatives()
	else:
		runner.create_individual_representatives()

if __name__ == "__main__":
	main()

