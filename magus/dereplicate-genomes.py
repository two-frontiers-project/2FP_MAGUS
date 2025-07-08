#!/usr/bin/env python3

import os
import argparse
import shutil
import subprocess
from pathlib import Path
import glob

class Dereplicator:
	def __init__(self, mag_glob, tmp, threads, extension, wildcard, output, kmer_size, individual_reps):
		self.input_paths = [Path(p).resolve() for p in glob.glob(mag_glob, recursive=True)]
		self.tmp = Path(tmp)
		self.threads = threads
		self.extension = extension
		self.wildcard = wildcard
		self.derep_out = output
		self.kmer_size = kmer_size
		self.individual_reps = individual_reps
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

	def trim_large_contigs(self):
		"""Trim contigs larger than 2 gigabase pairs."""
		if not os.path.exists(self.linearized_fa):
			print(f"Warning: Linearized file {self.linearized_fa} not found. Skipping contig trimming.")
			return
		
		trimmed_fa = self.derep_tmp / "lin_trimmed.fa"
		max_bases = 2_000_000_000  # 2 gigabase pairs
		
		print(f"Trimming contigs larger than {max_bases:,} bases...")
		
		with open(self.linearized_fa, 'r') as infile, open(trimmed_fa, 'w') as outfile:
			current_header = None
			current_sequence = []
			
			for line in infile:
				line = line.strip()
				if line.startswith('>'):
					# Process previous contig if exists
					if current_header and current_sequence:
						sequence = ''.join(current_sequence)
						if len(sequence) > max_bases:
							print(f"Trimming {current_header} from {len(sequence):,} to {max_bases:,} bases")
							sequence = sequence[:max_bases]
						outfile.write(f"{current_header}\n")
						outfile.write(f"{sequence}\n")
					
					# Start new contig
					current_header = line
					current_sequence = []
				else:
					current_sequence.append(line)
			
			# Don't forget the last contig
			if current_header and current_sequence:
				sequence = ''.join(current_sequence)
				if len(sequence) > max_bases:
					print(f"Trimming {current_header} from {len(sequence):,} to {max_bases:,} bases")
					sequence = sequence[:max_bases]
				outfile.write(f"{current_header}\n")
				outfile.write(f"{sequence}\n")
		
		# Replace original file with trimmed version
		shutil.move(trimmed_fa, self.linearized_fa)
		print("Contig trimming completed.")

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
		"""Parse canolax5 output and symlink individual representative genome files."""
		if not self.individual_reps:
			return
			
		if not os.path.exists(self.output_log):
			print(f"Warning: Cluster log file {self.output_log} not found. Skipping individual representatives.")
			return
		
		# First, parse the linearized FASTA to get the mapping from genome names to file paths
		genome_to_file_map = {}
		if os.path.exists(self.linearized_fa):
			with open(self.linearized_fa, 'r') as f:
				for line in f:
					if line.startswith('>'):
						# lingenome format: >genome_name;FILENAME=actual_filename
						parts = line.strip().split(';')
						if len(parts) >= 2:
							genome_name = parts[0].lstrip('>')
							filename_part = parts[1]
							if filename_part.startswith('FILENAME='):
								actual_filename = filename_part.replace('FILENAME=', '')
								genome_to_file_map[genome_name] = actual_filename
		
		print(f"Found {len(genome_to_file_map)} genome name mappings")
		
		representatives = set()
		
		# Parse the canolax5 log file to identify representatives
		with open(self.output_log, 'r') as f:
			for line in f:
				line = line.strip()
				if line and not line.startswith('#'):
					# The format may vary, but typically the first column is the representative
					# and subsequent columns are cluster members
					parts = line.split('\t')
					if parts:
						rep_name = parts[0]
						representatives.add(rep_name)
		
		print(f"Found {len(representatives)} representatives: {list(representatives)}")
		
		# Symlink original files for each representative
		symlinked_count = 0
		for rep_name in representatives:
			if rep_name in genome_to_file_map:
				actual_filename = genome_to_file_map[rep_name]
				if actual_filename in self.original_file_map:
					original_file = self.original_file_map[actual_filename]
					dest_file = os.path.join(self.individual_reps_dir, actual_filename)
					try:
						# Remove existing file/symlink if it exists
						if os.path.exists(dest_file):
							os.remove(dest_file)
						# Create symlink to original file
						os.symlink(original_file, dest_file)
						symlinked_count += 1
						print(f"Symlinked representative: {rep_name} -> {actual_filename}")
					except Exception as e:
						print(f"Error symlinking {rep_name}: {e}")
				else:
					print(f"Warning: Original file not found for {actual_filename} (representative: {rep_name})")
			else:
				print(f"Warning: Genome name {rep_name} not found in genome mapping")
		
		print(f"Created {symlinked_count} individual representative genome symlinks in {self.individual_reps_dir}")

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

	args = parser.parse_args()

	runner = Dereplicator(
		args.mag_dir,
		args.tmp,
		args.threads,
		args.extension,
		args.wildcard,
		args.output,
		args.kmer_size,
		args.individual_reps
	)
	runner.symlink_bins()
	runner.run_lingenome()
	runner.trim_large_contigs() # Added this line
	runner.run_canolax5()
	runner.create_individual_representatives()

if __name__ == "__main__":
	main()

