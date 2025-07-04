#!/usr/bin/env python3

import os
import argparse
import shutil
import subprocess
from pathlib import Path
import glob

class Dereplicator:
	def __init__(self, mag_glob, tmp, threads, extension, wildcard, output, kmer_size):
		self.input_paths = [Path(p).resolve() for p in glob.glob(mag_glob, recursive=True)]
		self.tmp = Path(tmp)
		self.threads = threads
		self.extension = extension
		self.wildcard = wildcard
		self.derep_out = output
		self.kmer_size = kmer_size
		self.derep_tmp = self.tmp / "dereplicate_tmp"
		self.tmp_input_bins = self.derep_tmp / "input_bins"
		self.linearized_fa = self.derep_tmp / "lin.fa"
		self.output_reps = f"{self.derep_out}/reps.fa"
		self.output_log = f"{self.derep_out}/clus.tsv"

		self._init_tmp_dir()

	def _init_tmp_dir(self):
		if self.derep_tmp.exists():
			shutil.rmtree(self.derep_tmp)
		os.makedirs(self.tmp_input_bins, exist_ok=True)
		os.makedirs(self.derep_out, exist_ok=True)

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

def main():
	parser = argparse.ArgumentParser(description="Dereplicate MAGs using lingenome and canolax5.")
	parser.add_argument("-m", "--mag_dir", type=str, required=True, help="Path or glob to MAGs (e.g. asm/*/bins).")
	parser.add_argument("--tmp", type=str, default="tmp", help="Temporary working directory.")
	parser.add_argument("--threads", type=int, default=4, help="Number of threads for canolax5.")
	parser.add_argument("--extension", type=str, default="fa", help="File extension of MAGs (default: fa).")
	parser.add_argument("-w", "--wildcard", type=str, default="", help="Pattern to match anywhere in MAG path.")
	parser.add_argument("-o", "--output", type=str, default="dereplicated_genomes", help="Output directory (default: dereplicated_genomes).")
	parser.add_argument("-k", "--kmer_size", type=int, default=16, help="K-mer size for canolax5 (default: 16).")

	args = parser.parse_args()

	runner = Dereplicator(
		args.mag_dir,
		args.tmp,
		args.threads,
		args.extension,
		args.wildcard,
		args.output,
		args.kmer_size
	)
	runner.symlink_bins()
	runner.run_lingenome()
	runner.run_canolax5()

if __name__ == "__main__":
	main()

