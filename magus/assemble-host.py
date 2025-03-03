#!/usr/bin/env python3

import os
import argparse
import subprocess
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

class AssembleHost:
    def __init__(self, config_file, threads, max_workers, output_config,tmpdir,ksize):
        self.config_file = config_file
        self.output_config = output_config
        self.threads = threads
        self.max_workers = max_workers
        self.tmp_dir = tmpdir + '/hostfilt'
        self.ind_files_dir = os.path.join(self.tmp_dir, "indfiles")
        self.linearized_file = os.path.join(self.tmp_dir, "linearized.fa")
        self.meta_file = os.path.join(self.tmp_dir, "meta.dm")
        self.config_data = self.load_config(config_file)
        self.ksize = ksize

    def load_config(self, config_path):
        """Loads the config file and returns relevant fields as a dictionary."""
        config_df = pd.read_csv(config_path, sep='\t')
        if 'pe2' not in config_df.columns:
            config_df['pe2'] = None
        return config_df.to_dict(orient='records')

    def process_sample(self, entry):
        os.system(f'rm -rf {self.tmp_dir}; mkdir -p {self.tmp_dir}')
        """Processes a single sample using the filename and pe1 columns, handling compressed files if needed."""
        sample_id = os.path.basename(entry['filename'])
        sample_dir = os.path.join(self.ind_files_dir)
        os.makedirs(sample_dir, exist_ok=True)
        sample_fa = os.path.join(sample_dir, f"{sample_id}.fa")
        
        pe1_file = entry['pe1']
        if pe1_file[-3:] == '.gz':
            command = f"zcat {pe1_file} | head -n 2000000 > {sample_fa}"

        else:
            command = f"head -n 2000000 {pe1_file}> {sample_fa}"

        print(f"Processing sample: {sample_id}, Command: {command}")
        self.run_command(command)
    
    def run_lingenome(self):
        self.run_command(f"lingenome {self.ind_files_dir} {self.linearized_file} HEADFIX FILENAME; sleep 3")
    
    def run_akmer(self):
        self.run_command(f"OMP_NUM_THREADS={self.threads} akmer102 {self.linearized_file} {self.meta_file} 16 ANI CHANCE GC GLOBAL RC; sleep 3")
    
    def run_spamw(self):
        self.run_command(f"spamw2 {self.meta_file} {self.tmp_dir}/clustout {self.ksize} {self.threads} ALL WEIGHTED D3 ")

    def generate_host_assembly_config(self):
        clustout_file = f"{self.tmp_dir}/clustout.txt"
        """Generates a new config file based on cluster output."""
        if not os.path.exists(clustout_file):
            print(f"Cluster output file {clustout_file} not found.")
            return

        clust_df = pd.read_csv(clustout_file, sep='\t', usecols=["medoids"])
        clust_df.dropna(inplace=True)
        clust_df.rename(columns={"medoids": "sample_id"}, inplace=True)

        config_df = pd.DataFrame(self.config_data)
        new_config = config_df[config_df['filename'].isin(clust_df['sample_id'])]
        new_config.to_csv(self.output_config, sep='\t', index=False)
        print(f"Generated new config file: {self.output_config}")
    
    def run_command(self, command):
        """Executes a shell command."""
        print(f"Running: {command}")
        subprocess.run(command, shell=True, check=True)
    
    def run(self):
        """Main function to execute all steps with parallel processing."""
        for entry in self.config_data:
            self.process_sample(entry)
        self.run_lingenome()
        self.run_akmer()
        self.run_spamw()
        self.generate_host_assembly_config()

def main():
    parser = argparse.ArgumentParser(description="Assemble host/large genomes")
    parser.add_argument('--config', type=str, required=True, help='Path to the configuration TSV file')
    parser.add_argument('--max_workers', type=int, default=1, help='Number of parallel workers (default: 1)')
    parser.add_argument('--threads', type=int, default=14, help='Number of threads for processing (default: 14)')
    parser.add_argument('--ksize', type=int, default=0, help='Number of clusters to seek (default 0, so let MAGUS decide based on silhouette score).')
    parser.add_argument('--output_config', type=str, default="configs/host_assembly_config", help='Path to the output configuration file')
    parser.add_argument('--tmpdir', type=str, default="tmp", help='Directory where temp files should be written (default: ./tmp)')
    args = parser.parse_args()

    assembler = AssembleHost(args.config, args.threads, args.max_workers, args.output_config, args.tmpdir ,args.ksize)
    assembler.run()
if __name__ == "__main__":
    main()




