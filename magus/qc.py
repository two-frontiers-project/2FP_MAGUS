
import subprocess
import os
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import math

class QualityControl:
    def __init__(self, config, max_workers=1, output="configs/post_qc_config", mode="local", slurm_config=None):
        self.config = self.load_config(config)
        self.config_path = config
        self.qc_dir = "qc"
        self.output = output
        self.qc_results = []
        self.max_workers = max_workers
        self.mode = mode
        self.slurm_config = self.load_slurm_config(slurm_config) if slurm_config else None
        os.makedirs(self.qc_dir, exist_ok=True)

    def load_config(self, config_path):
        config_df = pd.read_csv(config_path, sep='\t')
        return config_df.to_dict(orient='records')
    
    def load_slurm_config(self, slurm_config_path):
        slurm_config_df = pd.read_csv(slurm_config_path, sep=',', header=None)
        slurm_config_dict = dict(zip(slurm_config_df[0], slurm_config_df[1]))
        return slurm_config_dict

    def split_config(self):
        # Split the original config into smaller batch config files based on max_workers
        config_df = pd.read_csv(self.config_path, sep='\t')
        batch_size = math.ceil(len(config_df) / self.max_workers)
        batch_files = []
        
        for i in range(self.max_workers):
            batch_df = config_df.iloc[i * batch_size: (i + 1) * batch_size]
            batch_file = os.path.join(os.path.dirname(self.config_path), f"batchconfig_qc_{i+1}.tsv")
            batch_df.to_csv(batch_file, sep='\t', index=False)
            batch_files.append(batch_file)
            print(f"Created batch configuration file: {batch_file}")

        return batch_files

    def write_deploy_script(self, batch_files):
        deploy_script_path = "deploy_qc.sh"
        
        with open(deploy_script_path, 'w') as f:
            f.write("#!/bin/bash\n\n")
            
            for batch_file in batch_files:
                # Adjust sbatch command to run QC on each batch
                sbatch_cmd = (
                    f"sbatch -p {self.slurm_config['queue']} "
                    f"-t {self.slurm_config['time']} "
                    f"--mem={self.slurm_config['memory']} "
                    f"-c {self.slurm_config['cpus_per_task']} "
                    f"--wrap='python qc.py --config {batch_file} --mode local'"
                )
                f.write(sbatch_cmd + "\n")
        
        print(f"SLURM deployment script written to {deploy_script_path}.")

    def run_local_mode(self):
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            executor.map(self.run_qc, self.config)

    def run_slurm_mode(self):
        batch_files = self.split_config()
        self.write_deploy_script(batch_files)
        print("Setup complete. Run 'deploy_qc.sh' to submit jobs to Slurm.")

    def run_qc(self, sample_config):
        # Implement the actual quality control steps here.
        sample_name = sample_config['sample_name']
        fastq_file = sample_config['fastq_path']
        
        # Example of running a fastp quality control command
        qc_cmd = f"fastp -i {fastq_file} -o {os.path.join(self.qc_dir, sample_name + '_qc.fastq')}"
        subprocess.run(qc_cmd, shell=True, check=True)
        print(f"Quality control complete for {sample_name}.")

    def run(self):
        if self.mode == "slurm":
            self.run_slurm_mode()
        else:
            self.run_local_mode()

def main():
    parser = argparse.ArgumentParser(description="Perform quality control on sequencing data.")
    parser.add_argument('--config', type=str, required=True, help='Path to the configuration TSV file')
    parser.add_argument('--slurm_config', type=str, help='Path to the Slurm configuration TSV file')
    parser.add_argument('--max_workers', type=int, default=1, help='Number of parallel workers (default: 1)')
    parser.add_argument('--mode', type=str, default="local", help="Execution mode: local or slurm (default: local)")

    args = parser.parse_args()

    qc = QualityControl(
        config=args.config,
        max_workers=args.max_workers,
        mode=args.mode,
        slurm_config=args.slurm_config
    )
    qc.run()

if __name__ == '__main__':
    main()

