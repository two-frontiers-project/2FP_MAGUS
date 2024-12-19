
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

    def run_slurm_mode(self):
        batch_files = self.split_config()
        self.write_deploy_script(batch_files)
        print("Setup complete. Run 'deploy_qc.sh' to submit jobs to Slurm.")


    def split_config(self):
        # Split the original config into smaller batch config files based on max_workers
        config_df = pd.read_csv(self.config_path, sep='\t')
        batch_size = math.ceil(len(config_df) / self.max_workers)
        batch_files = []
        
        for i in range(self.max_workers):
            batch_df = config_df.iloc[i * batch_size: (i + 1) * batch_size]
            batch_df.columns = ['filename', 'pe1', 'pe2']  # Set column names to filename, pe1, pe2
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
                    f"--mem={self.slurm_config['mem']} "
                    f"-c {self.slurm_config['cpus_per_task']} "
                    f"magus_qc_helper.sh {batch_file}"
                )
                f.write(sbatch_cmd + "\n")
        
        print(f"SLURM deployment script written to {deploy_script_path}.")

    def run_shi7_trimmer(self, sample):
        sample_name = sample['filename']
        r1 = sample['pe1']
        r2 = sample['pe2']
      
        cmd = f"shi7_trimmer {r1} {self.qc_dir}/{sample_name} 75 12 FLOOR 4 ASS_QUALITY 20 CASTN 0 STRIP ADAP2 CTGTCTCTTATACA R2 {r2} OUTFASTA"
        print(f"Running shi7_trimmer on {sample_name}")
        subprocess.run(cmd, shell=True)
        
        self.compress_output(sample_name)
        self.qc_results.append({'filename': sample_name, 'pe1': f"{self.qc_dir}/{sample_name}.R1.fa.gz", 'pe2': f"{self.qc_dir}/{sample_name}.R2.fa.gz"})

    def compress_output(self, sample_name):
        for file in os.listdir(self.qc_dir):
            if file.startswith(sample_name):
                file_path = os.path.join(self.qc_dir, file)
                cmd = f"minigzip -4 {file_path}"
                print(f"Compressing {file_path}")
                subprocess.run(cmd, shell=True)

    def write_output_config(self):
        output_dir = os.path.dirname(self.output)
        os.makedirs(output_dir, exist_ok=True)
        df = pd.DataFrame(self.qc_results)
        df.to_csv(self.output, sep='\t', index=False)


    def run(self):
        if self.mode == "slurm":
            self.run_slurm_mode()
        else:
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                executor.map(self.run_shi7_trimmer, self.config)
            self.write_output_config()

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

