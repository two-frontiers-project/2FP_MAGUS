import subprocess
import os
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

class QualityControl:
    def __init__(self, config, qc_dir="qc", max_workers=1, mode="local", slurm_config=None):
        self.config_path = config
        self.config = self.load_config(config)
        self.qc_dir = qc_dir
        self.max_workers = max_workers
        self.mode = mode
        self.slurm_config = self.load_slurm_config(slurm_config) if slurm_config else None
        self.qc_results = []
        os.makedirs(self.qc_dir, exist_ok=True)

    def load_config(self, config_path):
        config_df = pd.read_csv(config_path, sep='\t')
        if 'pe2' in config_df.columns:
            return config_df.to_dict(orient='records')
        else:
            return config_df.assign(pe2=None).to_dict(orient='records')
    
    def run_shi7_trimmer(self, sample):
        sample_name = sample['filename']
        r1 = sample['pe1']
        r2 = sample.get('pe2', None)
        
        if pd.isna(r2) or r2 is None:
            cmd = f"shi7_trimmer {r1} {self.qc_dir}/{sample_name} 75 12 FLOOR 4 ASS_QUALITY 20 CASTN 0 STRIP ADAP2 CTGTCTCTTATACA OUTFASTA"
            print(f"Running shi7_trimmer in single-end mode on {sample_name}")
            subprocess.run(cmd, shell=True)
            self.compress_output(sample_name, r2)
            self.qc_results.append({'filename': sample_name, 'pe1': f"{self.qc_dir}/{sample_name}.fa.gz", 'pe2': f"{self.qc_dir}/{sample_name}.R2.fa.gz" if r2 else None})

        else:
            cmd = f"shi7_trimmer {r1} {self.qc_dir}/{sample_name} 75 12 FLOOR 4 ASS_QUALITY 20 CASTN 0 STRIP ADAP2 CTGTCTCTTATACA R2 {r2} OUTFASTA"
            print(f"Running shi7_trimmer in paired-end mode on {sample_name}")
            subprocess.run(cmd, shell=True)
            self.compress_output(sample_name, r2)
            self.qc_results.append({'filename': sample_name, 'pe1': f"{self.qc_dir}/{sample_name}.R1.fa.gz", 'pe2': f"{self.qc_dir}/{sample_name}.R2.fa.gz" if r2 else None})

    def compress_output(self, sample_name, r2):
        for file in os.listdir(self.qc_dir):
            if file.startswith(sample_name):
                file_path = os.path.join(self.qc_dir, file)
                cmd = f"minigzip -4 {file_path}"
                print(f"Compressing {file_path}")
                subprocess.run(cmd, shell=True)

    def write_output_config(self):
        output_dir = os.path.dirname(self.config_path)
        os.makedirs(output_dir, exist_ok=True)
        df = pd.DataFrame(self.qc_results)
        df.to_csv(output_dir + '/post_qc_config', sep='\t', index=False)    

    def run(self):
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
