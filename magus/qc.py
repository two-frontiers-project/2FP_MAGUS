import subprocess
import os
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

class QualityControl:
    def __init__(self, config, max_workers=1):
        self.config = self.load_config(config)
        self.qc_dir = "qc2"
        self.max_workers = max_workers
        os.makedirs(self.qc_dir, exist_ok=True)
    
    def load_config(self, config_path):
        # Load the TSV configuration file into a DataFrame and convert to a list of dictionaries
        config_df = pd.read_csv(config_path, sep='\t')
        return config_df.to_dict(orient='records')

    def run_shi7_trimmer(self, sample):
        sample_name = sample['filename']
        r1 = sample['pe1']
        r2 = sample['pe2']
        cmd = f"shi7_trimmer {r1} {self.qc_dir}/{sample_name} 75 12 FLOOR 4 ASS_QUALITY 20 CASTN 0 STRIP ADAP2 CTGTCTCTTATACA R2 {r2} OUTFASTA"
        print(f"Running shi7_trimmer on {sample_name}")
        subprocess.run(cmd, shell=True)
        self.compress_output(sample_name)

    def compress_output(self, sample_name):
        for file in os.listdir(self.qc_dir):
            if file.startswith(sample_name):
                file_path = os.path.join(self.qc_dir, file)
                cmd = f"minigzip -4 {file_path}"
                print(f"Compressing {file_path}")
                subprocess.run(cmd, shell=True)

    def run(self):
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            executor.map(self.run_shi7_trimmer, self.config)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Run quality control with shi7_trimmer on genomic data.")
    parser.add_argument('--config', type=str, required=True, help='Path to the configuration TSV file')
    parser.add_argument('--max_workers', type=int, default=1, help='Max number of workers (default: 1)')
    
    # Parse arguments
    args = parser.parse_args()

    # Initialize and run the QualityControl instance
    qc = QualityControl(config=args.config, max_workers=args.max_workers)
    qc.run()

if __name__ == '__main__':
    main()
