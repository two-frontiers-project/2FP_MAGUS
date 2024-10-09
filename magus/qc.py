import subprocess
import os
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

class QualityControl:
    def __init__(self, config, max_workers=1, output="configs/post_qc_config"):
        self.config = self.load_config(config)
        self.qc_dir = "qc"
        self.output = output
        self.qc_results = []
        self.max_workers = max_workers
        os.makedirs(self.qc_dir, exist_ok=True)
    
    def load_config(self, config_path):
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
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            executor.map(self.run_shi7_trimmer, self.config)
        self.write_output_config()

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Run quality control with shi7_trimmer on genomic data.")
    parser.add_argument('--config', type=str, required=True, help='Path to the configuration TSV file')
    parser.add_argument('--max_workers', type=int, default=1, help='Max number of workers (default: 1)')
    parser.add_argument('--output', type=str, default="output/qc_config.tsv", help='Output file path for the new QC config file (default: output/qc_config.tsv)')

    # Parse arguments
    args = parser.parse_args()

    # Initialize and run the QualityControl instance
    qc = QualityControl(config=args.config, max_workers=args.max_workers, output=args.output)
    qc.run()

if __name__ == '__main__':
    main()
