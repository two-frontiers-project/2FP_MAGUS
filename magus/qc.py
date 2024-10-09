import subprocess
import os
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

class QualityControl:
    def __init__(self, config, max_workers=1, output="output/qc_config.tsv"):
        self.config = self.load_config(config)
        self.qc_dir = "qc2"
        self.output = output
        self.qc_results = []
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
        qc_r1 = os.path.join(self.qc_dir, f"{sample_name}_R1_qc.fq")
        qc_r2 = os.path.join(self.qc_dir, f"{sample_name}_R2_qc.fq")
        
        cmd = f"shi7_trimmer {r1} {qc_r1} 75 12 FLOOR 4 ASS_QUALITY 20 CASTN 0 STRIP ADAP2 CTGTCTCTTATACA R2 {r2} OUTFASTA {qc_r2}"
        print(f"Running shi7_trimmer on {sample_name}")
        subprocess.run(cmd, shell=True)
        
        self.compress_output(sample_name, qc_r1, qc_r2)
        self.qc_results.append({'filename': sample_name, 'qc_r1': f"{qc_r1}.gz", 'qc_r2': f"{qc_r2}.gz"})

    def compress_output(self, sample_name, qc_r1, qc_r2):
        for qc_file in [qc_r1, qc_r2]:
            cmd = f"minigzip -4 {qc_file}"
            print(f"Compressing {qc_file}")
            subprocess.run(cmd, shell=True)

    def write_output_config(self):
        # Ensure the output directory exists
        output_dir = os.path.dirname(self.output)
        os.makedirs(output_dir, exist_ok=True)
        
        # Write the QC results to the specified output file
        df = pd.DataFrame(self.qc_results)
        df.to_csv(self.output, sep='\t', index=False)
        print(f"Output configuration written to {self.output}")

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
