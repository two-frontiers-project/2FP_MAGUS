import subprocess
import os
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

class QualityControl:
    def __init__(self, config, qc_dir="qc", max_workers=1, mode="local", slurm_config=None, seqtype="short"):
        self.config_path = config
        self.config = self.load_config(config)
        self.qc_dir = qc_dir
        self.max_workers = max_workers
        self.mode = mode
        self.slurm_config = self.load_slurm_config(slurm_config) if slurm_config else None
        self.seqtype = seqtype
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
        
        if self.seqtype == "long":
            cmd = f"shi7_trimmer {r1} {self.qc_dir}/{sample_name} 500 10 FLOOR 4 ASS_QUALITY 13"
            print(f"Running shi7_trimmer in long-read mode on {sample_name}")
            subprocess.run(cmd, shell=True)
            compressed = self.compress_long_output(sample_name)
            self.qc_results.append({'filename': sample_name, 'pe1': compressed, 'pe2': None})
            return

        if pd.isna(r2) or r2 is None:
            cmd = f"shi7_trimmer {r1} {self.qc_dir}/{sample_name} 75 12 FLOOR 4 ASS_QUALITY 20 CASTN 0 STRIP ADAP2 CTGTCTCTTATACA OUTFASTA"
            print(f"Running shi7_trimmer in single-end mode on {sample_name}")
            subprocess.run(cmd, shell=True)
            compressed = self.compress_short_output(sample_name, has_r2=False)
            self.qc_results.append({'filename': sample_name, 'pe1': compressed, 'pe2': None})

        else:
            cmd = f"shi7_trimmer {r1} {self.qc_dir}/{sample_name} 75 12 FLOOR 4 ASS_QUALITY 20 CASTN 0 STRIP ADAP2 CTGTCTCTTATACA R2 {r2} OUTFASTA"
            print(f"Running shi7_trimmer in paired-end mode on {sample_name}")
            subprocess.run(cmd, shell=True)
            r1_compressed, r2_compressed = self.compress_short_output(sample_name, has_r2=True)
            self.qc_results.append({'filename': sample_name, 'pe1': r1_compressed, 'pe2': r2_compressed})

    def compress_long_output(self, sample_name):
        fq_file = os.path.join(self.qc_dir, f"{sample_name}.fq")
        if os.path.exists(fq_file):
            cmd = f"minigzip -4 {fq_file}"
            print(f"Compressing {fq_file}")
            subprocess.run(cmd, shell=True)
        return f"{fq_file}.gz"

    def compress_short_output(self, sample_name, has_r2=False):
        if has_r2:
            r1_file = os.path.join(self.qc_dir, f"{sample_name}.R1.fa")
            r2_file = os.path.join(self.qc_dir, f"{sample_name}.R2.fa")

            if os.path.exists(r1_file):
                cmd = f"minigzip -4 {r1_file}"
                print(f"Compressing {r1_file}")
                subprocess.run(cmd, shell=True)

            if os.path.exists(r2_file):
                cmd = f"minigzip -4 {r2_file}"
                print(f"Compressing {r2_file}")
                subprocess.run(cmd, shell=True)

            return f"{r1_file}.gz", f"{r2_file}.gz"

        single_file = os.path.join(self.qc_dir, f"{sample_name}.fa")

        if os.path.exists(single_file):
            cmd = f"minigzip -4 {single_file}"
            print(f"Compressing {single_file}")
            subprocess.run(cmd, shell=True)

        return f"{single_file}.gz"

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
    parser.add_argument('--seqtype', type=str, choices=['long', 'short'], default='short',
                        help='Sequencing data type: long or short reads (default: short)')
    parser.add_argument('--outdir', dest='qc_dir', type=str, default='qc',
                        help='Directory where QC outputs are written (default: qc)')
    
    args = parser.parse_args()
    
    qc = QualityControl(
        config=args.config,
        max_workers=args.max_workers,
        mode=args.mode,
        slurm_config=args.slurm_config,
        seqtype=args.seqtype,
        qc_dir=args.qc_dir
    )
    qc.run()

if __name__ == '__main__':
    main()
