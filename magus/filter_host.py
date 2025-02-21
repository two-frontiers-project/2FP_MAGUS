import os
import subprocess
import argparse
import pandas as pd
from pathlib import Path

class FilterHost:
    def __init__(self, config, xtree_db, xtree_outdir="xtree_host_out", outdir="configs", threads=14):
        self.config = self.load_config(config)
        self.xtree_db = xtree_db
        self.xtree_outdir = xtree_outdir
        self.outdir = outdir
        self.threads = threads
        
        os.makedirs(self.xtree_outdir, exist_ok=True)
        os.makedirs(self.outdir, exist_ok=True)
    
    def load_config(self, config_path):
        config_df = pd.read_csv(config_path, sep='\t')
        if 'pe2' in config_df.columns:
            return config_df.to_dict(orient='records')
        else:
            return config_df.assign(pe2=None).to_dict(orient='records')

    def run_xtree_alignment(self, file):
        """Run XTREE alignment on a given file."""
        xtree_cmd = (
            f"xtree ALIGN --seqs {file} --doforage --threads {self.threads} "
            f"--db {self.xtree_db} --ref-out {self.xtree_outdir}/{Path(file).stem}.ref "
            f"--cov-out {self.xtree_outdir}/{Path(file).stem}.cov "
            f"--perq-out {self.xtree_outdir}/{Path(file).stem}.perq --redistribute"
        )
        print(f"Running: {xtree_cmd}")
        subprocess.run(xtree_cmd, shell=True, check=True)

    def filter_reads(self, perq_file, input_file, output_file):
        """Filter reads from input file based on XTREE perq output."""
        to_remove = set()
        with open(perq_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4 and int(parts[3]) > 25:
                    to_remove.add(parts[0])
        
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            write = True
            for line in infile:
                if line.startswith('>'):
                    read_name = line.strip().split()[0][1:]
                    write = read_name not in to_remove
                if write:
                    outfile.write(line)

    def process_samples(self):
        """Process all samples by running XTREE and filtering reads."""
        new_config = []
        
        for row in self.config:
            sample = row['filename']
            pe1 = row['pe1']
            pe2 = row['pe2'] if pd.notna(row['pe2']) else None
            
            self.run_xtree_alignment(pe1)
            perq_file = f"{self.xtree_outdir}/{Path(pe1).stem}.perq"
            
            filtered_pe1 = f"{Path(pe1).parent}/{Path(pe1).stem}_filtered{Path(pe1).suffix}"
            self.filter_reads(perq_file, pe1, filtered_pe1)
            
            if pe2:
                self.run_xtree_alignment(pe2)
                perq_file2 = f"{self.xtree_outdir}/{Path(pe2).stem}.perq"
                filtered_pe2 = f"{Path(pe2).parent}/{Path(pe2).stem}_filtered{Path(pe2).suffix}"
                self.filter_reads(perq_file2, pe2, filtered_pe2)
            else:
                filtered_pe2 = ""
            
            new_config.append([sample, filtered_pe1, filtered_pe2])
        
        new_config_df = pd.DataFrame(new_config, columns=['filename', 'pe1', 'pe2'])
        new_config_path = os.path.join(self.outdir, "filtered_reads_config.tsv")
        new_config_df.to_csv(new_config_path, sep='\t', index=False)
        print(f"Filtered config written to {new_config_path}")
    
    def run(self):
        self.process_samples()


def main():
    parser = argparse.ArgumentParser(description="Filter host reads using XTREE alignment.")
    parser.add_argument('--config', type=str, default='configs/post_qc_config', help='Path to post QC config file')
    parser.add_argument('--xtree_db', type=str, required=True, help='Path to an XTREE database for alignment')
    parser.add_argument('--xtree_outdir', type=str, default='xtree_host_out', help='Directory for XTREE output')
    parser.add_argument('--outdir', type=str, default='configs', help='Output directory for filtered config')
    parser.add_argument('--threads', type=int, default=14, help='Number of threads for XTREE (default: 14)')
    
    args = parser.parse_args()
    
    filter_host = FilterHost(
        config=args.config,
        xtree_db=args.xtree_db,
        xtree_outdir=args.xtree_outdir,
        outdir=args.outdir,
        threads=args.threads
    )
    
    filter_host.run()

if __name__ == '__main__':
    main()
