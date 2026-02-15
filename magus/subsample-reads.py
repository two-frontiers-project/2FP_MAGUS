import os
import subprocess
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

class ReadSubsampler:
    def __init__(self, config, outdir, out_config, depth, threads=4, max_workers=4):
        self.config = self.load_config(config)
        self.outdir = outdir
        self.out_config = out_config
        self.depth = depth
        self.threads = threads
        self.max_workers = max_workers
        self.seqtk_path = "seqtk"

        os.makedirs(self.outdir, exist_ok=True)
        os.makedirs(os.path.dirname(self.out_config), exist_ok=True)

    def load_config(self, config_file):
        """Load the sample configuration file into a DataFrame."""
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"Config file {config_file} not found.")

        df = pd.read_csv(config_file, sep="\t", dtype=str).fillna("")  # Ensure all values are strings, fill empty values with ""
        if df.shape[1] < 2 or df.shape[1] > 3:
            raise ValueError("Config file must have 2 or 3 columns: filename, pe1, [optional] pe2.")

        return df

    def detect_file_extension(self, filename):
        """Detect whether the file is FASTA (.fa) or FASTQ (.fq) based on its extension."""
        if filename.endswith(".fq.gz") or filename.endswith(".fastq.gz"):
            return "fq"
        elif filename.endswith(".fa.gz") or filename.endswith(".fasta.gz"):
            return "fa"
        elif filename.endswith(".fq") or filename.endswith(".fastq"):
            return "fq"
        elif filename.endswith(".fa") or filename.endswith(".fasta"):
            return "fa"
        else:
            raise ValueError(f"Unknown file format for {filename}. Must be .fq, .fa, .fastq, .fasta (compressed or uncompressed).")

    def subsample_reads(self, filename, pe1, pe2):
        """Subsample reads using seqtk. Handles single-end and paired-end cases."""
        file_extension = self.detect_file_extension(pe1)  # Detect whether the file is FASTA or FASTQ
        pe1_subsampled = os.path.join(self.outdir, f"{filename}_1.{file_extension}.gz")

        cmd1 = f"{self.seqtk_path} sample -s100 {pe1} {self.depth} | gzip > {pe1_subsampled}"

        try:
            subprocess.run(cmd1, shell=True, check=True)

            if pe2.strip():  # If pe2 is NOT empty, process it as paired-end
                pe2_subsampled = os.path.join(self.outdir, f"{filename}_2.{file_extension}.gz")
                cmd2 = f"{self.seqtk_path} sample -s100 {pe2} {self.depth} | gzip > {pe2_subsampled}"
                subprocess.run(cmd2, shell=True, check=True)
                return filename, pe1_subsampled, pe2_subsampled
            return filename, pe1_subsampled, ""  # If single-end, return empty pe2
        except subprocess.CalledProcessError as e:
            print(f"Error subsampling {filename}: {e}")
            return filename, None, None

    def run(self):
        """Run subsampling for all samples in parallel."""
        results = []
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(
                    self.subsample_reads,
                    row.iloc[0],  # filename
                    row.iloc[1],  # pe1
                    row.iloc[2] if len(row) > 2 else ""  # pe2 (if exists, otherwise empty)
                ): row.iloc[0]
                for _, row in self.config.iterrows()
            }

            for future in as_completed(futures):
                filename, pe1_sub, pe2_sub = future.result()
                if pe1_sub:
                    results.append([filename, pe1_sub, pe2_sub])

        # Save new config with headers
        new_config_df = pd.DataFrame(results, columns=["filename", "pe1", "pe2"])
        new_config_df.to_csv(self.out_config, sep="\t", index=False, header=True)

        print(f"Subsampling complete. New config saved to {self.out_config}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subsample sequencing reads using seqtk.")
    parser.add_argument("--config", required=True, help="Path to the input config file.")
    parser.add_argument("--outdir", default="subsampled_reads", help="Directory to store subsampled reads.")
    parser.add_argument("--out-config", "--out_config", dest="out_config", default="configs/subsampled_reads_config", help="Path to the output config file.")
    parser.add_argument("--depth", type=int, default=100000000, help="Number of reads to subsample to (default = 200000000).")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads.")
    parser.add_argument("--max-workers", "--max_workers", dest="max_workers", type=int, default=4, help="Maximum number of concurrent workers.")

    args = parser.parse_args()

    subsampler = ReadSubsampler(
        config=args.config,
        outdir=args.outdir,
        out_config=args.out_config,
        depth=args.depth,
        threads=args.threads,
        max_workers=args.max_workers
    )
    subsampler.run()
