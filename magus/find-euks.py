import os
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import glob
import pandas as pd

class EukRepRunner:
    def __init__(self, coasm_dir,skip_eukrep,eukrepenv,skip_eukcc, asm_dir, size_threshold, euk_binning_outputdir, dblocs, max_workers=4, threads=8):
        self.coasm_dir = coasm_dir
        self.asm_dir = asm_dir
        self.skip_eukrep = skip_eukrep
        self.skip_eukcc = skip_eukcc
        self.size_threshold = size_threshold
        self.euk_binning_outputdir = euk_binning_outputdir
        self.eukcc_db = self.get_db_location(dblocs, 'eukccdb')
        self.max_workers = max_workers
        self.threads = threads
        self.eukrepenv = eukrepenv
        self.input_bins_dir = os.path.join(euk_binning_outputdir, "input_bins")
        os.makedirs(self.input_bins_dir, exist_ok=True)

    def get_db_location(self, dblocs, db_name):
        """Retrieve the path of the specified database from the dblocs configuration file."""
        db_df = pd.read_csv(dblocs, header=None, index_col=0)
        return db_df.loc[db_name, 1]

    def find_bins(self):
        """Locate all bins in the asm and coasm directories, symlink them, and store bin sizes."""
        self.bin_sizes = {}  # bin_name : bin_size

        bin_paths = glob.glob(os.path.join(self.coasm_dir, "*/bins/*fa")) + glob.glob(os.path.join(self.asm_dir, "*/bins/*fa"))
        good_paths = glob.glob(os.path.join(self.coasm_dir, "*/good/*fa")) + glob.glob(os.path.join(self.asm_dir, "*/good/*fa"))
        good_bin_names = {os.path.basename(path) for path in good_paths}

        for bin_path in bin_paths:
            if self.coasm_dir in bin_path:
                sample_base = os.path.relpath(bin_path, self.coasm_dir)
                sample_output_dir = os.path.join(self.input_bins_dir, "coasm", os.path.dirname(sample_base))
            else:
                sample_base = os.path.relpath(bin_path, self.asm_dir)
                sample_output_dir = os.path.join(self.input_bins_dir, "asm", os.path.dirname(sample_base))
            os.makedirs(sample_output_dir, exist_ok=True)

            bin_name = os.path.basename(bin_path)
            is_large, bin_size = self.is_bin_large(bin_path)
            self.bin_sizes[bin_name] = bin_size  # store for output

            if bin_name not in good_bin_names or is_large:
                symlink_source = os.path.abspath(bin_path)
                symlink_dest = os.path.join(sample_output_dir, bin_name)
                if not os.path.exists(symlink_dest):
                    os.symlink(symlink_source, symlink_dest)
                    print(f"Symlinked {symlink_source} to {symlink_dest}")


    def is_bin_large(self, bin_path):
        """Check bin size and return (is_large, size_in_bp)."""
        bin_size = 0
        with open(bin_path, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    bin_size += len(line.strip())
        is_large = bin_size >= self.size_threshold
        return is_large, bin_size

    def run_eukrep(self):
        """Run EukRep on all symlinked bin files in the input_bins directory."""
        futures = []
        print('Running EukRep...')
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            for bin_path in glob.glob(os.path.join(self.input_bins_dir, "*/*/bins/*")):
                parts = bin_path.split(os.sep)
                assembly_type = 'coassembly' if 'coasm' in parts else 'singleassembly'
                sample = parts[-3]
                bin_name = os.path.basename(bin_path)#.split(".")[0]

                os.makedirs(os.path.join(self.euk_binning_outputdir,f"{assembly_type}_{sample}_{bin_name}"), exist_ok=True)
                
                output_file = os.path.join(
                    self.euk_binning_outputdir,
                    f"{assembly_type}_{sample}_{bin_name}/EUKREP_{assembly_type}_{sample}_{bin_name}_eukrepcontigs.fa"
                )

                # Run EukRep within the specified conda environment
                cmd = f"bash -c 'source activate {self.eukrepenv} && EukRep -i {bin_path} -o {output_file}'"
                #print(f"Running EukRep: {cmd}")
                futures.append(executor.submit(subprocess.run, cmd, shell=True))

            for future in as_completed(futures):
                future.result()
        print("EukRep processing completed for all bins.")


    def run_eukcc(self):
        """Run EukCC on all EukRep output files."""
        futures = []
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            for bin_path in glob.glob(os.path.join(self.input_bins_dir, "*/*/bins/*")):
                parts = bin_path.split(os.sep)
                assembly_type = 'coassembly' if 'coasm' in parts else 'singleassembly'
                sample = parts[-3]
                bin_name = os.path.basename(bin_path)#.split(".")[0]

                os.makedirs(os.path.join(self.euk_binning_outputdir,f"{assembly_type}_{sample}_{bin_name}"), exist_ok=True)

                output_dir = os.path.join(
                    self.euk_binning_outputdir,
                    f"{assembly_type}_{sample}_{bin_name}/eukcc"
                )

                cmd = f"eukcc single --out {output_dir} --threads {self.threads} --db {self.eukcc_db} {bin_path}"
                print(f"Running EukCC: {cmd}")
                futures.append(executor.submit(subprocess.run, cmd, shell=True))

            for future in as_completed(futures):
                future.result()
        print("EukCC processing completed for all bins.")

    def process_euk_output(self):
        summary_data = []
        contigs_found = False  # Track if there are any eukrep results with contigs

        sample_dirs = glob.glob(os.path.join(self.euk_binning_outputdir, '*assembly*'))
        for sample_dir in sample_dirs:
            parts = sample_dir.split(os.sep)
            assembly_sample_bin = parts[1]

            assembly_type, sample_id_and_bin = assembly_sample_bin.split('_', 1)
            sample_id, bin_id = sample_id_and_bin.rsplit('_bin', 1)
            bin_id = f"bin{bin_id}"

            eukcc_file = os.path.join(sample_dir, 'eukcc/eukcc.csv')

            contig_file = os.path.join(
                self.euk_binning_outputdir,
                assembly_sample_bin,
                f"EUKREP_{assembly_type}_{sample_id}_{bin_id}_eukrepcontigs.fa"
            )

            # original bin input file
            input_bin_file = os.path.join(
                self.input_bins_dir,
                'coasm' if assembly_type == 'coassembly' else 'asm',
                sample_id,
                'bins',
                f"{bin_id}.fa"
            )

            completeness = contamination = float('nan')
            
            if os.path.exists(eukcc_file) and os.path.getsize(eukcc_file) > 0:
                eukcc_data = pd.read_csv(eukcc_file, sep='\t')
                if 'completeness' in eukcc_data.columns:
                    completeness = eukcc_data['completeness'].iloc[0]
                if 'contamination' in eukcc_data.columns:
                    contamination = eukcc_data['contamination'].iloc[0]

            # Count contigs in input bin
            bin_contig_count = 0
            if os.path.exists(input_bin_file):
                with open(input_bin_file, 'r') as f:
                    bin_contig_count = sum(1 for line in f if line.startswith('>'))

            # Count contigs in eukrep output
            eukrep_contig_count = 0
            if os.path.exists(contig_file):
                with open(contig_file, 'r') as f:
                    eukrep_contig_count = sum(1 for line in f if line.startswith('>'))

            eukrep_contig_perc = (eukrep_contig_count / bin_contig_count) * 100 if bin_contig_count > 0 else float('nan')

            row = {
                'assembly_type': assembly_type,
                'sample_id': sample_id,
                'bin_id': bin_id,
                'bin_size': self.bin_sizes.get(f"{bin_id}.fa", float('nan')),
                'completeness': completeness,
                'contamination': contamination,
                'bin_contig_count': bin_contig_count,
                'eukrep_contig_count': eukrep_contig_count,
                'eukrep_contig_perc': eukrep_contig_perc
            }

            if os.path.exists(contig_file):
                contigs_exist = os.path.getsize(contig_file) > 0
                row['eukrep_contigs_present'] = 'yes' if contigs_exist else 'no'
                contigs_found = True

            summary_data.append(row)

        if summary_data:
            summary_df = pd.DataFrame(summary_data)

            if not contigs_found:
                summary_df = summary_df.drop(columns=['eukrep_contigs_present'], errors='ignore')

            output_file = os.path.join(self.euk_binning_outputdir, 'eukaryotic_summary_table.csv')
            summary_df.to_csv(output_file, index=False)
            print(f"Summary table saved to {output_file}")
        else:
            print("No data found for summary table.")

    def run(self):
        self.find_bins()
        if not self.skip_eukrep:
            self.run_eukrep()
        if not self.skip_eukcc:
            self.run_eukcc()
        self.process_euk_output()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run EukRep and EukCC on bins from coassembly and single assembly")
    parser.add_argument("--coasm_dir", default='coasm', type=str, help="Directory containing coassembly bins")
    parser.add_argument("--asm_dir", default='asm', type=str, help="Directory containing single assembly bins")
    parser.add_argument("--size_threshold", type=int, default=10000000, help="Size threshold for bins in base pairs (default: 10,000,000)")
    parser.add_argument("--euk_binning_outputdir", default='magus_output/magus_euks', type=str, help="Directory to save EukRep and EukCC outputs")
    parser.add_argument("--dblocs", type=str, required=True, help="Path to the dblocs configuration file")
    parser.add_argument("--max_workers", type=int, default=1, help="Maximum number of parallel jobs for EukRep and EukCC")
    parser.add_argument("--threads", type=str, default=8, help="Number of threads per EukCC job")
    parser.add_argument("--skip_eukrep", type=str, default=True, help="True or False, skip eukrep step")
    parser.add_argument("--eukrep_env", type=str, default='/mnt/win/braden/.conda/envs/mamba/envs/eukrep-env', help="FULL PATH to the EukRep conda environment (eg /home/user/conda/envs/eukrep-env)")
    parser.add_argument("--skip_eukcc", type=str, default=False, help="True or False, skip eukcc step")
    args = parser.parse_args()

    skip_eukrep = args.skip_eukrep == 'true'
    skip_eukcc = args.skip_eukcc == 'true'

    runner = EukRepRunner(
        coasm_dir=args.coasm_dir,
        asm_dir=args.asm_dir,
        size_threshold=args.size_threshold,
        euk_binning_outputdir=args.euk_binning_outputdir,
        dblocs=args.dblocs,
        max_workers=args.max_workers,
        threads=args.threads,
        skip_eukrep=skip_eukrep,
        skip_eukcc=skip_eukcc,
        eukrepenv = args.eukrep_env
    )
    runner.run()
