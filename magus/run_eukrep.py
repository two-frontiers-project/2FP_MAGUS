import os
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import glob

class EukRepRunner:
    def __init__(self, coasm_dir, asm_dir, size_threshold, eukrep_output, max_workers=4):
        self.coasm_dir = coasm_dir
        self.asm_dir = asm_dir
        self.size_threshold = size_threshold
        self.eukrep_output = eukrep_output
        self.input_bins_dir = os.path.join(eukrep_output, "input_bins")
        self.max_workers = max_workers
        os.makedirs(self.input_bins_dir, exist_ok=True)

    def find_bins(self):
        """Locate all bins in the asm and coasm directories and symlink appropriate bins for EukRep."""
        # Collect all bins from asm and coasm directories
        bin_paths = glob.glob(os.path.join(self.coasm_dir, "*/bins/*fa")) + glob.glob(os.path.join(self.asm_dir, "*/bins/*fa"))
        good_paths = glob.glob(os.path.join(self.coasm_dir, "*/good/*fa")) + glob.glob(os.path.join(self.asm_dir, "*/good/*fa"))
        
        # Create a set of all good bins for easy reference
        good_bin_names = {os.path.basename(path) for path in good_paths}

        # Iterate over all bins and determine if they need to be symlinked
        for bin_path in bin_paths:
            # Determine if it's in asm or coasm and set up the output path accordingly
            if self.coasm_dir in bin_path:
                sample_base = os.path.relpath(bin_path, self.coasm_dir)
                sample_output_dir = os.path.join(self.input_bins_dir, "coasm", os.path.dirname(sample_base))
            else:
                sample_base = os.path.relpath(bin_path, self.asm_dir)
                sample_output_dir = os.path.join(self.input_bins_dir, "asm", os.path.dirname(sample_base))
            
            # Ensure the output directory exists
            os.makedirs(sample_output_dir, exist_ok=True)
            
            bin_name = os.path.basename(bin_path)
            # Symlink bins not in good or those exceeding the size limit
            if bin_name not in good_bin_names or self.is_bin_large(bin_path):
                symlink_source = os.path.abspath(bin_path)
                symlink_dest = os.path.join(sample_output_dir, bin_name)
                if not os.path.exists(symlink_dest):
                    os.symlink(symlink_source, symlink_dest)
                    print(f"Symlinked {symlink_source} to {symlink_dest}")


    def is_bin_large(self, bin_path):
        """Check if the bin is larger than the size threshold."""
        size_check_cmd = f"awk '{{if(/^>/) {{if(seqlen >= {self.size_threshold}) print seq}}; seqlen=0; seq=\"\"}} !/^>/ {{seqlen+=length($0); seq=seq$0}} END {{if(seqlen >= {self.size_threshold}) print seq}}' {bin_path}"
        size_check_result = subprocess.run(size_check_cmd, shell=True, capture_output=True, text=True)
        return bool(size_check_result.stdout.strip())

    def run_eukrep(self):
        """Run EukRep on all symlinked bin files in the input_bins directory."""
        futures = []
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            for bin_path in glob.glob(os.path.join(self.input_bins_dir, "*/*/bins/*")):
                # Extract folder names to create output file name
                parts = bin_path.split(os.sep)
                assembly_type = 'coassembly' if 'coasm' in parts else 'singleassembly'
                sample = parts[-3]
                bin_name = os.path.basename(bin_path).split(".")[0]
                
                # Define the output file name
                output_file = os.path.join(
                    self.eukrep_output,
                    f"{assembly_type}_{sample}_{bin_name}_eukrepcontigs.fa"
                )

                # Prepare and submit the EukRep command
                cmd = f"EukRep -i {bin_path} -o {output_file}"
                print(f"Running: {cmd}")
                futures.append(executor.submit(subprocess.run, cmd, shell=True))

            # Wait for all jobs to complete
            for future in as_completed(futures):
                future.result()
        print("EukRep processing completed for all bins.")


    def run(self):
        """Main function to prepare bins and run EukRep."""
        self.find_bins()
        self.run_eukrep()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run EukRep on bins from coassembly and single assembly")
    parser.add_argument("--coasm_dir", default='coasm',type=str, required=False, help="Directory containing coassembly bins")
    parser.add_argument("--asm_dir",default = "asm", type=str, required=False, help="Directory containing single assembly bins")
    parser.add_argument("--size_threshold", type=int, default=10000000, help="Size threshold for bins in base pairs (default: 10,000,000)")
    parser.add_argument("--eukrep_output",default='eukrep_output', type=str, required=False, help="Directory to save EukRep outputs")
    parser.add_argument("--max_workers", type=int, default=1, help="Maximum number of parallel jobs for EukRep")
    args = parser.parse_args()

    runner = EukRepRunner(
        coasm_dir=args.coasm_dir,
        asm_dir=args.asm_dir,
        size_threshold=args.size_threshold,
        eukrep_output=args.eukrep_output,
        max_workers=args.max_workers
    )
    runner.run()
