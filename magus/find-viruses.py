import os
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import glob
from pathlib import Path
import csv

class CheckVRunner:
    def __init__(self, asm_paths, checkv_db, combined_contig_file, filtered_contig_file, min_length, max_length, threads, quality, tmp_dir="tmp/run_checkv", config_file=None):
        self.asm_paths = asm_paths
        self.combined_contig_file = os.path.join(tmp_dir, combined_contig_file)
        self.filtered_contig_file = os.path.join(tmp_dir, filtered_contig_file)
        self.min_length = min_length
        self.max_length = max_length
        self.threads = threads
        self.tmp_dir = tmp_dir
        self.quality = set(quality)
        self.virus_dir = os.path.join("magus_viruses")
        self.checkv_db = checkv_db
        self.config_file = config_file
        self.contig_files = []
        os.makedirs(self.virus_dir, exist_ok=True)
        os.makedirs(tmp_dir, exist_ok=True)

    def find_contig_files_from_config(self):
        """Find contig files from a config file similar to call_orfs.py."""
        if not self.config_file:
            return []
        
        contig_files = []
        with open(self.config_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if not row or row[0].startswith('#') or len(row) < 2:
                    continue
                sample_id, contig_path = row[0], row[1]
                if os.path.exists(contig_path):
                    # Create unique identifier from sample_id and path for safety
                    path_part = Path(contig_path).parent.name
                    unique_id = f"{sample_id}____{path_part}".replace('/', '____').replace('\\', '____')
                    contig_files.append((unique_id, contig_path))
                else:
                    print(f"Warning: Contig file {contig_path} for sample {sample_id} not found.")
        return contig_files

    def find_contig_files_from_paths(self):
        """Find final.contigs.fa files from assembly paths, handling wildcards like dereplicate-genomes.py."""
        all_matches = []
        
        # Handle pipe-separated paths
        path_list = self.asm_paths.split('|') if '|' in self.asm_paths else [self.asm_paths]
        
        for path_pattern in path_list:
            path_pattern = path_pattern.strip()
            
            # Check if it's a single directory with all assemblies
            if os.path.isdir(path_pattern) and not any(wildcard in path_pattern for wildcard in ['*', '?', '[']):
                # Single directory - look for final.contigs.fa files recursively
                for contig_file in Path(path_pattern).rglob("final.contigs.fa"):
                    # Create unique identifier from path, using last 3 directory levels to keep it manageable
                    dir_path = contig_file.parent.resolve()
                    # Use last 3 directory levels for uniqueness while keeping it readable  
                    path_parts = dir_path.parts[-3:] if len(dir_path.parts) >= 3 else dir_path.parts
                    unique_id = '____'.join(path_parts).replace('/', '____').replace('\\', '____')
                    # Fallback if we get an empty unique_id
                    if not unique_id:
                        unique_id = contig_file.parent.name
                    all_matches.append((unique_id, str(contig_file.resolve())))
            else:
                # Wildcard pattern - expand and find contig files
                expanded_paths = glob.glob(path_pattern, recursive=True)
                for expanded_path in expanded_paths:
                    if os.path.isdir(expanded_path):
                        # Look for final.contigs.fa in this directory
                        contig_candidates = list(Path(expanded_path).glob("final.contigs.fa"))
                        if not contig_candidates:
                            # Try looking one level deeper
                            contig_candidates = list(Path(expanded_path).glob("*/final.contigs.fa"))
                        
                        for contig_file in contig_candidates:
                            # Create unique identifier from path, using last 3 directory levels to keep it manageable
                            full_path = str(contig_file.resolve())
                            dir_path = contig_file.parent.resolve()
                            # Use last 3 directory levels for uniqueness while keeping it readable
                            path_parts = dir_path.parts[-3:] if len(dir_path.parts) >= 3 else dir_path.parts
                            unique_id = '____'.join(path_parts).replace('/', '____').replace('\\', '____')
                            all_matches.append((unique_id, full_path))
        
        if not all_matches:
            print("No final.contigs.fa files found in the specified paths.")
        
        return all_matches

    def merge_contig_files(self):
        """Concatenate all contig files into a single merged file, with proper naming."""
        if self.config_file:
            self.contig_files = self.find_contig_files_from_config()
        else:
            self.contig_files = self.find_contig_files_from_paths()
        
        if not self.contig_files:
            print("No contig files found.")
            return
        
        print(f"Found {len(self.contig_files)} contig files to merge:")
        for unique_id, contig_path in self.contig_files:
            print(f"  {unique_id}: {contig_path}")
        
        # Create temporary files with renamed headers
        temp_files = []
        for unique_id, contig_path in self.contig_files:
            temp_file = os.path.join(self.tmp_dir, f"{unique_id}_renamed.fa")
            temp_files.append(temp_file)
            
            # Rename headers to include unique path-based identifier
            with open(contig_path, 'r') as infile, open(temp_file, 'w') as outfile:
                for line in infile:
                    if line.startswith('>'):
                        # Prepend unique identifier to contig header
                        header = line.strip()[1:]  # Remove '>'
                        outfile.write(f">{unique_id}___{header}\n")
                    else:
                        outfile.write(line)
        
        # Merge all temporary files
        cat_command = f"cat {' '.join(temp_files)} > {self.combined_contig_file}"
        subprocess.run(cat_command, shell=True, check=True)
        
        # Clean up temporary files
        for temp_file in temp_files:
            os.remove(temp_file)
        
        print(f"All contigs merged into {self.combined_contig_file}")

    def filter_contigs(self, filtered_file):
        """Filter merged contigs by length and save to a new file using `awk`."""
        awk_command = (
            f"awk '/^>/ {{if (seqlen >= {self.min_length} && seqlen <= {self.max_length}) "
            f"print header\"\\n\"seq; header=$0; seq=\"\"; seqlen=0; next}} "
            f"{{seq=seq$0; seqlen+=length($0)}} "
            f"END {{if (seqlen >= {self.min_length} && seqlen <= {self.max_length}) print header\"\\n\"seq}}' "
            f"{self.combined_contig_file} > {filtered_file}"
        )
        subprocess.run(awk_command, shell=True, check=True)
        print(f"Filtered contigs saved to {filtered_file}")

    def run_checkv_single(self, filtered_file):
        """Run CheckV on a filtered contig file."""
        output_dir = os.path.join(self.virus_dir, 'checkv_output')
        os.makedirs(output_dir, exist_ok=True)
        cmd = f"checkv end_to_end {filtered_file} {output_dir} -d {self.checkv_db} -t {self.threads}"
        subprocess.run(cmd, shell=True)
        print(f"CheckV analysis completed for {filtered_file}")
        return output_dir

    def process_results(self):
        quality_levels = {
            'C': "Complete",
            'H': "High-quality",
            'M': "Medium-quality",
            'L': "Low-quality"
        }
        
        selected_qualities = {quality_levels[q] for q in self.quality if q in quality_levels}
        print(f"Filtering CheckV output for quality levels: {selected_qualities}")
        
        binids = []
        for result_dir in os.listdir(self.virus_dir):
            result_path = os.path.join(self.virus_dir, result_dir)
            summary_file = os.path.join(result_path, "quality_summary.tsv")
            if os.path.exists(summary_file):
                good_file = os.path.join(result_path, "checkv_quality_summary_filtered.tsv")
                with open(good_file, "w") as good_out:
                    with open(summary_file) as summary_in:
                        header = summary_in.readline().strip()
                        good_out.write(f"{header}\tRepresentative\tUnique_ID\n")  # Add new columns
                        for line in summary_in:
                            columns = line.strip().split("\t")
                            contig_id, quality = columns[0], columns[7]
                            if quality in selected_qualities:
                                binids.append(contig_id)
                                # Extract unique identifier from contig ID if it has our naming convention
                                unique_id = contig_id.split('___')[0] if '___' in contig_id else 'unknown'
                                good_out.write(f"{line.strip()}\tNO\t{unique_id}\n")  # Initialize as 'NO'
        print(f"Processed CheckV results and saved filtered summary with selected qualities.")
        return binids

    def dereplicate_viruses(self, binids):
        """Subset viral contigs, run canolax5 to dereplicate, and update good.tsv with Representative sequences."""
        viruses_fna = os.path.join(self.virus_dir, "checkv_output/viruses.fna")
        tobe_dereplicated_viruses = os.path.join(self.virus_dir, "tobe_dereplicated_viruses.fasta")
        good_summary_output = os.path.join(self.virus_dir, "good_dereplicated_viruses.tsv")
        
        with open(viruses_fna, 'r') as infile, open(tobe_dereplicated_viruses, 'w') as outfile:
            write_flag = False
            for line in infile:
                if line.startswith('>'):
                    contig_id = line[1:].strip().split()[0]
                    write_flag = contig_id in binids
                if write_flag:
                    outfile.write(line)

        # Run canolax5 for dereplication
        canola_cmd = (
            f"canolax5 -db {tobe_dereplicated_viruses} -q {tobe_dereplicated_viruses} -o {self.virus_dir}/dereplicated_viruses.fasta -k 14 -t {self.threads} -fitC 0.02; rm {tobe_dereplicated_viruses}"
        )
        subprocess.run(canola_cmd, shell=True, check=True)

        # Collect dereplicated IDs
        dereplicated_ids = set()
        with open(os.path.join(self.virus_dir, "dereplicated_viruses.fasta"), 'r') as derep_in:
            for line in derep_in:
                if line.startswith('>'):
                    dereplicated_ids.add(line[1:].strip().split()[0])

        # Update good.tsv to mark representative sequences
        with open(os.path.join(self.virus_dir, "checkv_output/checkv_quality_summary_filtered.tsv"), 'r') as good_in, \
             open(good_summary_output, 'w') as good_out:
            for line in good_in:
                if line.startswith("contig_id"):  # Write header
                    good_out.write(line)
                else:
                    columns = line.strip().split("\t")
                    contig_id = columns[0]
                    representative = "YES" if contig_id in dereplicated_ids else "NO"
                    # Update the representative column
                    columns[-2] = representative  # Assuming Representative is second to last column
                    good_out.write("\t".join(columns) + "\n")
        print(f"Updated {good_summary_output} with representative information.")

    def run(self):
        self.merge_contig_files()
        self.filter_contigs(self.filtered_contig_file)
        self.run_checkv_single(self.filtered_contig_file)
        binids = self.process_results()
        self.dereplicate_viruses(binids)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CheckV on merged and filtered contigs from assembly directories")
    
    # Make asm_paths and config mutually exclusive
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--asm_paths", type=str, help="Pipe-separated list of assembly directory patterns with wildcards (e.g., 'path/to/*/asm1/|other/path/*/asm2') OR single directory with all assemblies")
    input_group.add_argument("--config", type=str, help="Tab-delimited config file with sample_id<TAB>contig_file_path")
    
    parser.add_argument("--combined_contig_file", type=str, default='all_contigs.fasta', help="Filename for merged contigs file")
    parser.add_argument("--filtered_contig_file", type=str, default='filtered_all_contigs.fasta', help="Filename for length-filtered merged contigs file")
    parser.add_argument("--min_length", type=int, default=500, help="Minimum length of contigs to include")
    parser.add_argument("--max_length", type=int, default=1000000000, help="Maximum length of contigs to include")
    parser.add_argument("--threads", type=int, default=28, help="Number of threads for CheckV")
    parser.add_argument("--quality", type=str, default="CHM", help="Viral contig levels to include (C [Complete], H [High], M [Medium], L [Low])")
    parser.add_argument("--tmp_dir", type=str, default="tmp/run_checkv", help="Temporary directory for storing intermediate files")
    parser.add_argument("--checkv_db", type=str, required=True, help="Path to the CheckV database directory")
    parser.add_argument(
        "--restart",
        choices=["cleanup"],
        help=(
            "Restart a step using existing file structures. Only 'cleanup' is supported, "
            "which resumes immediately after CheckV has run by processing results and dereplicating."
        ),
    )
    
    args = parser.parse_args()

    # Determine input method
    asm_paths = args.asm_paths if args.asm_paths else None
    config_file = args.config if args.config else None

    runner = CheckVRunner(
        asm_paths=asm_paths,
        checkv_db=args.checkv_db,
        combined_contig_file=args.combined_contig_file,
        filtered_contig_file=args.filtered_contig_file,
        min_length=args.min_length,
        max_length=args.max_length,
        threads=args.threads,
        quality=args.quality,
        tmp_dir=args.tmp_dir,
        config_file=config_file
    )
    # If restarting from cleanup, validate CheckV outputs exist and resume post-CheckV steps
    if args.restart == "cleanup":
        checkv_output_dir = os.path.join(runner.virus_dir, "checkv_output")
        quality_summary = os.path.join(checkv_output_dir, "quality_summary.tsv")
        viruses_fna = os.path.join(checkv_output_dir, "viruses.fna")

        if not os.path.exists(quality_summary) or not os.path.exists(viruses_fna):
            print(
                "Error: CheckV outputs not found. Expected files: "
                f"{quality_summary} and {viruses_fna}. Ensure CheckV has completed before restarting with --restart cleanup."
            )
            raise SystemExit(1)

        binids = runner.process_results()
        runner.dereplicate_viruses(binids)
    else:
        runner.run()
