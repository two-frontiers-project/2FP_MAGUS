import argparse
import os
import subprocess

import pandas as pd


class ContigClustering:
    def __init__(
        self,
        config,
        asmdir,
        magdir,
        tmpdir,
        contig_dir="contigs",
        combined_output="Contigs.fasta",
        threads=28,
        ksize=0,
        assemblies_only=False,
    ):
        self.config = self.load_config(config)
        self.threads = threads
        self.ksize = ksize
        self.tmp_dir = tmpdir
        self.asmdir = asmdir
        self.magdir = magdir
        self.assemblies_only = assemblies_only
        self.contig_dir = self.tmp_dir + '/' + contig_dir
        self.combined_output = self.tmp_dir + '/' + combined_output
        # Create the directory for storing the contigs
        os.makedirs(self.tmp_dir, exist_ok=True)
        os.system('rm -rf %s' % self.contig_dir)
        os.makedirs(self.contig_dir, exist_ok=True)

    def load_config(self, config_path):
        config_df = pd.read_csv(config_path, sep='\t')
        return config_df.to_dict(orient='records')

    def collect_filtered_contigs(self):
        MAX_BASES = 1_000_000_000
        large_mag_contigs = {}  # Dictionary to store large MAG contigs
        if not self.assemblies_only:
            print('Getting genome bin sizes')
            # Step 1: Scan self.magdir for large MAGs and store contigs in dictionary
            for mag_file in os.listdir(self.magdir):
                mag_path = os.path.join(self.magdir, mag_file)

                # Use Bash to count bases in the MAG
                count_cmd = f"grep -v '^>' {mag_path} | tr -d '\\n' | wc -c"
                base_count = int(subprocess.check_output(count_cmd, shell=True).strip())

                if base_count > MAX_BASES:
                    # Extract the sample name (everything including & after the last underscore in the filename)
                    sample_name = mag_file.rsplit("_", 1)[-1].rsplit(".", 1)[0]  # Removes file extension too

                    # Extract contig IDs from the MAG
                    contig_ids_cmd = f"grep '^>' {mag_path} | sed 's/^>//'"
                    contig_ids = set(subprocess.check_output(contig_ids_cmd, shell=True).decode().splitlines())

                    # Store in dictionary
                    if sample_name not in large_mag_contigs:
                        large_mag_contigs[sample_name] = set()
                    large_mag_contigs[sample_name].update(contig_ids)
        else:
            print('Assemblies-only mode: skipping MAG size scan and filtering')

        # Step 2: Process each sample's final contigs
        for sample in self.config:
            sample_name = sample['filename']
            final_contig_path = os.path.abspath(f"{self.asmdir}/{sample_name}/final.contigs.fa")
            output_contig_file = os.path.abspath(f"{self.contig_dir}/{sample_name}.contigs.fa")

            if os.path.exists(final_contig_path) and os.path.getsize(final_contig_path) > 0:
                # Create symlink first
                cmd = f"ln -s {final_contig_path} {output_contig_file}"
                subprocess.run(cmd, shell=True)
                # Step 3: Filter out contigs from large MAGs before applying base limit filter
                if not self.assemblies_only and sample_name in large_mag_contigs:
                    print(f"Filtering out large MAG contigs from {sample_name}...")
                    filtered_output = f"{output_contig_file}.filtered"

                    # Write contig IDs to a temp file for grep exclusion
                    temp_contig_list = f"{output_contig_file}.contigs_to_remove.txt"
                    with open(temp_contig_list, "w") as f:
                        f.write("\n".join(large_mag_contigs[sample_name]) + "\n")

                    # Use Bash to remove the identified contigs from the sample
                    filter_cmd = f"awk 'BEGIN{{while((getline k < \"{temp_contig_list}\")>0) c[k]=1}} /^>/ {{print ($1 in c)?\"\":$0; next}} 1' {final_contig_path} > {filtered_output}"
                    subprocess.run(filter_cmd, shell=True)

                    os.replace(filtered_output, output_contig_file)  # Replace symlink with filtered file

                # Step 4: Check if the total base count exceeds the limit
                count_cmd = f"grep -v '^>' {output_contig_file} | tr -d '\\n' | wc -c"
                base_count = int(subprocess.check_output(count_cmd, shell=True).strip())

                if base_count > MAX_BASES:
                    print(f"File {output_contig_file} exceeds 1 gigabase ({base_count} bases). Trimming...")

                    temp_output = f"{output_contig_file}.trimmed"
                    trim_cmd = f"awk 'BEGIN {{bases=0}} !/^>/ {{bases += length($0); if (bases > {MAX_BASES}) exit}} {{print}}' {output_contig_file} > {temp_output}"
                    subprocess.run(trim_cmd, shell=True)

                    os.replace(temp_output, output_contig_file)  # Replace symlinked file with trimmed file

    def run_lingenome(self):
        # Run lingenome to combine all contigs into a single file
        cmd = f"lingenome {self.contig_dir}/ {self.combined_output} FILENAME"
        print(f"Running lingenome to generate combined contigs file: {self.combined_output}")
        subprocess.run(cmd, shell=True)

    def run_akmer102(self):
        # Run akmer102 to generate the distance matrix
        distance_matrix = f"{self.combined_output.replace('.fasta', '')}.dm"
        cmd_akmer = f"OMP_NUM_THREADS={self.threads} akmer102 {self.combined_output} {distance_matrix} 16 ANI CHANCE GC LOCAL RC"
        print(f"Running akmer102 to generate distance matrix: {distance_matrix}")
        subprocess.run(cmd_akmer, shell=True)
        return distance_matrix

    def run_clustering(self):
        # Run spamw2 for clustering and optionally bestmag for selecting representative clusters
        distance_matrix = self.run_akmer102()  # Generate the distance matrix
        cluster_output = f"{self.tmp_dir}/clus/S"
        coasm_output = f"{self.tmp_dir}/coasm.reps"

        os.makedirs(f"{self.tmp_dir}/clus/", exist_ok=True)

        # Cluster the contigs. For manual K, spamw2 is invoked directly for that K, which writes S.txt.
        # Silhouette/auto mode still uses spamw2's ALL iteration to explore multiple Ks and pick the
        # K that maximizes the silhouette score.
        if self.ksize > 0:
            cmd_spamw = f"spamw2 {distance_matrix} {cluster_output} {self.ksize} {self.threads} WEIGHTED NO2 D4"
        else:
            cmd_spamw = f"spamw2 {distance_matrix} {cluster_output} 0 {self.threads} ALL NO2 WEIGHTED D4"
        print(f"Running spamw2 clustering on {distance_matrix}")
        subprocess.run(cmd_spamw, shell=True)

        # Decide which spamw2 cluster file to pass along
        cluster_txt = f"{cluster_output}.txt"
        if self.ksize > 0:
            if not os.path.exists(cluster_txt):
                raise FileNotFoundError(
                    f"Requested K={self.ksize}, but expected spamw2 output {cluster_txt} was not found. "
                    "Check the spamw2 logs to confirm clustering finished successfully."
                )
            print(f"Using spamw2 cluster assignments for K={self.ksize}: {os.path.basename(cluster_txt)}")
        
        # Select best bins with bestmag using the default NOSTAT heuristics.
        bestmag_txt = f"{cluster_output}.bestmag.txt"
        cmd_bestmag = f"bestmag2 {distance_matrix} {cluster_txt} NOSTAT {bestmag_txt} REPS {coasm_output}"
        print(f"Running bestmag to select best bins before coassembly")
        subprocess.run(cmd_bestmag, shell=True)

        # Create the coassembly task list
        cmd_generate_coasm = f"sed 's/\\.contigs\\t/\\t/g' {coasm_output} | sed 's/\\t$//' | grep -P '\\t' > coasm.todo"
        print(f"Generating coassembly task list in coasm.todo")
        subprocess.run(cmd_generate_coasm, shell=True)

    def run(self):
        self.collect_filtered_contigs()
        self.run_lingenome()
        self.run_clustering()


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Run contig clustering pipeline on genomic data.")
    parser.add_argument('--config', type=str, required=True, help='Path to the configuration TSV file')
    parser.add_argument('--threads', type=int, default=28, help='Number of threads for tools (default: 28)')
    parser.add_argument('--contig_dir', type=str, default="contigs", help='Directory for storing contigs (default: contigs)')
    parser.add_argument('--combined_output', type=str, default="Contigs.fasta", help='Output file for combined contigs (default: Contigs.fasta)')
    parser.add_argument('--asmdir', type=str, default='asm/', help='Directory with assemblies. Default is ./asm')
    parser.add_argument('--magdir', type=str, default='asm/mags/', help='Directory with single assembled MAGs. Default is ./asm/mags')
    parser.add_argument('--tmpdir', type=str, default='tmp/cluster-contigs', help='Temp directory. Default tmp/cluster-contigs.')
    parser.add_argument('--ksize', type=int, default=0, help='Number of clusters to request from spamw2 (0 keeps silhouette mode).')
    parser.add_argument('--assemblies_only', action='store_true', help='Skip MAG filtering and cluster the final.contigs.fa files directly.')

    # Parse arguments
    args = parser.parse_args()

    # Initialize and run the ContigClustering instance
    clustering = ContigClustering(
        config=args.config,
        contig_dir=args.contig_dir,
        asmdir=args.asmdir,
        magdir=args.magdir,
        combined_output=args.combined_output,
        threads=args.threads,
        tmpdir=args.tmpdir,
        ksize=args.ksize,
        assemblies_only=args.assemblies_only,
    )
    clustering.run()


if __name__ == '__main__':
    main()
