import subprocess
import os
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

class XTreeAligner:
    def __init__(self, dblocs, viral_file, euk_file, bacarc_file, outdir, threads=28, max_workers=1, taxonomy_db="databases/xtree_taxonomy_db"):
        self.viral_file = viral_file
        self.euk_file = euk_file
        self.bacarc_file = bacarc_file
        self.viral_db = self.get_db_location(dblocs, 'xtree_vir')
        self.euk_db = self.get_db_location(dblocs, 'xtree_euk')
        self.bacarc_db = self.get_db_location(dblocs, 'xtree_bac')
        self.outdir = outdir
        self.threads = threads
        self.max_workers = max_workers
        self.taxonomy_db = self.load_taxonomy_db(taxonomy_db)

        # Create necessary directories
        os.makedirs(self.outdir, exist_ok=True)

    def get_db_location(self, dblocs, db_name):
        """Retrieve the path of the specified database from the dblocs configuration file."""
        db_df = pd.read_csv(dblocs, header=None, index_col=0)
        return db_df.loc[db_name, 1]

    def load_taxonomy_db(self, taxonomy_db):
        taxonomy_mapping = {}
        with open(taxonomy_db, 'r') as f:
            for line in f:
                x = line.rstrip().split('\t')
                taxonomy_mapping[x[0]] = x[1]
        return taxonomy_mapping

    def run_xtree(self, db, seqs, output_folder):
        """Runs the xtree command with the given database and sequences"""
        ref_out = os.path.join(output_folder, 'xtree.ref')
        cov_out = os.path.join(output_folder, 'xtree.cov')
        perq_out = os.path.join(output_folder, 'xtree.perq')
        
        cmd = f"xtree --db {db} --seqs {seqs} --threads {self.threads} --ref-out {ref_out} --cov-out {cov_out} --perq-out {perq_out}"
        print(f"Running command: {cmd}")
        subprocess.run(cmd, shell=True, check=True)

    def run_classification(self, bintype):
        """Runs the xtree alignment for a given bintype (bacarc, viral, euk)"""
        if bintype == 'bacarc':
            print("Aligning bacterial and archaeal MAGs...")
            self.run_xtree(self.bacarc_db, self.bacarc_file, f"{self.outdir}/magus_bacteria_archaea")

        elif bintype == 'viral':
            print("Aligning viral MAGs...")
            self.run_xtree(self.viral_db, self.viral_file, f"{self.outdir}/magus_viruses")

        elif bintype == 'euk':
            print("Aligning eukaryotic MAGs...")
            self.run_xtree(self.euk_db, self.euk_file, f"{self.outdir}/magus_euks")

    def run(self, bintypes):
        """Run all alignments in parallel for the specified bintypes."""
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {executor.submit(self.run_classification, bintype): bintype for bintype in bintypes}
            for future in as_completed(futures):
                bintype = futures[future]
                try:
                    future.result()
                    print(f"Processing complete for {bintype}")
                except Exception as exc:
                    print(f"Error processing {bintype}: {exc}")

    def find_perq_files(self):
        """Find all .perq files in the output directory structure."""
        perq_files = []
        for root, dirs, files in os.walk(self.outdir):
            for file in files:
                if file.endswith('.perq'):
                    perq_files.append(os.path.join(root, file))
        return perq_files

    def parse_xtree_output(self, perq_files):
        species_dict = {}

        for perq_file in perq_files:
            with open(perq_file, 'r') as file:
                for line in file:
                    parts = line.split('\t')
                    contig = parts[0]
                    genome_id = ' '.join(parts[1])#.split('.')[0]  # Extract genome identifier
                    print(genome_id)
                    try:
                        species_name = self.taxonomy_db[genome_id]
                    except:
                        species_name = genome_id
                    try:
                        kmer_count = int(line.split()[-1])
                    except:
                        continue
                    if contig not in species_dict:
                        species_dict[contig] = []
                    species_dict[contig].append((species_name, kmer_count))

        # For each contig, sort by kmer count and keep top 5 species
        summary = []
        for contig, species_hits in species_dict.items():
            top_species = sorted(species_hits, key=lambda x: x[1], reverse=True)[:5]
            top_species_data = [f"{species} ({kmers} kmers)" for species, kmers in top_species]
            while len(top_species_data) < 5:
                top_species_data.append('None')  # Fill in with 'None' if fewer than 5 species
            summary.append([contig] + top_species_data)

        # Create dataframe
        df = pd.DataFrame(summary, columns=['sample_id', 'topspecies1', 'topspecies2', 'topspecies3', 'topspecies4', 'topspecies5'])
        return df


def main():
    parser = argparse.ArgumentParser(description="Run xtree alignment for viral, bacterial/archaeal, and eukaryotic MAGs.")
    parser.add_argument('--viral_file', type=str, default="magus_output/magus_viruses/dereplicated_viruses.fasta", help='Path to the viral MAG file')
    parser.add_argument('--euk_file', type=str, default="magus_output/magus_euks/potential_euks.fa", help='Path to the eukaryotic MAG file')
    parser.add_argument('--bacarc_file', type=str, default="magus_output/magus_bacteria_archaea/magus_consolidated_bac_arc_mags.fasta", help='Path to the bacterial/archaeal MAG file')
    parser.add_argument('--outdir', type=str, default="magus_output", help='Output directory for the results')
    parser.add_argument('--threads', type=int, default=28, help='Number of threads for xtree')
    parser.add_argument('--max-workers', '--max_workers', dest='max_workers', type=int, default=1, help='Maximum number of parallel workers')
    parser.add_argument('--bintypes', nargs='+', choices=['bacarc', 'viral', 'euk'], default=['bacarc', 'viral', 'euk'], help='Bintypes to align')
    parser.add_argument("--dblocs", type=str, default='configs/db_locs', help="Path to the dblocs configuration file")
    parser.add_argument("--taxonomy_db", type=str, default='databases/xtree_taxonomy_db', help="Path to the taxonomy database")

    args = parser.parse_args()

    aligner = XTreeAligner(
        viral_file=args.viral_file,
        euk_file=args.euk_file,
        bacarc_file=args.bacarc_file,
        outdir=args.outdir,
        threads=args.threads,
        max_workers=args.max_workers,
        dblocs=args.dblocs,
        taxonomy_db=args.taxonomy_db
    )

    # Run alignments
    #aligner.run(bintypes=args.bintypes)

    # Find and parse perq files
    perq_files = aligner.find_perq_files()
    if perq_files:
        df = aligner.parse_xtree_output(perq_files)
        df.to_csv(os.path.join(args.outdir, 'mag_classification_summary.csv'), index=False)
        print(df)
    else:
        print("No perq files found for parsing.")


if __name__ == '__main__':
    main()
