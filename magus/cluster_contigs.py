import subprocess
import os

class ContigClustering:
	def __init__(self, config, contig_dir="contigs", combined_output="Contigs.fasta", threads=28):
		self.config = config
		self.contig_dir = contig_dir
		self.combined_output = combined_output
		self.threads = threads
		
		# Create the directory for storing the contigs
		os.makedirs(self.contig_dir, exist_ok=True)

	def collect_filtered_contigs(self):
		# Collect and rename filtered contigs (final.contigs.fa) for each sample
		for sample in self.config:
			sample_name = sample['filename']
			final_contig_path = f"asm/{sample_name}/final.contigs.fa"  # Path to the final contigs file
			output_contig_file = f"{self.contig_dir}/{sample_name}.contigs.fa"
			
			# Add the sample prefix to each contig
			cmd = f"sed 's/^>/>{sample_name}-/' {final_contig_path} > {output_contig_file}"
			print(f"Collecting final contigs for {sample_name}")
			subprocess.run(cmd, shell=True)

	def run_lingenome(self):
		# Run lingenome to combine all contigs into a single file
		cmd = f"lingenome {self.contig_dir}/ {self.combined_output} FILENAME"
		print(f"Running lingenome to generate combined contigs file: {self.combined_output}")
		subprocess.run(cmd, shell=True)

	def run_akmer100b(self):
		# Run akmer100b to generate the distance matrix
		distance_matrix = f"{self.combined_output.replace('.fasta', '')}.dm"
		cmd_akmer = f"OMP_NUM_THREADS={self.threads} akmer100b {self.combined_output} {distance_matrix} 16 ANI CHANCE GC LOCAL RC"
		print(f"Running akmer100b to generate distance matrix: {distance_matrix}")
		subprocess.run(cmd_akmer, shell=True)
		return distance_matrix

	def run_clustering(self):
		# Run spamw2 for clustering and bestmag for selecting best bins
		distance_matrix = self.run_akmer100b()  # Generate the distance matrix
		cluster_output = "clus/S"
		coasm_output = "coasm.reps"
		
		# Cluster the contigs
		cmd_spamw = f"spamw2 {distance_matrix} {cluster_output} 0 60 ALL NO2 WEIGHTED D2"
		print(f"Running spamw2 clustering on {distance_matrix}")
		subprocess.run(cmd_spamw, shell=True)
		
		# Select best bins with bestmag
		cmd_bestmag = f"bestmag {distance_matrix} {cluster_output}.txt NOSTAT REPS {coasm_output}"
		print(f"Running bestmag to select best bins")
		subprocess.run(cmd_bestmag, shell=True)
		
		# Create the coassembly task list
		cmd_generate_coasm = f"sed 's/\\.contigs\\t/\\t/g' {coasm_output} | sed 's/\\t$//' | grep -P '\\t' > coasm.todo"
		print(f"Generating coassembly task list in coasm.todo")
		subprocess.run(cmd_generate_coasm, shell=True)

	def run(self):
		self.collect_filtered_contigs()
		self.run_lingenome()
		self.run_clustering()

# Example usage for the samples:
config = [
    {
        'filename': 'SRR15373729', 
        'pe1': './qc/SRR15373729_1.fa.gz', 
        'pe2': './qc/SRR15373729_2.fa.gz'
    },
    {
        'filename': 'SRR15373733', 
        'pe1': './qc/SRR15373733_1.fa.gz', 
        'pe2': './qc/SRR15373733_2.fa.gz'
    }
]

clustering = ContigClustering(config)
clustering.run()
