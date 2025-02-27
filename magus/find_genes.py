import os
import subprocess
import pandas as pd
import argparse
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict

class CallAnnotateORFs:
	def __init__(self,database, assembly_dir = 'asm', config_dir = 'configs', output_dir = 'magus_output/annotation_output', tmp_dir = 'tmp', num_workers = 1, threads = 1):
		self.assembly_dir = assembly_dir
		self.config_dir = config_dir
		self.output_dir = output_dir
		self.database = database
		self.tmp_dir = tmp_dir
		self.num_workers = num_workers
		self.threads = threads
		
		self.annot_dir = os.path.join(self.tmp_dir, "scours/annot")
		self.manicure_dir = os.path.join(self.tmp_dir, "scours/manicure")
		self.hmmsearch_dir = os.path.join(self.tmp_dir, "scours")
		
		os.makedirs(self.annot_dir, exist_ok=True)
		os.makedirs(self.manicure_dir, exist_ok=True)
		os.makedirs(self.hmmsearch_dir, exist_ok=True)
		os.makedirs(self.output_dir, exist_ok=True)
		
		self.sample_map = self.load_sample_config()
	
	def load_sample_config(self):
		config_path = os.path.join(self.config_dir, "post_qc_config")
		config_df = pd.read_csv(config_path, sep='\t')
		return {row['filename']: os.path.join(self.assembly_dir, row['filename'], "final.contigs.fa") for row in config_df.to_dict(orient='records')}
	
	def call_orfs(self):
		def run_prodigal(sample, f):
			print('Running calling ORFs for %s'%f)
			fn = os.path.basename(f).replace('.fa', '')
			ffn_out = os.path.join(self.annot_dir, f"{sample}.ffn")
			faa_out = os.path.join(self.annot_dir, f"{sample}.faa")
			
			cmd = f"prodigal -p meta -q -i {f} -d {ffn_out} -a {faa_out} -o /dev/null"
			subprocess.run(cmd, shell=True, check=True)
			
			# Cleaning up the output
			with open(faa_out, 'r') as faa_file:
				faa_data = faa_file.read().replace('*', 'X').replace(' # ', '+').replace(' # ', '-')
			with open(os.path.join(self.manicure_dir, f"{sample}.faa"), 'w') as out_faa:
				out_faa.write(faa_data)
			
			with open(ffn_out, 'r') as ffn_file:
				ffn_data = ffn_file.read().replace(' # ', '+').replace(' # ', '-')
			with open(os.path.join(self.manicure_dir, f"{sample}.ffn"), 'w') as out_ffn:
				out_ffn.write(ffn_data)
			
			os.remove(ffn_out)
			os.remove(faa_out)
		
		with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
			executor.map(lambda item: run_prodigal(*item), self.sample_map.items())
	
	def filter_orfs(self):
		# Placeholder for pseudogene filtering implementation
		pass
	
	def annotate_orfs(self):
		def run_hmmsearch(sample, f):
			print('Running hmmsearch for %s'%f)
			tbl_out = os.path.join(self.hmmsearch_dir, f"{sample}.pfam")
			
			#cmd = f"hmmsearch -Z 16143 -o /dev/null --notextw --noali --nobias --cpu 1 --nseq_buffer 100000 --nhmm_buffer 1000 --tblout {tbl_out} /dev/shm/Pfam-A.hmm {f}"
			cmd = f"hmmsearch -Z 16143 -o /dev/null --notextw --noali --nobias --cpu 1 --tblout {tbl_out} /dev/shm/Pfam-A.hmm {f}"
			subprocess.run(cmd, shell=True, check=True)
		
		faa_files = {sample: os.path.join(self.manicure_dir, f"{sample}.faa") for sample in self.sample_map}
		
		with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
			executor.map(lambda item: run_hmmsearch(*item), faa_files.items())
	
	def process_hmmsearch_results(self):
		scores_dir = os.path.join(self.tmp_dir, "scours")
		os.makedirs(scores_dir, exist_ok=True)
		
		def process_pfam_file(sample, f):
			output_file = os.path.join(scores_dir, f"{sample}.pfam")
			
			with open(f, 'r') as pfam_file:
				lines = [line.strip().split() for line in pfam_file if not line.startswith('#')]
			
			sorted_lines = sorted(lines, key=lambda x: float(x[4]), reverse=True)  # Sort by e-value column
			
			with open(output_file, 'w') as out:
				for line in sorted_lines:
					out.write("\t".join([line[0], line[2], line[4], line[5], line[8]]) + "\n")
		
		pfam_files = {sample: os.path.join(self.hmmsearch_dir, f"{sample}.pfam") for sample in self.sample_map}
		
		with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
			executor.map(lambda item: process_pfam_file(*item), pfam_files.items())
	
	def generate_final_output(self):
		output_path = os.path.join(self.output_dir, "annotation_output.tsv")
		results = []
		
		for sample in self.sample_map:
			pfam_file = os.path.join(self.hmmsearch_dir, f"{sample}.pfam")
			if not os.path.exists(pfam_file) or os.stat(pfam_file).st_size == 0:
				continue
			
			gene_annotations = defaultdict(lambda: {"Sample": sample, "Contig_ID": None, "Gene_ID": None, "Annotations": [], "E-values": []})
			
			with open(pfam_file, 'r') as f:
				for line in f:
					if line.startswith('#'):
						continue
					parts = line.strip().split()
					gene_id = parts[0]
					annotation = parts[1]
					e_value = float(parts[4])
					
					contig_id = gene_id.split("+")[0]
					
					gene_annotations[gene_id]["Gene_ID"] = gene_id
					gene_annotations[gene_id]["Contig_ID"] = contig_id
					gene_annotations[gene_id]["Annotations"].append(annotation)
					gene_annotations[gene_id]["E-values"].append(e_value)
			
			for gene_id, data in gene_annotations.items():
				annotations_sorted = [x for _, x in sorted(zip(data["E-values"], data["Annotations"]))]
				evalues_sorted = sorted(data["E-values"])
				row = [
					data["Sample"],
					data["Contig_ID"],
					data["Gene_ID"],
					annotations_sorted[0] if len(annotations_sorted) > 0 else 'NA',
					evalues_sorted[0] if len(evalues_sorted) > 0 else 'NA',
					annotations_sorted[1] if len(annotations_sorted) > 1 else 'NA',
					evalues_sorted[1] if len(evalues_sorted) > 1 else 'NA',
					annotations_sorted[2] if len(annotations_sorted) > 2 else 'NA',
					evalues_sorted[2] if len(evalues_sorted) > 2 else 'NA'
				]
				results.append(row)
		
		df = pd.DataFrame(results, columns=[
			"Sample", "Contig_ID", "Gene_ID", "Top_Annotation", "E-value_1",
			"Second_Annotation", "E-value_2", "Third_Annotation", "E-value_3"
		])
		df.to_csv(output_path, sep="\t", index=False)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--assembly_dir", default="asm")
	parser.add_argument("--hmm_db", required=True)
	parser.add_argument("--config_dir", default="configs")
	parser.add_argument("--output_dir", default="magus/annotation_output")
	parser.add_argument("--tmp_dir", default="tmp")
	parser.add_argument("--num_workers", type=int, default=8)
  #  parser.add_argument("--threads", type=int, default=1)
	args = parser.parse_args()
	
	annotator = CallAnnotateORFs( args.hmm_db,args.assembly_dir,args.config_dir, args.output_dir, args.tmp_dir, args.num_workers)
	annotator.call_orfs()
	annotator.filter_orfs()
	annotator.annotate_orfs()
	annotator.process_hmmsearch_results()
	annotator.generate_final_output()

if __name__ == "__main__":
	main()
