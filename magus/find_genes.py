import os
import subprocess
import pandas as pd
import argparse
import logging
from concurrent.futures import ThreadPoolExecutor

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class CallAnnotateORFs:
    def __init__(self, assembly_dir, config_dir, output_dir, tmp_dir, num_workers, threads, hmm_file=None):
        self.assembly_dir = assembly_dir
        self.config_dir = config_dir
        self.output_dir = output_dir
        self.tmp_dir = tmp_dir
        self.num_workers = num_workers
        self.threads = threads
        self.hmm_file = hmm_file
        
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
        try:
            config_df = pd.read_csv(config_path, sep='\t')
            return {row['sample_id']: os.path.join(self.assembly_dir, row['sample_id'], "final.contigs.fa") for row in config_df.to_dict(orient='records')}
        except Exception as e:
            logging.error(f"Error loading sample config from {config_path}: {str(e)}")
            raise
    
    def call_orfs(self):
        logging.info("Calling ORFs with Prodigal...")
        def run_prodigal(sample, f):
            try:
                fn = os.path.basename(f).replace('.fa', '')
                ffn_out = os.path.join(self.annot_dir, f"{sample}.ffn")
                faa_out = os.path.join(self.annot_dir, f"{sample}.faa")
                
                cmd = f"prodigal -p meta -q -i {f} -d {ffn_out} -a {faa_out} -o /dev/null -c {self.threads}"
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
                logging.info(f"Processed ORFs for sample {sample}")
            except Exception as e:
                logging.error(f"Error processing ORFs for sample {sample}: {str(e)}")
                raise
        
        with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
            executor.map(lambda item: run_prodigal(*item), self.sample_map.items())
        logging.info("ORF calling completed")
    
    def filter_orfs(self):
        # Placeholder for pseudogene filtering implementation
        logging.info("Skipping ORF filtering (not implemented)")
        pass
    
    def annotate_orfs(self):
        if not self.hmm_file or not os.path.exists(self.hmm_file):
            logging.info(f"HMM file not provided or not found at {self.hmm_file}. Skipping annotation.")
            return
        logging.info("Annotating ORFs with HMMsearch...")
        def run_hmmsearch(sample, f):
            try:
                tbl_out = os.path.join(self.hmmsearch_dir, f"{sample}.pfam")
                cmd = f"./hmmsearch-g -Z 16143 -o /dev/null --notextw --noali --nobias --cpu 1 --nseq_buffer 100000 --nhmm_buffer 1000 --tblout {tbl_out} {self.hmm_file} {f}"
                subprocess.run(cmd, shell=True, check=True)
                logging.info(f"Annotated ORFs for sample {sample}")
            except Exception as e:
                logging.error(f"Error annotating ORFs for sample {sample}: {str(e)}")
                raise
        
        faa_files = {sample: os.path.join(self.manicure_dir, f"{sample}.faa") for sample in self.sample_map}
        
        with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
            executor.map(lambda item: run_hmmsearch(*item), faa_files.items())
        logging.info("ORF annotation completed")
    
    def process_hmmsearch_results(self):
        logging.info("Processing HMMsearch results...")
        scores_dir = os.path.join(self.tmp_dir, "scours")
        os.makedirs(scores_dir, exist_ok=True)
        
        def process_pfam_file(sample, f):
            try:
                output_file = os.path.join(scores_dir, f"{sample}.pfam")
                
                with open(f, 'r') as pfam_file:
                    lines = [line.strip().split() for line in pfam_file if not line.startswith('#')]
                
                sorted_lines = sorted(lines, key=lambda x: float(x[4]), reverse=True)  # Sort by e-value column
                
                with open(output_file, 'w') as out:
                    for line in sorted_lines:
                        out.write("\t".join([line[0], line[2], line[4], line[5], line[8]]) + "\n")
                logging.info(f"Processed HMMsearch results for sample {sample}")
            except Exception as e:
                logging.error(f"Error processing HMMsearch results for sample {sample}: {str(e)}")
                raise
        
        pfam_files = {sample: os.path.join(self.hmmsearch_dir, f"{sample}.pfam") for sample in self.sample_map}
        
        with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
            executor.map(lambda item: process_pfam_file(*item), pfam_files.items())
        logging.info("HMMsearch results processing completed")
    
    def generate_final_output(self):
        logging.info("Generating final output...")
        output_path = os.path.join(self.output_dir, "annotation_output.tsv")
        results = []
        
        for sample in self.sample_map:
            pfam_file = os.path.join(self.hmmsearch_dir, f"{sample}.pfam")
            if os.path.exists(pfam_file):
                try:
                    with open(pfam_file, 'r') as f:
                        lines = [line.strip().split() for line in f]
                    
                    annotations = lines[:3] if len(lines) >= 3 else lines
                    
                    row = [sample]
                    for ann in annotations:
                        row.extend([ann[1], ann[4]])  # Gene ID, e-value
                    while len(row) < 10:
                        row.append('NA')  # Fill missing values
                    
                    results.append(row)
                except Exception as e:
                    logging.error(f"Error processing final output for sample {sample}: {str(e)}")
                    raise
        
        df = pd.DataFrame(results, columns=["Sample", "Gene_ID", "Top_Annotation", "E-value_1", "Second_Annotation", "E-value_2", "Third_Annotation", "E-value_3", "Gene_Coordinates"])
        df.to_csv(output_path, sep="\t", index=False)
        logging.info(f"Final output saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Call and annotate ORFs in assembled contigs using Prodigal and HMMsearch.")
    parser.add_argument("--assembly_dir", type=str, default="asm", help="Directory containing assembled contigs (default: asm)")
    parser.add_argument("--config_dir", type=str, default="configs", help="Directory containing configuration files (default: configs)")
    parser.add_argument("--output_dir", type=str, default="magus/annotation_output", help="Directory to store output files (default: magus/annotation_output)")
    parser.add_argument("--tmp_dir", type=str, default="tmp", help="Temporary directory for intermediate files (default: tmp)")
    parser.add_argument("--max_workers", type=int, default=8, help="Maximum number of parallel workers (default: 8)")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads per worker (default: 1)")
    parser.add_argument("--hmm_file", type=str, default=None, help="Path to HMM file for annotation. If not provided or not found, annotation will be skipped.")
    
    args = parser.parse_args()
    
    try:
        annotator = CallAnnotateORFs(
            assembly_dir=args.assembly_dir,
            config_dir=args.config_dir,
            output_dir=args.output_dir,
            tmp_dir=args.tmp_dir,
            num_workers=args.max_workers,
            threads=args.threads,
            hmm_file=args.hmm_file
        )
        annotator.call_orfs()
        annotator.filter_orfs()
        annotator.annotate_orfs()
        annotator.process_hmmsearch_results()
        annotator.generate_final_output()
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        raise

if __name__ == "__main__":
    main()
