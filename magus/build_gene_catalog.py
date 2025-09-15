#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path

def parse_summary_file(summary_file, evalue_cutoff=0.01):
    annotated_genes = set()
    unannotated_genes = set()
    
    with open(summary_file, 'r') as f:
        headers = f.readline().strip().split('\t')
        query_col = 'query_accession' if 'query_accession' in headers else 'query_name'
        
        for line in f:
            parts = line.strip().split('\t')
            while len(parts) < len(headers):
                parts.append('')
            
            row = dict(zip(headers, parts))
            sample_id = row['sample_id']
            gene_id = row['sequence_id'] if 'sequence_id' in headers else row['target_name']
            query_name = row.get(query_col, '')
            full_evalue = row.get('full_evalue', '')
            
            full_gene_id = f"{sample_id}-----{gene_id}"
            
            if query_name and query_name.strip():
                try:
                    evalue = float(full_evalue)
                    if evalue <= evalue_cutoff:
                        annotated_genes.add(full_gene_id)
                    else:
                        unannotated_genes.add(full_gene_id)
                except:
                    unannotated_genes.add(full_gene_id)
            else:
                unannotated_genes.add(full_gene_id)
    
    print(f"Annotated genes: {len(annotated_genes)}")
    print(f"Unannotated genes: {len(unannotated_genes)}")
    return annotated_genes, unannotated_genes

def create_tocluster_faa(faa_dir, output_dir, unannotated_genes, identity_only=False):
    tocluster_faa = Path(output_dir) / "tocluster.faa"
    
    if identity_only:
        subprocess.run(f"cat {faa_dir}/*.faa > {tocluster_faa}", shell=True)
    else:
        gene_ids_file = Path(output_dir) / "gene_ids.txt"
        with open(gene_ids_file, 'w') as f:
            for gene_id in unannotated_genes:
                f.write(f"{gene_id}\n")
        
        with open(tocluster_faa, 'w') as out_f:
            for faa_file in Path(faa_dir).glob("*.faa"):
                cmd = f"grep -F -f {gene_ids_file} -A 1 {faa_file} | grep -v '^--$'"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                if result.stdout:
                    out_f.write(result.stdout)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--summary-file', required=True)
    parser.add_argument('--faa-dir', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--evalue-cutoff', type=float, default=0.01)
    parser.add_argument('--identity-threshold', type=float, default=0.9)
    parser.add_argument('--coverage-threshold', type=float, default=0.8)
    parser.add_argument('--identity-only', action='store_true')
    parser.add_argument('--multi-sample-single-copy', action='store_true')
    parser.add_argument('--tmp-dir', default='./tmp/')
    
    args = parser.parse_args()
    
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    annotated_genes, unannotated_genes = parse_summary_file(args.summary_file, args.evalue_cutoff)
    
    if args.identity_only:
        create_tocluster_faa(args.faa_dir, args.output_dir, unannotated_genes, identity_only=True)
    else:
        create_tocluster_faa(args.faa_dir, args.output_dir, unannotated_genes)

if __name__ == '__main__':
    main()