import os
import glob
import argparse
import pandas as pd

def merge_alignments(input_dir, output_file, taxmap=None, coverage_cutoff=0.05):
    cov_files = glob.glob(os.path.join(input_dir, "*.cov"))

    if not cov_files:
        print(f"No .cov files found in {input_dir}")
        return

    dfs = []
    for file in sorted(cov_files):
        df = pd.read_csv(file, sep='\t')
        if df.empty:
            continue
        sample_name = os.path.splitext(os.path.basename(file))[0]
        df['sample'] = sample_name
        df['min_cov'] = df[['Adamantium_prop', 'Unique_proportion_covered', 'Proportion_covered']].min(axis=1)
        df['RA'] = df['Coverage_est'] / df['Coverage_est'].sum()
        dfs.append(df[['sample', 'Reference', 'min_cov', 'RA']])

    if not dfs:
        print("All .cov files were empty.")
        return

    merged = pd.concat(dfs, ignore_index=True)
    merged = merged[merged['min_cov'] >= coverage_cutoff]

    if merged.empty:
        print(f"No entries passed coverage cutoff of {coverage_cutoff}")
        return

    tax_path = taxmap
    if tax_path and os.path.exists(tax_path):
        tax_df = pd.read_csv(tax_path, sep='\t', header=None, names=['Reference', 'taxonomy'])
        # Strip RS_/GB_ prefixes in case the file has them
        tax_df['Reference'] = tax_df['Reference'].str.replace(r'^(RS_|GB_)', '', regex=True)
        merged = merged.merge(tax_df, on='Reference', how='left')
    elif tax_path:
        print(f"Warning: taxmap not found at {tax_path}, skipping taxonomy merge")

    col_order = ['sample', 'Reference', 'taxonomy', 'min_cov', 'RA'] if 'taxonomy' in merged.columns else ['sample', 'Reference', 'min_cov', 'RA']
    merged = merged[col_order].sort_values(['sample', 'RA'], ascending=[True, False])

    merged.to_csv(output_file, sep='\t', index=False)
    print(f"Merged alignments written to {output_file} ({len(merged)} rows, {merged['sample'].nunique()} samples)")


def parse_args():
    parser = argparse.ArgumentParser(description="Merge .cov files from xtree into a long-form abundance table.")
    parser.add_argument("--input-dir", required=True, help="Directory containing .cov files")
    parser.add_argument("--output", default="merged_alignments.tsv", help="Output TSV file path (default: merged_alignments.tsv)")
    parser.add_argument("--taxmap", default=None, help="Path to taxonomy TSV (no header; columns: Reference, taxonomy)")
    parser.add_argument("--coverage-cutoff", type=float, default=0.05, help="Minimum min_cov threshold (default 0.05)")
    return parser.parse_args()


def main():
    args = parse_args()
    merge_alignments(
        input_dir=args.input_dir,
        output_file=args.output,
        taxmap=args.taxmap,
        coverage_cutoff=args.coverage_cutoff,
    )


if __name__ == "__main__":
    main()
