#!/usr/bin/env python3

import argparse
import csv
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_gene_catalog(catalog_file):
    """Load gene catalog into memory."""
    catalog = {}
    header = None
    with open(catalog_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        for row in reader:
            if len(row) < 3:
                continue
            sample = row[0]
            gene = row[1]
            key = (sample, gene)
            catalog[key] = row
    return catalog, header

def load_annotations(annotation_file):
    """Load annotations into memory, keyed by (sample_id, gene)."""
    annotations = {}
    with open(annotation_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sample_id = row.get('sample_id', '')
            gene = row.get('gene', '')
            key = (sample_id, gene)
            annotations[key] = row
    return annotations

def find_merge_column_index(header, merge_column):
    """Find the index of the merge column in the gene catalog header."""
    if merge_column is None:
        # Default: last column
        return len(header) - 1
    
    # Try to match by threshold value (e.g., "0.9")
    for idx, col in enumerate(header):
        if col == merge_column or col == str(merge_column):
            return idx
    
    # If not found, try to find column that contains the threshold
    for idx, col in enumerate(header):
        if merge_column in col or str(merge_column) in col:
            return idx
    
    logger.error(f"Could not find merge column '{merge_column}' in header: {header}")
    return None

def consolidate_catalog(catalog_file, annotation_file, output_file, merge_column=None, harmonize=False):
    """
    Consolidate gene catalog with annotations.
    
    Args:
        catalog_file: Path to gene catalog TSV
        annotation_file: Path to merged annotations TSV
        output_file: Path to output consolidated catalog
        merge_column: Column name or threshold value to merge annotations into (default: last column)
        harmonize: If True, create a new harmonized column
    """
    logger.info(f"Loading gene catalog from {catalog_file}")
    catalog, catalog_header = load_gene_catalog(catalog_file)
    
    logger.info(f"Loading annotations from {annotation_file}")
    annotations = load_annotations(annotation_file)
    
    # Find merge column index
    merge_col_idx = find_merge_column_index(catalog_header, merge_column)
    if merge_col_idx is None:
        logger.error("Failed to find merge column")
        return
    
    logger.info(f"Merging annotations into column '{catalog_header[merge_col_idx]}' (index {merge_col_idx})")
    
    # Build output header
    output_header = catalog_header.copy()
    
    # Add annotation columns after the merge column
    annotation_cols = ['match', 'evalue', 'seq_score', 'domain_score', 'label', 'gene_symbol', 
                       'product_name', 'ec_numbers', 'go_terms', 'for_AMRFinder', 'comment', 'annotation_source']
    
    # Insert annotation columns after merge column
    insert_pos = merge_col_idx + 1
    for col in annotation_cols:
        output_header.insert(insert_pos, col)
        insert_pos += 1
    
    # Add harmonized column if requested
    if harmonize:
        output_header.append('harmonized_identity')
    
    # Process catalog entries
    output_rows = []
    for key, catalog_row in catalog.items():
        sample, gene = key
        
        # Create output row starting with catalog row
        output_row = catalog_row.copy()
        
        # Get annotation if exists
        annotation = annotations.get(key, None)
        
        # Insert annotation columns at the right position
        if annotation:
            annotation_values = [
                annotation.get('match', ''),
                annotation.get('evalue', ''),
                annotation.get('seq_score', ''),
                annotation.get('domain_score', ''),
                annotation.get('label', ''),
                annotation.get('gene_symbol', ''),
                annotation.get('product_name', ''),
                annotation.get('ec_numbers', ''),
                annotation.get('go_terms', ''),
                annotation.get('for_AMRFinder', ''),
                annotation.get('comment', ''),
                annotation.get('annotation_source', '')
            ]
        else:
            annotation_values = [''] * len(annotation_cols)
        
        # Insert annotation values after merge column
        insert_pos = merge_col_idx + 1
        for val in annotation_values:
            output_row.insert(insert_pos, val)
            insert_pos += 1
        
        # Add harmonized column if requested
        if harmonize:
            # Find highest resolution clustering column (lowest threshold value)
            # Columns 0 and 1 are sample and gene, columns 2+ are clustering thresholds
            highest_res_rep = ''
            if len(catalog_row) > 2:
                # Start from the last column (highest resolution) and work backwards
                for idx in range(len(catalog_row) - 1, 1, -1):
                    if idx < len(catalog_row) and catalog_row[idx]:
                        highest_res_rep = catalog_row[idx]
                        break
            
            # Use annotation if available (treating it as even higher resolution)
            if annotation and annotation.get('label'):
                harmonized_val = annotation.get('label')
            elif annotation and annotation.get('product_name'):
                harmonized_val = annotation.get('product_name')
            else:
                harmonized_val = highest_res_rep
            output_row.append(harmonized_val)
        
        output_rows.append(output_row)
    
    # Write output
    logger.info(f"Writing consolidated catalog to {output_file}")
    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(output_header)
        writer.writerows(output_rows)
    
    logger.info(f"Consolidated catalog written with {len(output_rows)} entries")

def main():
    parser = argparse.ArgumentParser(
        description="Consolidate gene catalog with annotations"
    )
    
    parser.add_argument(
        '--gene-catalog',
        required=True,
        help='Path to gene catalog TSV file (from build-gene-catalog)'
    )
    
    parser.add_argument(
        '--annotations',
        required=True,
        help='Path to merged annotations TSV file (from annotate)'
    )
    
    parser.add_argument(
        '--output',
        required=True,
        help='Path to output consolidated catalog TSV file'
    )
    
    parser.add_argument(
        '--annotation-merge-column',
        default=None,
        help='Column name or threshold value to merge annotations into (default: last column in catalog)'
    )
    
    parser.add_argument(
        '--harmonize-annotations',
        action='store_true',
        help='Create a new harmonized column where annotations replace rep gene names'
    )
    
    args = parser.parse_args()
    
    consolidate_catalog(
        catalog_file=args.gene_catalog,
        annotation_file=args.annotations,
        output_file=args.output,
        merge_column=args.annotation_merge_column,
        harmonize=args.harmonize_annotations
    )

if __name__ == '__main__':
    main()

