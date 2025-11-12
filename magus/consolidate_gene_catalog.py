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
    """
    Load annotations into memory, keyed by clustering rep extracted from gene name.
    Gene names are like: GCF_050908935.1_ASM5090893v1_genomic-----NZ_CP184181.1_592-----...
    We extract the clustering rep (NZ_CP184181.1_592) which is between the first two '-----'
    """
    annotations = {}
    with open(annotation_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row.get('gene', '')
            if gene:
                # Extract clustering rep from gene name (between first two '-----')
                # Format: sample-----clustering_rep-----rest
                parts = gene.split('-----')
                if len(parts) >= 2:
                    clustering_rep = parts[1]
                    # Key by clustering rep - if multiple annotations for same rep, last one wins
                    annotations[clustering_rep] = row
                else:
                    # Fallback: use gene name as-is if format doesn't match
                    annotations[gene] = row
    return annotations

def load_gene_metadata(metadata_file, gene_column='gene'):
    """
    Load additional gene metadata into memory, keyed by gene name.
    
    Args:
        metadata_file: Path to metadata TSV file
        gene_column: Name of the gene column in metadata file (default: 'gene')
    
    Returns:
        tuple: (metadata dict keyed by gene, metadata header list)
    """
    metadata = {}
    header = None
    
    with open(metadata_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        header = reader.fieldnames
        if not header:
            logger.warning(f"Metadata file {metadata_file} has no header")
            return metadata, []
        
        # Find gene column
        if gene_column not in header:
            # Try first column if 'gene' not found
            logger.warning(f"Gene column '{gene_column}' not found in metadata header. Using first column.")
            gene_column = header[0]
        
        for row in reader:
            gene = row.get(gene_column, '')
            if gene:
                metadata[gene] = row
    
    # Remove gene column from header (already in catalog)
    metadata_header = [col for col in header if col != gene_column]
    
    return metadata, metadata_header

def load_sample_metadata(metadata_file, sample_column='sample'):
    """
    Load additional sample metadata into memory, keyed by sample name.
    
    Args:
        metadata_file: Path to metadata TSV file
        sample_column: Name of the sample column in metadata file (default: 'sample')
    
    Returns:
        tuple: (metadata dict keyed by sample, metadata header list)
    """
    metadata = {}
    header = None
    
    with open(metadata_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        header = reader.fieldnames
        if not header:
            logger.warning(f"Metadata file {metadata_file} has no header")
            return metadata, []
        
        # Find sample column
        if sample_column not in header:
            # Try first column if 'sample' not found
            logger.warning(f"Sample column '{sample_column}' not found in metadata header. Using first column.")
            sample_column = header[0]
        
        for row in reader:
            sample = row.get(sample_column, '')
            if sample:
                metadata[sample] = row
    
    # Remove sample column from header (already in catalog)
    metadata_header = [col for col in header if col != sample_column]
    
    return metadata, metadata_header

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

def consolidate_catalog(catalog_file, annotation_file, output_file, merge_column=None, harmonize=False, metadata_file=None, sample_metadata_file=None):
    """
    Consolidate gene catalog with annotations.
    
    Args:
        catalog_file: Path to gene catalog TSV
        annotation_file: Path to merged annotations TSV
        output_file: Path to output consolidated catalog
        merge_column: Column name or threshold value to merge annotations into (default: last column)
        harmonize: If True, create a new harmonized column
        metadata_file: Path to additional gene metadata TSV file (optional)
        sample_metadata_file: Path to additional sample metadata TSV file (optional)
    """
    logger.info(f"Loading gene catalog from {catalog_file}")
    catalog, catalog_header = load_gene_catalog(catalog_file)
    
    logger.info(f"Loading annotations from {annotation_file}")
    annotations = load_annotations(annotation_file)
    logger.info(f"Loaded {len(annotations)} unique gene annotations")
    
    # Debug: show first few annotation clustering reps
    if annotations:
        sample_reps = list(annotations.keys())[:5]
        logger.info(f"Sample annotation clustering reps: {sample_reps}")
    
    # Load additional gene metadata if provided
    gene_metadata = {}
    gene_metadata_header = []
    if metadata_file:
        logger.info(f"Loading additional gene metadata from {metadata_file}")
        gene_metadata, gene_metadata_header = load_gene_metadata(metadata_file)
    
    # Load additional sample metadata if provided
    sample_metadata = {}
    sample_metadata_header = []
    if sample_metadata_file:
        logger.info(f"Loading additional sample metadata from {sample_metadata_file}")
        sample_metadata, sample_metadata_header = load_sample_metadata(sample_metadata_file)
    
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
    
    # Add harmonized column if requested (RIGHT AFTER annotations, before metadata)
    if harmonize:
        output_header.insert(insert_pos, 'harmonized_identity')
        insert_pos += 1
    
    # Add gene metadata columns after harmonized column
    if gene_metadata_header:
        for col in gene_metadata_header:
            output_header.insert(insert_pos, col)
            insert_pos += 1
    
    # Add sample metadata columns after gene metadata
    if sample_metadata_header:
        for col in sample_metadata_header:
            output_header.insert(insert_pos, col)
            insert_pos += 1
    
    # Process catalog entries
    output_rows = []
    annotated_count = 0
    sample_catalog_reps = []
    for key, catalog_row in catalog.items():
        sample, gene = key
        
        # Debug: collect sample clustering reps from merge column
        if len(sample_catalog_reps) < 5 and merge_col_idx < len(catalog_row):
            sample_catalog_reps.append(catalog_row[merge_col_idx])
        
        # Create output row starting with catalog row
        output_row = list(catalog_row)  # Make sure it's a list copy
        
        # Get clustering rep from merge column (user-specified OR last column)
        clustering_rep = catalog_row[merge_col_idx] if merge_col_idx < len(catalog_row) else ''
        
        # Get annotation if exists (left join on clustering rep from merge column)
        annotation = annotations.get(clustering_rep, None)
        if annotation:
            annotated_count += 1
        
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
        
        # Add harmonized column if requested (RIGHT AFTER annotations, before metadata)
        if harmonize:
            # Get clustering rep from merge column
            clustering_rep_val = catalog_row[merge_col_idx] if merge_col_idx < len(catalog_row) else ''
            
            # Use annotation if available (treating it as even higher resolution)
            if annotation and annotation.get('label'):
                harmonized_val = annotation.get('label')
            elif annotation and annotation.get('product_name'):
                harmonized_val = annotation.get('product_name')
            else:
                harmonized_val = clustering_rep_val
            output_row.insert(insert_pos, harmonized_val)
            insert_pos += 1
        
        # Insert gene metadata values after harmonized column (left join on root gene column, column 1)
        if gene_metadata_header:
            root_gene = catalog_row[1] if len(catalog_row) > 1 else ''
            metadata_row = gene_metadata.get(root_gene, {})
            for col in gene_metadata_header:
                metadata_val = metadata_row.get(col, '')
                output_row.insert(insert_pos, metadata_val)
                insert_pos += 1
        
        # Insert sample metadata values after gene metadata (left join on sample column, column 0)
        if sample_metadata_header:
            sample_name = catalog_row[0] if len(catalog_row) > 0 else ''
            sample_metadata_row = sample_metadata.get(sample_name, {})
            for col in sample_metadata_header:
                metadata_val = sample_metadata_row.get(col, '')
                output_row.insert(insert_pos, metadata_val)
                insert_pos += 1
        
        output_rows.append(output_row)
    
    # Debug: show sample catalog clustering reps
    if sample_catalog_reps:
        logger.info(f"Sample catalog clustering reps (merge column '{catalog_header[merge_col_idx]}'): {sample_catalog_reps}")
    
    logger.info(f"Processed {len(output_rows)} catalog entries, {annotated_count} had annotations")
    
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
    
    parser.add_argument(
        '--additional-gene-metadata',
        default=None,
        help='Path to additional gene metadata TSV file (one row per gene, left joined on root gene column)'
    )
    
    parser.add_argument(
        '--sample-metadata',
        default=None,
        help='Path to additional sample metadata TSV file (one row per sample, left joined on sample column)'
    )
    
    args = parser.parse_args()
    
    consolidate_catalog(
        catalog_file=args.gene_catalog,
        annotation_file=args.annotations,
        output_file=args.output,
        merge_column=args.annotation_merge_column,
        harmonize=args.harmonize_annotations,
        metadata_file=args.additional_gene_metadata,
        sample_metadata_file=args.sample_metadata
    )

if __name__ == '__main__':
    main()

