import os
import subprocess
import random
import csv
import gzip
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil  # New import for removing directories

def parse_config(config_file):
    genomes = []
    with open(config_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            genome_info = {
                'type': row['type'],
                'identifier': row['identifier'],
                'frequency': float(row['frequency']),
                'coverage_min': float(row['coverage_min']),
                'coverage_max': float(row['coverage_max'])
            }
            genomes.append(genome_info)
    return genomes

def generate_coverage(genome_path, sample_dir, genome_basename, coverage, read_length=150, fragment_length=300, std_dev=10):
    output_prefix = os.path.join(sample_dir, genome_basename)
    print(f"Generating {coverage}x coverage for {genome_basename} in {sample_dir}")

    dwgsim_cmd = [
        'dwgsim',
        '-e', '0.0001',
        '-E', '0.0001',
        '-1', str(read_length),
        '-2', str(read_length),
        '-C', str(coverage),
        '-r', '0',
        '-X', '0',
        '-d', str(fragment_length),
        '-s', str(std_dev),
        genome_path,
        output_prefix
    ]
    
    with open(os.devnull, 'w') as devnull:
        subprocess.run(dwgsim_cmd, stdout=devnull, stderr=devnull, check=True)

    r1_path = f"{output_prefix}.bwa.read1.fastq.gz"
    r2_path = f"{output_prefix}.bwa.read2.fastq.gz"
    return r1_path, r2_path

def assign_genomes_to_samples(genomes, num_samples):
    sample_compositions = [[] for _ in range(num_samples)]
    actual_coverage_data = []

    for genome in genomes:
        num_samples_to_appear = int(genome['frequency'] * num_samples)
        selected_samples = random.sample(range(num_samples), num_samples_to_appear)
        
        for sample_idx in selected_samples:
            coverage = random.uniform(genome['coverage_min'], genome['coverage_max'])
            sample_compositions[sample_idx].append((genome, coverage))
            
            actual_coverage_data.append({
                'type': genome['type'],
                'identifier': genome['identifier'],
                'sample': f'sample_{sample_idx+1}',
                'actual_coverage': coverage
            })
    
    return sample_compositions, actual_coverage_data

def generate_sample_reads(sample, sample_idx, output_dir, read_length=150):
    sample_dir = os.path.join(output_dir, f'sample_{sample_idx+1}')
    os.makedirs(sample_dir, exist_ok=True)

    sample_r1 = os.path.join(output_dir, f'sample_{sample_idx+1}_R1.fq')
    sample_r2 = os.path.join(output_dir, f'sample_{sample_idx+1}_R2.fq')

    # Initialize empty files for concatenation
    open(sample_r1, 'wb').close()
    open(sample_r2, 'wb').close()

    for genome, coverage in sample:
        genome_basename = os.path.basename(genome['identifier']).replace('.fna', '').replace('.gz', '')
        genome_r1, genome_r2 = generate_coverage(genome['identifier'], sample_dir, genome_basename, coverage, read_length)

        # Properly handle gzipped files by decompressing them on the fly
        if genome_r1.endswith('.fastq.gz') and genome_r2.endswith('.fastq.gz'):
            with gzip.open(genome_r1, 'rb') as infile_r1, open(sample_r1, 'ab') as outfile_r1:
                shutil.copyfileobj(infile_r1, outfile_r1)
            with gzip.open(genome_r2, 'rb') as infile_r2, open(sample_r2, 'ab') as outfile_r2:
                shutil.copyfileobj(infile_r2, outfile_r2)

    # Clean up the entire sample directory after concatenation
    shutil.rmtree(sample_dir)
    print(f"Cleaned up {sample_dir}")


def generate_synthetic_reads(sample_compositions, output_dir, read_length=150, max_workers=4):
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(generate_sample_reads, sample, idx, output_dir, read_length)
            for idx, sample in enumerate(sample_compositions)
        ]
        for future in as_completed(futures):
            future.result()

def write_actual_coverage(output_file, actual_coverage_data):
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['type', 'identifier', 'sample', 'actual_coverage']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(actual_coverage_data)

def main(config_file, output_dir, num_samples, actual_coverage_output, max_workers=4):
    genomes = parse_config(config_file)
    sample_compositions, actual_coverage_data = assign_genomes_to_samples(genomes, num_samples)
    generate_synthetic_reads(sample_compositions, output_dir, max_workers=max_workers)
    write_actual_coverage(actual_coverage_output, actual_coverage_data)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Synthetic data generation for co-assembly pipeline.')
    parser.add_argument('config_file', type=str, help='Path to the configuration file.')
    parser.add_argument('output_dir', type=str, help='Directory to store synthetic reads.')
    parser.add_argument('num_samples', type=int, help='Number of samples to generate.')
    parser.add_argument('actual_coverage_output', type=str, help='File path to store the actual coverage CSV output.')
    parser.add_argument('--max_workers', type=int, default=4, help='Maximum number of parallel workers.')

    args = parser.parse_args()
    main(args.config_file, args.output_dir, args.num_samples, args.actual_coverage_output, max_workers=args.max_workers)

