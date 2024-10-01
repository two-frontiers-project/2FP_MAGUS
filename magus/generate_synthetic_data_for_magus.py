import os
import subprocess
import random
import argparse
import gzip

def parse_config(config_file):
    """Parses the configuration file into a list of dictionaries."""
    genomes = []
    with open(config_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split(',')
                genome_info = {
                    'type': parts[0],
                    'identifier': parts[1],
                    'frequency': float(parts[2]),
                    'coverage': float(parts[3])  # Directly read coverage value
                }
                genomes.append(genome_info)
    return genomes

def select_genomes(genomes, num_samples):
    """Selects genomes for each sample based on frequency."""
    selected_genomes = []
    for _ in range(num_samples):
        sample_genomes = []
        for genome in genomes:
            if random.random() < genome['frequency']:
                sample_genomes.append(genome)
        selected_genomes.append(sample_genomes)
    return selected_genomes

def generate_1x_coverage(genome_path, sample_dir, genome_basename, read_length=150):
    """Generates 1x coverage reads for the genome using ART."""
    output_prefix = os.path.join(sample_dir, genome_basename)
    
    # Command for ART to simulate 1x coverage reads
    art_cmd = [
        'art_illumina',
        '-ss', 'HS25',  # HiSeq 2500 as an example
        '-i', genome_path,
        '-p', '-l', str(read_length),
        '-f', '1',  # Set coverage to 1x for now
        '-o', output_prefix
    ]
    
    # Run ART command to generate 1x coverage reads
    subprocess.run(art_cmd)
    
    # Return paths to generated read files
    return f"{output_prefix}1.fq", f"{output_prefix}2.fq"

def concat_reads(output_file, input_file, coverage):
    """Concatenates input file content into the output file for the specified coverage."""
    with open(output_file, 'ab') as outfile:
        for _ in range(int(coverage)):
            with open(input_file, 'rb') as infile:
                outfile.write(infile.read())

def generate_synthetic_reads(genomes, output_dir, num_samples, read_length=150):
    """Generates synthetic reads for each genome in each sample based on the specified coverage."""
    for i in range(num_samples):
        sample_genomes = genomes[i]
        sample_dir = os.path.join(output_dir, f'sample_{i+1}')
        os.makedirs(sample_dir, exist_ok=True)

        sample_r1 = os.path.join(sample_dir, f'sample_{i+1}_R1.fq')
        sample_r2 = os.path.join(sample_dir, f'sample_{i+1}_R2.fq')

        # Initialize empty files to concatenate reads into
        open(sample_r1, 'wb').close()
        open(sample_r2, 'wb').close()

        for genome in sample_genomes:
            genome_path = genome['identifier']
            coverage = genome['coverage']
            genome_basename = os.path.basename(genome_path).replace('.fasta.gz', '')

            # Generate 1x coverage reads
            r1_path, r2_path = generate_1x_coverage(genome_path, sample_dir, genome_basename, read_length)

            # Concatenate reads into the final sample file based on the required coverage
            concat_reads(sample_r1, r1_path, coverage)
            concat_reads(sample_r2, r2_path, coverage)

            # Optionally remove intermediate 1x coverage files
            os.remove(r1_path)
            os.remove(r2_path)

def main(config_file, output_dir, num_samples):
    # Parsing configuration file
    genomes = parse_config(config_file)

    # Select genomes for each sample based on frequency
    selected_genomes = select_genomes(genomes, num_samples)

    # Generate synthetic reads for each sample
    generate_synthetic_reads(selected_genomes, output_dir, num_samples)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Synthetic data generator for co-assembly pipeline.')
    parser.add_argument('config_file', type=str, help='Path to the configuration file.')
    parser.add_argument('output_dir', type=str, help='Directory to store synthetic reads.')
    parser.add_argument('num_samples', type=int, help='Number of samples to generate.')
    
    args = parser.parse_args()
    
    # Run the main function
    main(args.config_file, args.output_dir, args.num_samples)

