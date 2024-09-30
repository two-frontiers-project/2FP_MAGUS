import os
import subprocess
import sys
import pandas as pd

class SequencingFileLoader:
    """Class to handle loading of sequencing files."""
    
    def __init__(self, seq_type, filepath1, filepath2=None):
        self.seq_type = seq_type  # 'single' or 'paired'
        self.filepath1 = filepath1
        self.filepath2 = filepath2
    
    def load(self):
        if not os.path.exists(self.filepath1):
            raise FileNotFoundError(f"Sequencing file {self.filepath1} not found.")
        
        if self.seq_type == 'paired':
            if not self.filepath2 or not os.path.exists(self.filepath2):
                raise FileNotFoundError(f"Paired-end file {self.filepath2} not found.")
            print(f"Loading paired-end sequencing files: {self.filepath1}, {self.filepath2}")
            return [self.filepath1, self.filepath2]
        else:
            print(f"Loading single-end sequencing file: {self.filepath1}")
            return self.filepath1

class KrakenTaxonomy:
    """Class to run Kraken taxonomy analysis."""
    
    def __init__(self, input_files, output_dir, taxonomic_level, config_file='config/kraken_config'):
        self.input_files = input_files  # single file or list of paired-end files
        self.output_dir = output_dir
        self.taxonomic_level = taxonomic_level  # e.g., 'phylum', 'species', 'genus'
        self.config = self._load_config(config_file)
        self.taxonomy_map = {
            'kingdom': 'K',
            'phylum': 'P',
            'class': 'C',
            'order': 'O',
            'family': 'F',
            'genus': 'G',
            'species': 'S'
        }
    
    def _load_config(self, config_file):
        config = {}
        with open(config_file, 'r') as f:
            for line in f:
                if not line.startswith("#") and line.strip():
                    key, value = line.strip().split('\t')
                    config[key] = value
        return config
    
    def run(self):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        kraken_db = self.config.get('kraken_database', 'default_db')
        confidence = self.config.get('kraken_confidence', 0.1)

        print(f"Running Kraken on {self.input_files} with confidence {confidence} using database {kraken_db}")
        
        output_file = os.path.join(self.output_dir, f"kraken_output.txt")
        report_file = os.path.join(self.output_dir, f"kraken_report.txt")

        command = [
            "kraken2",
            "--db", kraken_db,
            "--confidence", str(confidence),
            "--report", report_file,
            "--output", output_file,
            "--report-minimizer-data"
        ]

        # Add paired-end or single-end files
        if isinstance(self.input_files, list):
            command.extend(["--paired", self.input_files[0], self.input_files[1]])
        else:
            command.append(self.input_files)

        # Run Kraken
        print("Executing Kraken2 command: " + " ".join(command))
        subprocess.run(command, check=True)

        # Process the report file to aggregate counts at the requested taxonomic level
        return self._filter_report_by_taxonomic_level(report_file)
    
    def _filter_report_by_taxonomic_level(self, report_file):
        # Get the Kraken2 code for the requested taxonomic level (e.g., 'P' for Phylum)
        tax_level_code = self.taxonomy_map.get(self.taxonomic_level.lower())
        if tax_level_code is None:
            raise ValueError(f"Unsupported taxonomic level: {self.taxonomic_level}")

        # Data structure to hold counts for the taxonomic level
        tax_counts = {}

        # Open the Kraken2 report file
        with open(report_file, 'r') as report:
            for line in report:
                fields = line.strip().split('\t')
                
                # Ensure the line has at least 6 fields (matching the example you provided)
                if len(fields) < 6:
                    continue

                # Column 5 is the taxonomic rank, column 4 is the count
                taxonomic_rank = fields[5]
                
                # Only proceed if the current line matches the requested taxonomic rank
                if taxonomic_rank == tax_level_code:
                    read_count = int(fields[3])  # Column 4 contains the counts
                    
                    # Extract the taxon name by splitting the line after the NCBI taxid (column 6)
                    taxon_name = " ".join(fields[6:]).strip().split(' ')[-1]  # All remaining text after the taxid

                    # Store the counts in a dictionary keyed by the taxon name
                    if taxon_name not in tax_counts:
                        tax_counts[taxon_name] = read_count
                    else:
                        tax_counts[taxon_name] += read_count

        # Convert the counts into a DataFrame for wide format output
        tax_df = pd.DataFrame(list(tax_counts.items()), columns=['Taxon', 'Count'])
        tax_df.set_index('Taxon', inplace=True)

        # Save the wide-form matrix to a CSV file
        filtered_report = os.path.join(self.output_dir, f"kraken_report_{self.taxonomic_level}.csv")
        tax_df.to_csv(filtered_report)

        return filtered_report

class TaxonomyPipeline:
    """Main class to run the taxonomy pipeline."""
    
    def __init__(self, seq_type, output_dir, taxonomic_level, config_file, seq_file1, seq_file2=None):
        self.seq_type = seq_type  # 'single' or 'paired'
        self.output_dir = output_dir
        self.taxonomic_level = taxonomic_level
        self.config_file = config_file
        self.seq_file1 = seq_file1
        self.seq_file2 = seq_file2
    
    def run(self):
        loader = SequencingFileLoader(self.seq_type, self.seq_file1, self.seq_file2)
        input_file = loader.load()

        taxonomy = KrakenTaxonomy(input_file, self.output_dir, self.taxonomic_level, self.config_file)

        result_files = taxonomy.run()
        print(f"Taxonomy results saved to: {result_files}")
        return result_files


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: taxonomy.py <sequencing_type> <output_directory> <taxonomic_level> <config_file> <sequencing_file1> [<sequencing_file2>]")
        sys.exit(1)
    
    sequencing_type = sys.argv[1].lower()
    output_directory = sys.argv[2]
    taxonomic_level = sys.argv[3]
    config_file = sys.argv[4]
    sequencing_file1 = sys.argv[5]
    sequencing_file2 = sys.argv[6] if len(sys.argv) > 6 else None

    if sequencing_type not in ['single', 'paired']:
        print("Error: sequencing_type must be either 'single' or 'paired'")
        sys.exit(1)

    # Run the pipeline
    pipeline = TaxonomyPipeline(sequencing_type, output_directory, taxonomic_level, config_file, sequencing_file1, sequencing_file2)
    pipeline.run()
