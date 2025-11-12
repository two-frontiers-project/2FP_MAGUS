import sys
import subprocess
import os

def main():
    if len(sys.argv) < 2:
        print("Usage: magus <command> [<args>]")
        sys.exit(1)

    command = sys.argv[1]
    args = sys.argv[2:]

    # Mapping commands to scripts
    command_map = {
        'qc': 'qc.py',
        'assemble-hosts':'assemble-host.py',
        'subsample-reads':'subsample-reads.py',
        'filter-reads':'filter_reads.py',
        'taxonomy': 'taxonomy.py',
        'cluster-contigs': 'cluster-contigs.py',
        'single-assembly': 'single-assembly.py',
        'binning': 'single-binning.py',
        'cluster-contigs': 'cluster-contigs.py',
        'coassembly': 'coassembly.py',
        'coassembly-binning': 'coassembly-binning.py',
        'dereplicate': 'dereplicate-genomes.py',
        'find-viruses': 'find-viruses.py',
        'find-euks': 'find-euks.py',
        'finalize-bacterial-mags': 'finalize-bacterial-mags.py',
        'call-orfs': 'call_orfs2.py',
        'annotate': 'annotate.py',
        'build-gene-catalog': 'build_gene_catalog.py',
        'filter-mags': 'filter_mags.py',
        'build-tree': 'build_tree.py'
  }

    # Check for valid command
    if command not in command_map:
        print(f"Unknown command: {command}")
        print("Available commands:", ", ".join(command_map.keys()))
        sys.exit(1)

    # Get the absolute path to the script inside the package
    script = command_map[command]
    script_path = os.path.join(os.path.dirname(__file__), script)

    # Run the appropriate script
    subprocess.run(['python', script_path] + args)



if __name__ == '__main__':
    main()
