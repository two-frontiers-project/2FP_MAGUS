# magus_main.py
import sys
import subprocess

def main():
    if len(sys.argv) < 2:
        print("Usage: magus <command> [<args>]")
        sys.exit(1)

    command = sys.argv[1]
    args = sys.argv[2:]

    # Mapping commands to scripts
    command_map = {
        'qc': 'magus/qc.py',
        'cluster_contigs': 'magus/cluster_contigs.py',
        'coassembly': 'magus/coassembly.py',
        'coassembly_binning': 'magus/coassembly_binning.py',
        'run_checkv': 'magus/run_checkv.py',
        'run_eukrep': 'magus/run_eukrep.py',
        'select_final_bins': 'magus/select_final_bins.py',
        'single-assembly': 'magus/single-assembly.py',
        'single-binning': 'magus/single-binning.py'
    }

    # Check for valid command
    if command not in command_map:
        print(f"Unknown command: {command}")
        print("Available commands:", ", ".join(command_map.keys()))
        sys.exit(1)

    # Run the appropriate script
    script = command_map[command]
    subprocess.run(['python', script] + args)

if __name__ == '__main__':
    main()

