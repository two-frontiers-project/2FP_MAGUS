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
        'cluster-contigs': 'magus/cluster-contigs.py',
        'coassembly': 'magus/coassembly.py',
        'coassembly-binning': 'magus/coassembly-binning.py',
        'find-viruses': 'magus/find-viruses.py',
        'find-euks': 'magus/find-euks.py',
        'finalize-bacterial-mags': 'magus/finalize-bacterial-mags.py',
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

