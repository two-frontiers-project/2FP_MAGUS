import subprocess
import os

class QualityControl:
    def __init__(self, config):
        self.config = config
        self.qc_dir = "qc"
        os.makedirs(self.qc_dir, exist_ok=True)

    def run_shi7_trimmer(self, sample):
        sample_name = sample['filename']
        r1 = sample['pe1']
        r2 = sample['pe2']
        cmd = f"shi7_trimmer {r1} {self.qc_dir}/{sample_name} 75 12 FLOOR 4 ASS_QUALITY 20 CASTN 0 STRIP ADAP2 CTGTCTCTTATACA R2 {r2} OUTFASTA"
        print(f"Running shi7_trimmer on {sample_name}")
        subprocess.run(cmd, shell=True)
    
    def compress_output(self, sample_name):
        # List all files in the 'qc' directory that start with the sample prefix
        for file in os.listdir(self.qc_dir):
            if file.startswith(sample_name):
                file_path = os.path.join(self.qc_dir, file)
                cmd = f"minigzip -4 {file_path}"
                print(f"Compressing {file_path}")
                subprocess.run(cmd, shell=True)
    
    def run(self):
        for sample in self.config:
            self.run_shi7_trimmer(sample)
            self.compress_output(sample['filename'])


# Example usage:
#config = [{'filename': 'SRR30713783', 'pe1': 'SRR30713783.sra_1.fastq.gz', 'pe2': 'SRR30713783.sra_2.fastq.gz'}]
#config = [{'filename': 'testsample', 'pe1': 'subsample_1.fastq', 'pe2': 'subsample_2.fastq'},{'filename': 'testsample2', 'pe1': 'subsample2_1.fastq', 'pe2': 'subsample2_2.fastq'}]
config = [
    {
        'filename': 'SRR15373729', 
        'pe1': '/mnt/b/random_testing_data/PRJNA751553/SRR15373729_1.fastq.gz', 
        'pe2': '/mnt/b/random_testing_data/PRJNA751553/SRR15373729_2.fastq.gz'
    },
    {
        'filename': 'SRR15373731', 
        'pe1': '/mnt/b/random_testing_data/PRJNA751553/SRR15373731_1.fastq.gz', 
        'pe2': '/mnt/b/random_testing_data/PRJNA751553/SRR15373731_2.fastq.gz'
    },
    {
        'filename': 'SRR15373733', 
        'pe1': '/mnt/b/random_testing_data/PRJNA751553/SRR15373733_1.fastq.gz', 
        'pe2': '/mnt/b/random_testing_data/PRJNA751553/SRR15373733_2.fastq.gz'
    }
]
qc = QualityControl(config)
qc.run()
