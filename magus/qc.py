import subprocess
import os
from concurrent.futures import ThreadPoolExecutor

class QualityControl:
    def __init__(self, config, max_workers=4):
        self.config = config
        self.qc_dir = "qc2"
        self.max_workers = max_workers
        os.makedirs(self.qc_dir, exist_ok=True)

    def run_shi7_trimmer(self, sample):
        sample_name = sample['filename']
        r1 = sample['pe1']
        r2 = sample['pe2']
        cmd = f"shi7_trimmer {r1} {self.qc_dir}/{sample_name} 75 12 FLOOR 4 ASS_QUALITY 20 CASTN 0 STRIP ADAP2 CTGTCTCTTATACA R2 {r2} OUTFASTA"
        print(f"Running shi7_trimmer on {sample_name}")
        subprocess.run(cmd, shell=True)
        self.compress_output(sample_name)

    def compress_output(self, sample_name):
        for file in os.listdir(self.qc_dir):
            if file.startswith(sample_name):
                file_path = os.path.join(self.qc_dir, file)
                cmd = f"minigzip -4 {file_path}"
                print(f"Compressing {file_path}")
                subprocess.run(cmd, shell=True)

    def run(self):
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            executor.map(self.run_shi7_trimmer, self.config)
 
config = [
    {'filename': 'sample_1', 'pe1': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_1_R1.fq', 'pe2': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_1_R2.fq'},
    {'filename': 'sample_2', 'pe1': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_2_R1.fq', 'pe2': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_2_R2.fq'},
    {'filename': 'sample_3', 'pe1': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_3_R1.fq', 'pe2': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_3_R2.fq'},
    {'filename': 'sample_4', 'pe1': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_4_R1.fq', 'pe2': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_4_R2.fq'},
    {'filename': 'sample_5', 'pe1': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_5_R1.fq', 'pe2': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_5_R2.fq'},
    {'filename': 'sample_6', 'pe1': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_6_R1.fq', 'pe2': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_6_R2.fq'},
    {'filename': 'sample_7', 'pe1': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_7_R1.fq', 'pe2': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_7_R2.fq'},
    {'filename': 'sample_8', 'pe1': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_8_R1.fq', 'pe2': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_8_R2.fq'},
    {'filename': 'sample_9', 'pe1': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_9_R1.fq', 'pe2': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_9_R2.fq'},
    {'filename': 'sample_10', 'pe1': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_10_R1.fq', 'pe2': '/mnt/b/2FP_MAGUS/dev/syndata/synthetic_metagenomes/sample_10_R2.fq'}
]

# Example usage:
qc = QualityControl(config, max_workers=10)
qc.run()
