import unittest
from unittest.mock import patch
import tempfile
import os
import pandas as pd

# Import the module using importlib since the file name contains a hyphen
import importlib.util
module_path = os.path.join(os.path.dirname(__file__), '..', 'magus', 'subsample-reads.py')
spec = importlib.util.spec_from_file_location('subsample_reads', module_path)
subsample_reads = importlib.util.module_from_spec(spec)
spec.loader.exec_module(subsample_reads)

class TestSubsampleReads(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        config_path = os.path.join(self.temp_dir.name, 'config.tsv')
        df = pd.DataFrame([
            ['sample1', 'read1.fq.gz', 'read2.fq.gz']
        ], columns=['filename', 'pe1', 'pe2'])
        df.to_csv(config_path, sep='\t', index=False)
        out_config = os.path.join(self.temp_dir.name, 'out_config.tsv')
        self.subsampler = subsample_reads.ReadSubsampler(
            config=config_path,
            outdir=self.temp_dir.name,
            out_config=out_config,
            depth=10
        )

    def tearDown(self):
        self.temp_dir.cleanup()

    @patch('subprocess.run')
    def test_paired_end_output_has_gz(self, mock_run):
        _, _, pe2_out = self.subsampler.subsample_reads('sample1', 'read1.fq.gz', 'read2.fq.gz')
        self.assertTrue(pe2_out.endswith('.gz'))

if __name__ == '__main__':
    unittest.main()
