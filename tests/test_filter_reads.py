import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from magus import filter_reads


class TestParseConfig(unittest.TestCase):
    def test_parse_config(self):
        with tempfile.NamedTemporaryFile('w+', delete=False) as tmp:
            tmp.write('filename\tpe1\tpe2\n')
            tmp.write('sample1\treads_1.fq\treads_2.fq\n')
            tmp.write('sample2\tsingle.fq\n')
            tmp_path = tmp.name
        try:
            cfg = filter_reads.parse_config(tmp_path)
            self.assertIn('sample1', cfg)
            self.assertIn('sample2', cfg)
            self.assertEqual(cfg['sample1'], ('reads_1.fq', 'reads_2.fq'))
            self.assertEqual(cfg['sample2'], ('single.fq', None))
        finally:
            os.unlink(tmp_path)


if __name__ == '__main__':
    unittest.main()
