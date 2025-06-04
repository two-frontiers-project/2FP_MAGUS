import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

class TestAssemble(unittest.TestCase):
    def test_main(self):
        # Your test code
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()

