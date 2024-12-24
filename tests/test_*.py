import unittest
from src.qc.quality_control import QualityControl
import os

class TestQualityControl(unittest.TestCase):
    def setUp(self):
        self.qc = QualityControl('test_output')
        os.makedirs('test_output', exist_ok=True)

    def test_run_fastqc(self):
        test_file = 'test_data/test.fastq'
        self.qc.run_fastqc([test_file])
        self.assertTrue(os.path.exists('test_output/qc_results/test_fastqc.html'))

    def tearDown(self):
        import shutil
        shutil.rmtree('test_output')

if __name__ == '__main__':
    unittest.main()
