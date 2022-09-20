import os
import unittest

import auto_fastq_symlink.samplesheet as ss

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class Test(unittest.TestCase):
    def test_parse_samplesheet_miseq_valid_01(self):
        samplesheet_path = os.path.join(THIS_DIR, "data", "samplesheets", "SampleSheet_miseq_v1_valid_01.csv")
        samplesheet = ss.parse_samplesheet_miseq(samplesheet_path)
        num_samples = len(samplesheet['data'])

        self.assertEqual(21, num_samples)

    def test_parse_samplesheet_nextseq_valid_01(self):
       samplesheet_path = os.path.join(THIS_DIR, "data", "samplesheets", "SampleSheet_nextseq_v1_valid_01.csv")
       samplesheet = ss.parse_samplesheet_nextseq(samplesheet_path)
       num_samples = len(samplesheet['cloud_data'])

       self.assertEqual(46, num_samples)

if __name__ == '__main__':
    unittest.main()
