import os
import jsonschema
import logging
import unittest

import auto_fastq_symlink.samplesheet as ss

logging.disable(logging.CRITICAL)

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class Test(unittest.TestCase):
    def test_parse_samplesheet_miseq_valid_01_correct_num_samples(self):
        samplesheet_path = os.path.join(THIS_DIR, "data", "samplesheets", "SampleSheet_miseq_v1_valid_01.csv")
        samplesheet = ss.parse_samplesheet_miseq(samplesheet_path)
        num_samples = len(samplesheet['data'])

        self.assertEqual(21, num_samples)

    def test_parse_samplesheet_miseq_valid_01_all_samples_have_id(self):
        samplesheet_path = os.path.join(THIS_DIR, "data", "samplesheets", "SampleSheet_miseq_v1_valid_01.csv")
        samplesheet = ss.parse_samplesheet_miseq(samplesheet_path)
        libraries = samplesheet['data']
        has_sample_id = ['sample_id' in library for library in libraries]

        self.assertTrue(all(has_sample_id))


    def test_parse_samplesheet_miseq_valid_01_all_samples_have_project_id(self):
        samplesheet_path = os.path.join(THIS_DIR, "data", "samplesheets", "SampleSheet_miseq_v1_valid_01.csv")
        samplesheet = ss.parse_samplesheet_miseq(samplesheet_path)
        libraries = samplesheet['data']
        has_project_id = ['sample_project' in library for library in libraries]

        self.assertTrue(all(has_project_id))

    def test_parse_samplesheet_miseq_invalid_01_throws(self):
        samplesheet_path = os.path.join(THIS_DIR, "data", "samplesheets", "SampleSheet_miseq_v1_invalid_01.csv")

        self.assertRaises(jsonschema.ValidationError, ss.parse_samplesheet_miseq, samplesheet_path)

    def test_parse_samplesheet_nextseq_valid_01_correct_num_samples(self):
       samplesheet_path = os.path.join(THIS_DIR, "data", "samplesheets", "SampleSheet_nextseq_v1_valid_01.csv")
       samplesheet = ss.parse_samplesheet_nextseq(samplesheet_path)
       num_samples = len(samplesheet['cloud_data'])

       self.assertEqual(46, num_samples)

    def test_parse_samplesheet_nextseq_valid_01_all_samples_have_id(self):
        samplesheet_path = os.path.join(THIS_DIR, "data", "samplesheets", "SampleSheet_nextseq_v1_valid_01.csv")
        samplesheet = ss.parse_samplesheet_nextseq(samplesheet_path)
        libraries = samplesheet['cloud_data']
        has_sample_id = ['sample_id' in library for library in libraries]

        self.assertTrue(all(has_sample_id))

    def test_parse_samplesheet_nextseq_valid_01_all_samples_have_project_id(self):
        samplesheet_path = os.path.join(THIS_DIR, "data", "samplesheets", "SampleSheet_nextseq_v1_valid_01.csv")
        samplesheet = ss.parse_samplesheet_nextseq(samplesheet_path)
        libraries = samplesheet['cloud_data']
        has_project_id = ['project_name' in library for library in libraries]

        self.assertTrue(all(has_project_id))

    def test_parse_samplesheet_nextseq_invalid_01_throws(self):
        samplesheet_path = os.path.join(THIS_DIR, "data", "samplesheets", "SampleSheet_nextseq_v1_invalid_01.csv")

        self.assertRaises(jsonschema.ValidationError, ss.parse_samplesheet_nextseq, samplesheet_path)

if __name__ == '__main__':
    unittest.main()
