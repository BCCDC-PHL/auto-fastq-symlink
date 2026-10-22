import os
import shutil
import tempfile
import unittest

import jsonschema

import auto_fastq_symlink.samplesheet as ss


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


class TestParseSamplesheetNextseqSections(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.samplesheet_path = os.path.join(self.tmp_dir, 'SampleSheet.csv')

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def _parse(self, lines, encoding='utf-8'):
        with open(self.samplesheet_path, 'w', encoding=encoding) as f:
            f.write('\n'.join(lines) + '\n')
        return ss.parse_samplesheet(self.samplesheet_path, 'i100')

    def test_sections_not_separated_by_blank_lines(self):
        samplesheet = self._parse([
            '[Header]', 'FileFormatVersion,2',
            '[BCLConvert_Settings]', 'SoftwareVersion,4.3.13', 'FastqCompressionFormat,gzip',
            '[Cloud_Settings]', 'GeneratedVersion,1.0.0',
            '[Cloud_Data]', 'Sample_ID,ProjectName', 'LIB-A,routine_testing',
            '[Analysis_Settings]', 'Workflow,GenerateFASTQ',
        ])

        self.assertEqual({'software_version': '4.3.13', 'fastq_compression_format': 'gzip'}, samplesheet['bclconvert_settings'])
        self.assertEqual({'generated_version': '1.0.0'}, samplesheet['cloud_settings'])
        self.assertEqual([{'sample_id': 'LIB-A', 'project_name': 'routine_testing'}], samplesheet['cloud_data'])

    def test_only_cloud_data_required(self):
        samplesheet = self._parse(['[Cloud_Data]', 'Sample_ID,ProjectName', 'LIB-A,routine_testing'])

        self.assertEqual([], samplesheet['bclconvert_data'])
        self.assertEqual({}, samplesheet['bclconvert_settings'])
        self.assertEqual([{'sample_id': 'LIB-A', 'project_name': 'routine_testing'}], samplesheet['cloud_data'])

    def test_data_rows_with_empty_trailing_fields_padded(self):
        samplesheet = self._parse([
            '[BCLConvert_Data]', 'Sample_ID,Index,Index2,OverrideCycles,Lane', 'LIB-A,GAACTGAGCG,TCGTGGAGCG,,',
            '[Cloud_Data]', 'Sample_ID,ProjectName,LibraryName', 'LIB-A,routine_testing,', ',,',
        ])

        self.assertEqual([{'sample_id': 'LIB-A', 'index': 'GAACTGAGCG', 'index2': 'TCGTGGAGCG', 'override_cycles': '', 'lane': ''}], samplesheet['bclconvert_data'])
        self.assertEqual([{'sample_id': 'LIB-A', 'project_name': 'routine_testing', 'library_name': ''}], samplesheet['cloud_data'])

    def test_byte_order_mark_ignored(self):
        samplesheet = self._parse(['[Header]', 'FileFormatVersion,2', '[Cloud_Data]', 'Sample_ID,ProjectName', 'LIB-A,routine_testing'], encoding='utf-8-sig')

        self.assertEqual('2', samplesheet['header']['file_format_version'])

    def test_missing_cloud_data_project_name_throws(self):
        self.assertRaises(jsonschema.ValidationError, self._parse, ['[Cloud_Data]', 'Sample_ID,LibraryName', 'LIB-A,'])

if __name__ == '__main__':
    unittest.main()
