import json
import os
import shutil
import tempfile
import unittest

from auto_fastq_symlink.config import load_config


class TestLoadConfig(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.excluded_runs_path = os.path.join(self.tmp_dir, 'excluded_runs.csv')
        self.excluded_libraries_path = os.path.join(self.tmp_dir, 'excluded_libraries.csv')
        with open(self.excluded_runs_path, 'w') as f:
            f.write('220602_M00123_0300_000000000-Q5539\n')
        with open(self.excluded_libraries_path, 'w') as f:
            f.write('LIB-E\n')
        self.projects_path = os.path.join(self.tmp_dir, 'projects.csv')
        self._write_projects(self.excluded_runs_path, self.excluded_libraries_path)
        self.config = {
            'run_parent_dirs': [os.path.join(self.tmp_dir, 'runs')],
            'projects_definition_file': self.projects_path,
            'scan_interval_seconds': 60,
            'fastq_extensions': ['.fastq.gz'],
        }
        self.config_path = os.path.join(self.tmp_dir, 'config.json')

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def _write_projects(self, excluded_runs_path, excluded_libraries_path):
        with open(self.projects_path, 'w') as f:
            f.write('project_id,fastq_symlinks_dir,excluded_runs_list,excluded_libraries_list,simplify_symlink_filenames\n')
            f.write(','.join(['routine_testing', os.path.join(self.tmp_dir, 'symlinks'), excluded_runs_path, excluded_libraries_path, 'True']) + '\n')

    def _write_config(self):
        with open(self.config_path, 'w') as f:
            json.dump(self.config, f)

    def test_valid_config_loaded(self):
        self._write_config()

        config = load_config(self.config_path)

        project = config['projects']['routine_testing']
        self.assertEqual({'220602_M00123_0300_000000000-Q5539'}, project['excluded_runs'])
        self.assertEqual({'LIB-E'}, project['excluded_libraries'])
        self.assertTrue(project['simplify_symlink_filenames'])
        self.assertEqual({}, config['project_id_translation'])
        self.assertEqual(60.0, config['scan_interval_seconds'])

    def test_missing_excluded_runs_list_raises(self):
        self._write_projects(os.path.join(self.tmp_dir, 'missing.csv'), self.excluded_libraries_path)
        self._write_config()

        self.assertRaises(FileNotFoundError, load_config, self.config_path)

    def test_missing_excluded_libraries_list_raises(self):
        self._write_projects(self.excluded_runs_path, os.path.join(self.tmp_dir, 'missing.csv'))
        self._write_config()

        self.assertRaises(FileNotFoundError, load_config, self.config_path)

    def test_missing_required_key_raises(self):
        del self.config['fastq_extensions']
        self._write_config()

        with self.assertRaisesRegex(ValueError, 'fastq_extensions'):
            load_config(self.config_path)

    def test_invalid_scan_interval_raises(self):
        self.config['scan_interval_seconds'] = 'hourly'
        self._write_config()

        self.assertRaises(ValueError, load_config, self.config_path)

    def test_invalid_json_raises(self):
        with open(self.config_path, 'w') as f:
            f.write('{"run_parent_dirs": [')

        self.assertRaises(ValueError, load_config, self.config_path)


if __name__ == '__main__':
    unittest.main()
