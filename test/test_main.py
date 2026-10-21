import json
import os
import shutil
import tempfile
import unittest
from unittest import mock

import auto_fastq_symlink.__main__ as main


class TestMain(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_reload_config_keeps_last_valid_config_on_failure(self):
        config_path = os.path.join(self.tmp_dir, 'config.json')
        with open(config_path, 'w') as f:
            f.write('{"run_parent_dirs": [')
        last_valid_config = {'run_parent_dirs': []}

        with self.assertLogs('auto_fastq_symlink.__main__', level='ERROR') as logs:
            config = main.reload_config(config_path, last_valid_config)

        self.assertIs(last_valid_config, config)
        self.assertEqual('load_config_failed', logs.records[0].msg['event_type'])

    def test_symlink_run_errors_are_logged_not_raised(self):
        run = {'run_id': '220602_M00123_0300_000000000-Q5539'}
        for error, expect_traceback in [(PermissionError('denied'), False), (KeyError('libraries'), True)]:
            with self.subTest(error=error):
                with mock.patch('auto_fastq_symlink.core.symlink_run', side_effect=error):
                    with self.assertLogs('auto_fastq_symlink.__main__', level='ERROR') as logs:
                        main.symlink_run({}, run)

                [record] = logs.records
                self.assertEqual('symlink_run_failed', record.msg['event_type'])
                self.assertEqual(type(error).__name__, record.msg['error_type'])
                self.assertEqual(expect_traceback, record.exc_info is not None)


if __name__ == '__main__':
    unittest.main()
