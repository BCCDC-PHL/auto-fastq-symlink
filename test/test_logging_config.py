import json
import logging
import os
import shutil
import tempfile
import time
import unittest

from auto_fastq_symlink.logging_config import JSONFormatter, make_log_handler


def _make_logger(handler):
    handler.setFormatter(JSONFormatter())
    logger = logging.getLogger('test_logging_config')
    logger.handlers = [handler]
    logger.propagate = False
    logger.setLevel(logging.INFO)
    return logger


def _read_log_entries(path):
    with open(path) as f:
        return [json.loads(line) for line in f]


def _read_events(path):
    return [entry['message']['event_type'] for entry in _read_log_entries(path)]


class TestLoggingConfig(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.log_file = os.path.join(self.tmp_dir, 'auto-fastq-symlink.jsonl')

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_log_entries_are_json_with_dict_message(self):
        handler = make_log_handler(self.log_file)
        _make_logger(handler).info({'event_type': 'scan_start', 'excluded_runs': {'run_01'}})
        handler.close()

        [entry] = _read_log_entries(self.log_file)
        self.assertEqual({'timestamp', 'level', 'module', 'function_name', 'line_num', 'message'}, set(entry))
        self.assertEqual('INFO', entry['level'])
        self.assertEqual({'event_type': 'scan_start', 'excluded_runs': "{'run_01'}"}, entry['message'])

    def test_exception_traceback_included(self):
        handler = make_log_handler(self.log_file)
        logger = _make_logger(handler)
        try:
            raise RuntimeError('boom')
        except RuntimeError:
            logger.error({'event_type': 'symlink_run_failed'}, exc_info=True)
        handler.close()

        [entry] = _read_log_entries(self.log_file)
        self.assertIn('RuntimeError: boom', entry['exception'])

    def test_log_file_is_appended_to_on_restart(self):
        for event_type in ['first_start', 'second_start']:
            handler = make_log_handler(self.log_file)
            _make_logger(handler).info({'event_type': event_type})
            handler.close()

        self.assertEqual(['first_start', 'second_start'], _read_events(self.log_file))

    def test_log_file_is_rotated_by_date(self):
        handler = make_log_handler(self.log_file, log_retention_days=2)
        logger = _make_logger(handler)

        logger.info({'event_type': 'yesterday'})
        yesterday = time.time() - 24 * 60 * 60
        os.utime(self.log_file, (yesterday, yesterday))
        handler.rolloverAt = time.time() - 1
        logger.info({'event_type': 'today'})
        handler.close()

        rotated_log_file = self.log_file + '.' + time.strftime('%Y-%m-%d', time.localtime(yesterday))
        self.assertEqual(['yesterday'], _read_events(rotated_log_file))
        self.assertEqual(['today'], _read_events(self.log_file))

    def test_old_rotated_log_files_are_deleted(self):
        for day in ['2026-10-01', '2026-10-02', '2026-10-03']:
            open(self.log_file + '.' + day, 'w').close()
        handler = make_log_handler(self.log_file, log_retention_days=2)
        logger = _make_logger(handler)

        logger.info({'event_type': 'before_rotation'})
        handler.rolloverAt = time.time() - 1
        logger.info({'event_type': 'after_rotation'})
        handler.close()

        rotated = sorted(f for f in os.listdir(self.tmp_dir) if f.startswith('auto-fastq-symlink.jsonl.'))
        self.assertEqual(2, len(rotated))
        self.assertNotIn('auto-fastq-symlink.jsonl.2026-10-01', rotated)


if __name__ == '__main__':
    unittest.main()
