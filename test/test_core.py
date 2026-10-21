import json
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

import auto_fastq_symlink.core as core

logging.disable(logging.CRITICAL)

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLESHEETS_DIR = os.path.join(THIS_DIR, "data", "samplesheets")

MISEQ_RUN_ID = "220602_M00123_0300_000000000-Q5539"
NEXTSEQ_RUN_ID = "220602_VH00123_12_AAAAAAAAA"


def _write_miseq_samplesheet(path, libraries):
    """
    Write a MiSeq SampleSheet, using the header from a valid test SampleSheet.

    :param libraries: List of (library_id, project_id) tuples.
    """
    with open(os.path.join(SAMPLESHEETS_DIR, "SampleSheet_miseq_v1_valid_01.csv"), 'r') as f:
        lines = f.read().splitlines()
    data_header_index = lines.index('[Data]') + 1
    lines = lines[:data_header_index + 1]
    for library_id, project_id in libraries:
        lines.append(','.join([library_id, '', '', '', 'A01', 'UDP0001', 'GAACTGAGCG', 'UDP0001', 'TCGTGGAGCG', project_id, '']))
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def _write_nextseq_samplesheet(path, libraries):
    """
    Write a NextSeq SampleSheet, using the header from a valid test SampleSheet.

    :param libraries: List of (library_id, project_id) tuples.
    """
    with open(os.path.join(SAMPLESHEETS_DIR, "SampleSheet_nextseq_v1_valid_01.csv"), 'r') as f:
        lines = f.read().splitlines()
    bclconvert_data_index = lines.index('[BCLConvert_Data]')
    cloud_settings_index = lines.index('[Cloud_Settings]')
    cloud_data_index = lines.index('[Cloud_Data]')
    output_lines = lines[:bclconvert_data_index + 2]
    for library_id, _ in libraries:
        output_lines.append(','.join([library_id, 'GAACTGAGCG', 'TCGTGGAGCG']))
    output_lines.append('')
    output_lines += lines[cloud_settings_index:cloud_data_index + 2]
    for library_id, project_id in libraries:
        output_lines.append(','.join([library_id, project_id, '', '', '', '', '']))
    with open(path, 'w') as f:
        f.write('\n'.join(output_lines) + '\n')


def _write_fastqs(fastq_dir, library_ids):
    os.makedirs(fastq_dir, exist_ok=True)
    for sample_num, library_id in enumerate(library_ids, start=1):
        for read_type in ['R1', 'R2']:
            open(os.path.join(fastq_dir, library_id + '_S' + str(sample_num) + '_L001_' + read_type + '_001.fastq.gz'), 'w').close()


def _mark_run_ready(run_dir):
    open(os.path.join(run_dir, 'upload_complete.json'), 'w').close()
    with open(os.path.join(run_dir, 'qc_check_complete.json'), 'w') as f:
        json.dump({'overall_pass_fail': 'PASS'}, f)


class TestCore(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.runs_dir = os.path.join(self.tmp_dir, 'runs')
        self.symlinks_dir = os.path.join(self.tmp_dir, 'symlinks')
        os.makedirs(self.runs_dir)
        self.config = {
            'run_parent_dirs': [self.runs_dir],
            'fastq_extensions': ['.fastq.gz'],
            'project_id_translation': {},
            'projects': {
                'routine_testing': {
                    'fastq_symlinks_dir': os.path.join(self.symlinks_dir, 'routine_testing'),
                    'simplify_symlink_filenames': True,
                    'excluded_runs': set(),
                    'excluded_libraries': set(),
                },
                'assay_development': {
                    'fastq_symlinks_dir': os.path.join(self.symlinks_dir, 'assay_development'),
                    'simplify_symlink_filenames': False,
                    'excluded_runs': set(),
                    'excluded_libraries': set(),
                },
            },
        }

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def _scan_and_symlink(self):
        for run in core.scan(self.config):
            if run is not None:
                core.symlink_run(self.config, run)

    def _run_symlinks_dir(self, project_id, run_id):
        return os.path.join(self.config['projects'][project_id]['fastq_symlinks_dir'], run_id)

    def _symlink_targets(self, project_id, run_id):
        run_symlinks_dir = self._run_symlinks_dir(project_id, run_id)
        targets = {}
        for entry in os.scandir(run_symlinks_dir):
            if entry.is_symlink():
                targets[entry.name] = os.readlink(entry.path)
        return targets

    def _add_miseq_alignment(self, run_dir, alignment_num, timestamp, libraries):
        output_dir = os.path.join(run_dir, 'Alignment_' + str(alignment_num), timestamp)
        _write_fastqs(os.path.join(output_dir, 'Fastq'), [library_id for library_id, _ in libraries])
        _write_miseq_samplesheet(os.path.join(output_dir, 'SampleSheetUsed.csv'), libraries)
        return output_dir

    def _create_miseq_run(self, libraries):
        run_dir = os.path.join(self.runs_dir, MISEQ_RUN_ID)
        os.makedirs(run_dir)
        _write_miseq_samplesheet(os.path.join(run_dir, 'SampleSheet.csv'), libraries)
        _mark_run_ready(run_dir)
        return run_dir

    def test_miseq_v1_symlinks_created(self):
        libraries = [('LIB-A', 'routine_testing'), ('LIB-B', 'assay_development'), ('LIB-C', 'other')]
        run_dir = self._create_miseq_run(libraries)
        fastq_dir = os.path.join(run_dir, 'Data', 'Intensities', 'BaseCalls')
        _write_fastqs(fastq_dir, [library_id for library_id, _ in libraries])

        self._scan_and_symlink()

        self.assertEqual(
            {
                'LIB-A_R1.fastq.gz': os.path.join(fastq_dir, 'LIB-A_S1_L001_R1_001.fastq.gz'),
                'LIB-A_R2.fastq.gz': os.path.join(fastq_dir, 'LIB-A_S1_L001_R2_001.fastq.gz'),
            },
            self._symlink_targets('routine_testing', MISEQ_RUN_ID),
        )
        self.assertEqual(
            {
                'LIB-B_S2_L001_R1_001.fastq.gz': os.path.join(fastq_dir, 'LIB-B_S2_L001_R1_001.fastq.gz'),
                'LIB-B_S2_L001_R2_001.fastq.gz': os.path.join(fastq_dir, 'LIB-B_S2_L001_R2_001.fastq.gz'),
            },
            self._symlink_targets('assay_development', MISEQ_RUN_ID),
        )
        with open(os.path.join(self._run_symlinks_dir('routine_testing', MISEQ_RUN_ID), 'symlinks_complete.json')) as f:
            self.assertEqual(2, json.load(f)['num_symlinks_created'])
        self.assertFalse(os.path.exists(os.path.join(self.config['projects']['routine_testing']['fastq_symlinks_dir'], '.' + MISEQ_RUN_ID + '.tmp')))

    def test_existing_run_symlinks_dir_not_modified_after_redemultiplexing(self):
        libraries = [('LIB-A', 'routine_testing'), ('LIB-B', 'assay_development')]
        run_dir = self._create_miseq_run(libraries)
        first_output_dir = self._add_miseq_alignment(run_dir, 1, '20220603_101010', libraries)
        self._scan_and_symlink()
        targets_before = self._symlink_targets('routine_testing', MISEQ_RUN_ID)
        symlinks_complete_path = os.path.join(self._run_symlinks_dir('routine_testing', MISEQ_RUN_ID), 'symlinks_complete.json')
        with open(symlinks_complete_path) as f:
            symlinks_complete_before = f.read()

        new_libraries = [('LIB-A', 'routine_testing'), ('LIB-B', 'routine_testing'), ('LIB-D', 'routine_testing')]
        self._add_miseq_alignment(run_dir, 2, '20220610_101010', new_libraries)
        self._scan_and_symlink()

        self.assertEqual(targets_before, self._symlink_targets('routine_testing', MISEQ_RUN_ID))
        self.assertTrue(all(target.startswith(first_output_dir) for target in targets_before.values()))
        with open(symlinks_complete_path) as f:
            self.assertEqual(symlinks_complete_before, f.read())

    def test_deleted_run_symlinks_dir_recreated_from_latest_demultiplexing(self):
        libraries = [('LIB-A', 'routine_testing'), ('LIB-B', 'assay_development')]
        run_dir = self._create_miseq_run(libraries)
        self._add_miseq_alignment(run_dir, 1, '20220603_101010', libraries)
        self._scan_and_symlink()

        # New demultiplexing, with a SampleSheet that moves LIB-B to routine_testing and adds LIB-D
        new_libraries = [('LIB-A', 'routine_testing'), ('LIB-B', 'routine_testing'), ('LIB-D', 'routine_testing')]
        new_output_dir = self._add_miseq_alignment(run_dir, 2, '20220610_101010', new_libraries)
        shutil.rmtree(self._run_symlinks_dir('routine_testing', MISEQ_RUN_ID))
        self._scan_and_symlink()

        new_fastq_dir = os.path.join(new_output_dir, 'Fastq')
        self.assertEqual(
            {
                'LIB-A_R1.fastq.gz': os.path.join(new_fastq_dir, 'LIB-A_S1_L001_R1_001.fastq.gz'),
                'LIB-A_R2.fastq.gz': os.path.join(new_fastq_dir, 'LIB-A_S1_L001_R2_001.fastq.gz'),
                'LIB-B_R1.fastq.gz': os.path.join(new_fastq_dir, 'LIB-B_S2_L001_R1_001.fastq.gz'),
                'LIB-B_R2.fastq.gz': os.path.join(new_fastq_dir, 'LIB-B_S2_L001_R2_001.fastq.gz'),
                'LIB-D_R1.fastq.gz': os.path.join(new_fastq_dir, 'LIB-D_S3_L001_R1_001.fastq.gz'),
                'LIB-D_R2.fastq.gz': os.path.join(new_fastq_dir, 'LIB-D_S3_L001_R2_001.fastq.gz'),
            },
            self._symlink_targets('routine_testing', MISEQ_RUN_ID),
        )

    def test_miseq_latest_alignment_subdir_chosen(self):
        libraries = [('LIB-A', 'routine_testing')]
        run_dir = self._create_miseq_run(libraries)
        self._add_miseq_alignment(run_dir, 2, '20220603_101010', libraries)
        self._add_miseq_alignment(run_dir, 10, '20220603_101010', libraries)
        latest_output_dir = self._add_miseq_alignment(run_dir, 10, '20220610_101010', libraries)

        latest_demultiplexing = core._find_latest_demultiplexing(run_dir, 'miseq')

        self.assertEqual(latest_output_dir, latest_demultiplexing['output_directory'])
        self.assertEqual(os.path.join(latest_output_dir, 'Fastq'), latest_demultiplexing['fastq_directory'])

    def test_miseq_samplesheet_used_from_latest_alignment_preferred(self):
        libraries = [('LIB-A', 'routine_testing')]
        run_dir = self._create_miseq_run(libraries)
        output_dir = self._add_miseq_alignment(run_dir, 1, '20220603_101010', libraries)

        self.assertEqual(os.path.join(output_dir, 'SampleSheetUsed.csv'), core._find_samplesheet(run_dir, output_dir, 'miseq'))

    def test_miseq_samplesheet_used_preferred_over_competing_samplesheets_in_alignment_dir(self):
        libraries = [('LIB-A', 'routine_testing')]
        run_dir = self._create_miseq_run(libraries)
        output_dir = self._add_miseq_alignment(run_dir, 1, '20220603_101010', libraries)
        for competing_samplesheet in ['SampleSheet.csv', 'SampleSheet_old.csv']:
            _write_miseq_samplesheet(os.path.join(output_dir, competing_samplesheet), libraries)

        self.assertEqual(os.path.join(output_dir, 'SampleSheetUsed.csv'), core._find_samplesheet(run_dir, output_dir, 'miseq'))

    def test_miseq_alignment_dir_samplesheet_other_than_samplesheet_used_ignored(self):
        libraries = [('LIB-A', 'routine_testing')]
        run_dir = self._create_miseq_run(libraries)
        output_dir = self._add_miseq_alignment(run_dir, 1, '20220603_101010', libraries)
        os.rename(os.path.join(output_dir, 'SampleSheetUsed.csv'), os.path.join(output_dir, 'SampleSheet.csv'))

        self.assertEqual(os.path.join(run_dir, 'SampleSheet.csv'), core._find_samplesheet(run_dir, output_dir, 'miseq'))

    def test_miseq_falls_back_to_run_level_samplesheet(self):
        libraries = [('LIB-A', 'routine_testing')]
        run_dir = self._create_miseq_run(libraries)
        output_dir = self._add_miseq_alignment(run_dir, 1, '20220603_101010', libraries)
        os.remove(os.path.join(output_dir, 'SampleSheetUsed.csv'))

        self.assertEqual(os.path.join(run_dir, 'SampleSheet.csv'), core._find_samplesheet(run_dir, output_dir, 'miseq'))

    def test_nextseq_samplesheet_selection_unchanged(self):
        run_dir = os.path.join(self.runs_dir, NEXTSEQ_RUN_ID)
        output_dir = os.path.join(run_dir, 'Analysis', '1', 'Data')
        os.makedirs(output_dir)
        _write_nextseq_samplesheet(os.path.join(run_dir, 'SampleSheet.csv'), [('LIB-A', 'routine_testing')])
        _write_nextseq_samplesheet(os.path.join(output_dir, 'SampleSheet.csv'), [('LIB-A', 'routine_testing')])
        _write_nextseq_samplesheet(os.path.join(output_dir, 'SampleSheet_old.csv'), [('LIB-A', 'routine_testing')])

        self.assertEqual(os.path.join(output_dir, 'SampleSheet.csv'), core._find_samplesheet(run_dir, output_dir, 'nextseq'))

    def test_nextseq_latest_analysis_chosen_numerically(self):
        run_dir = os.path.join(self.runs_dir, NEXTSEQ_RUN_ID)
        os.makedirs(run_dir)
        _mark_run_ready(run_dir)
        for analysis_num, project_id in [(9, 'assay_development'), (10, 'routine_testing')]:
            output_dir = os.path.join(run_dir, 'Analysis', str(analysis_num), 'Data')
            _write_fastqs(os.path.join(output_dir, 'fastq'), ['LIB-A'])
            _write_nextseq_samplesheet(os.path.join(output_dir, 'SampleSheet.csv'), [('LIB-A', project_id)])

        self._scan_and_symlink()

        latest_fastq_dir = os.path.join(run_dir, 'Analysis', '10', 'Data', 'fastq')
        self.assertEqual(
            {
                'LIB-A_R1.fastq.gz': os.path.join(latest_fastq_dir, 'LIB-A_S1_L001_R1_001.fastq.gz'),
                'LIB-A_R2.fastq.gz': os.path.join(latest_fastq_dir, 'LIB-A_S1_L001_R2_001.fastq.gz'),
            },
            self._symlink_targets('routine_testing', NEXTSEQ_RUN_ID),
        )
        self.assertFalse(os.path.exists(self._run_symlinks_dir('assay_development', NEXTSEQ_RUN_ID)))

    def test_excluded_runs_and_libraries_not_symlinked(self):
        libraries = [('LIB-A', 'routine_testing'), ('LIB-E', 'routine_testing'), ('LIB-B', 'assay_development')]
        run_dir = self._create_miseq_run(libraries)
        self._add_miseq_alignment(run_dir, 1, '20220603_101010', libraries)
        self.config['projects']['routine_testing']['excluded_libraries'] = {'LIB-E'}
        self.config['projects']['assay_development']['excluded_runs'] = {MISEQ_RUN_ID}

        self._scan_and_symlink()

        self.assertEqual({'LIB-A_R1.fastq.gz', 'LIB-A_R2.fastq.gz'}, set(self._symlink_targets('routine_testing', MISEQ_RUN_ID)))
        self.assertFalse(os.path.exists(self._run_symlinks_dir('assay_development', MISEQ_RUN_ID)))

    def test_symlinks_complete_written_only_after_rename(self):
        libraries = [('LIB-A', 'routine_testing')]
        run_dir = self._create_miseq_run(libraries)
        self._add_miseq_alignment(run_dir, 1, '20220603_101010', libraries)
        tmp_run_symlinks_dir = os.path.join(self.config['projects']['routine_testing']['fastq_symlinks_dir'], '.' + MISEQ_RUN_ID + '.tmp')
        tmp_dir_contents_at_rename = []
        real_rename = os.rename

        def recording_rename(src, dst):
            if src == tmp_run_symlinks_dir:
                tmp_dir_contents_at_rename.extend(os.listdir(src))
            real_rename(src, dst)

        with mock.patch('auto_fastq_symlink.core.os.rename', side_effect=recording_rename):
            self._scan_and_symlink()

        self.assertEqual({'LIB-A_R1.fastq.gz', 'LIB-A_R2.fastq.gz'}, set(tmp_dir_contents_at_rename))
        self.assertEqual(
            {'LIB-A_R1.fastq.gz', 'LIB-A_R2.fastq.gz', 'symlinks_complete.json'},
            set(os.listdir(self._run_symlinks_dir('routine_testing', MISEQ_RUN_ID))),
        )

    def test_incomplete_tmp_dir_replaced(self):
        libraries = [('LIB-A', 'routine_testing')]
        run_dir = self._create_miseq_run(libraries)
        self._add_miseq_alignment(run_dir, 1, '20220603_101010', libraries)
        tmp_run_symlinks_dir = os.path.join(self.config['projects']['routine_testing']['fastq_symlinks_dir'], '.' + MISEQ_RUN_ID + '.tmp')
        os.makedirs(tmp_run_symlinks_dir)
        os.symlink('/nonexistent/stale_R1.fastq.gz', os.path.join(tmp_run_symlinks_dir, 'stale_R1.fastq.gz'))

        self._scan_and_symlink()

        self.assertFalse(os.path.lexists(tmp_run_symlinks_dir))
        self.assertEqual({'LIB-A_R1.fastq.gz', 'LIB-A_R2.fastq.gz'}, set(self._symlink_targets('routine_testing', MISEQ_RUN_ID)))


if __name__ == '__main__':
    unittest.main()
