#!/usr/bin/env python

import argparse
import glob
import os
import shutil


def main(args):
    """
    The run simulator only writes a SampleSheet at the top level of each run directory.
    Real runs also have a SampleSheet in each demultiplexing output directory, which is
    the only SampleSheet that auto-fastq-symlink uses for those directory structures.
    Copy the top-level SampleSheet to where the real demultiplexing output would put it.
    """
    run_dirs_glob = os.path.join(args.simulated_runs_dir, '*', '*', '*')
    run_dirs = glob.glob(run_dirs_glob)
    for run_dir in run_dirs:
        run_samplesheet = os.path.join(run_dir, 'SampleSheet.csv')
        if not os.path.exists(run_samplesheet):
            continue
        # MiSeq (new directory structure)
        demultiplexing_samplesheets = [os.path.join(d, 'SampleSheetUsed.csv') for d in glob.glob(os.path.join(run_dir, 'Alignment_*', '*'))]
        # NextSeq
        demultiplexing_samplesheets += [os.path.join(d, 'SampleSheet.csv') for d in glob.glob(os.path.join(run_dir, 'Analysis', '*', 'Data'))]
        for demultiplexing_samplesheet in demultiplexing_samplesheets:
            if not os.path.exists(demultiplexing_samplesheet):
                shutil.copyfile(run_samplesheet, demultiplexing_samplesheet)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--simulated-runs-dir')
    args = parser.parse_args()
    main(args)
