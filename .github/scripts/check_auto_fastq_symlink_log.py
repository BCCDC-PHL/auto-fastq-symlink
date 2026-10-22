#!/usr/bin/env python

import argparse
import json
import sys


def main(args):
    """
    Check the auto-fastq-symlink log from the integration test. Fail if anything was logged
    at ERROR level, or if any run whose demultiplexing output was found was never symlinked.

    Runs without demultiplexing output (eg. the last run, if the simulator was stopped part-way
    through writing it) are skipped by auto-fastq-symlink with a warning, and aren't expected to be symlinked.
    """
    with open(args.log) as f:
        entries = [json.loads(line) for line in f]

    errors = [entry for entry in entries if entry['level'] == 'ERROR']
    runs_found = {entry['message']['sequencing_run_id'] for entry in entries if entry['message'].get('event_type') == 'sequencing_run_found'}
    runs_symlinked = {entry['message']['sequencing_run_id'] for entry in entries if entry['message'].get('event_type') == 'symlink_run_complete'}
    runs_not_symlinked = sorted(runs_found - runs_symlinked)

    num_warnings = len([entry for entry in entries if entry['level'] == 'WARNING'])
    print(json.dumps({"num_runs_found": len(runs_found), "num_runs_symlinked": len(runs_symlinked), "num_errors": len(errors), "num_warnings": num_warnings}))
    for entry in errors[:10]:
        print(json.dumps(entry))
    for run_id in runs_not_symlinked:
        print(json.dumps({"run_not_symlinked": run_id}))

    if len(runs_found) == 0 or errors or runs_not_symlinked:
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--log')
    args = parser.parse_args()
    main(args)
