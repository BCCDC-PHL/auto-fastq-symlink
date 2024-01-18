#!/usr/bin/env python

import argparse
import glob
import os
import random
import json


def main(args):
    run_dirs_glob = os.path.join(args.simulated_runs_dir, '*', '*', '*')
    run_dirs = glob.glob(run_dirs_glob)
    for run_dir in run_dirs:
        qc_check_complete_file = os.path.join(run_dir, 'qc_check_complete.json')
        if not(os.path.exists(qc_check_complete_file)):
            qc_check = {}
            if random.random() < args.proportion_failed:
                qc_check['overall_pass_fail'] = "FAIL"
            else:
                qc_check['overall_pass_fail'] = "PASS"

            with open(qc_check_complete_file, 'w') as f:
                f.write(json.dumps(qc_check, indent=2))
                f.write('\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--simulated-runs-dir')
    parser.add_argument('--proportion-failed', default=0.1, type=float)
    args = parser.parse_args()
    main(args)
