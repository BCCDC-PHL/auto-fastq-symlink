#!/usr/bin/env python

import argparse
import json
import time
import logging
import os

import auto_fastq_symlink.config
import auto_fastq_symlink.core as core


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config')
    parser.add_argument('-i', '--scan-interval', default=10)
    parser.add_argument('--log-level', default="info")
    args = parser.parse_args()
    config = {}

    try:
        log_level = getattr(logging, args.log_level.upper())
    except AttributeError as e:
        log_level = logging.INFO

    logging.basicConfig(
        format='{"timestamp": "%(asctime)s.%(msecs)03d", "log_level": "%(levelname)s", "module", "%(module)s", "function_name": "%(funcName)s", "line_num", %(lineno)d, "message": %(message)s}',
        datefmt='%Y-%m-%dT%H:%M:%S',
        encoding='utf-8',
        level=log_level,
    )
    logging.debug(json.dumps({"event_type": "debug_logging_enabled"}))

    scan_interval = args.scan_interval

    # We'll trap any KeyboardInterrupt and toggle this to True,
    # then exit at a safe time (in between scan & symlink passes).
    quit_when_safe = False

    while(True):
        try:
            if quit_when_safe:
                exit(0)

            if args.config:
                logging.info(json.dumps({"event_type": "load_config_start", "config_file": os.path.abspath(args.config)}))
                try:
                    config = auto_fastq_symlink.config.load_config(args.config)
                    # Uncomment below to see the config on stdout each time it's reloaded
                    # print(json.dumps(auto_fastq_symlink.config.make_config_json_serializable(config), indent=2))
                except json.decoder.JSONDecodeError as e:
                    # If we fail to load the config file, we continue on with the
                    # last valid config that was loaded.
                    logging.error(json.dumps({"event_type": "load_config_failed", "config_file": os.path.abspath(args.config)}))

            # All of the action happens here.
            for run in core.scan(config):
                core.symlink_run(config, run)
                if quit_when_safe:
                    exit(0)            

            if "scan_interval_seconds" in config:
                scan_interval = config['scan_interval_seconds']
            time.sleep(scan_interval)
        except KeyboardInterrupt as e:
            logging.info(json.dumps({"event_type": "quit_when_safe_enabled"}))
            quit_when_safe = True


if __name__ == '__main__':
    main()
