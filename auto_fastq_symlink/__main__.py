#!/usr/bin/env python

import argparse
import datetime
import json
import time
import logging
import os

import auto_fastq_symlink.config
import auto_fastq_symlink.core as core

DEFAULT_SCAN_INTERVAL_SECONDS = 3600.0

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
        format='{"timestamp": "%(asctime)s.%(msecs)03d", "level": "%(levelname)s", "module", "%(module)s", "function_name": "%(funcName)s", "line_num", %(lineno)d, "message": %(message)s}',
        datefmt='%Y-%m-%dT%H:%M:%S',
        encoding='utf-8',
        level=log_level,
    )
    logging.debug(json.dumps({"event_type": "debug_logging_enabled"}))

    scan_interval = args.scan_interval

    # We'll trap any KeyboardInterrupt and toggle this to True,
    # then exit at a safe time (in between runs or at the end of a scan of all runs)
    quit_when_safe = False

    while(True):
        if quit_when_safe:
            exit(0)

        try:
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
            scan_start_timestamp = datetime.datetime.now()
            for run in core.scan(config):
                if run is not None:
                    core.symlink_run(config, run)
                if quit_when_safe:
                    exit(0)
            scan_complete_timestamp = datetime.datetime.now()
            scan_duration_delta = scan_complete_timestamp - scan_start_timestamp
            scan_duration_seconds = scan_duration_delta.total_seconds()
            scan_interval_seconds = config.get('scan_interval_seconds', None)
            if scan_interval_seconds:
                next_scan_timestamp = datetime.datetime.now() + datetime.timedelta(seconds=scan_interval_seconds)
            else:
                next_scan_timestamp = datetime.datetime.now() + datetime.timedelta(seconds=DEFAULT_SCAN_INTERVAL_SECONDS)

            logging.info(json.dumps({"event_type": "scan_complete", "scan_duration_seconds": scan_duration_seconds, "timestamp_next_scan": str(next_scan_timestamp.isoformat())}))
            
            if quit_when_safe:
                exit(0)

            if "scan_interval_seconds" in config:
                try:
                    scan_interval = float(str(config['scan_interval_seconds']))
                except ValueError as e:
                    scan_interval = DEFAULT_SCAN_INTERVAL_SECONDS
            time.sleep(scan_interval)
        except KeyboardInterrupt as e:
            logging.info(json.dumps({"event_type": "quit_when_safe_enabled"}))
            quit_when_safe = True


if __name__ == '__main__':
    main()
