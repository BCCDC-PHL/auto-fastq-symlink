#!/usr/bin/env python

import argparse
import datetime
import logging
import os
import signal
import threading

import auto_fastq_symlink.core as core
from auto_fastq_symlink.config import load_config
from auto_fastq_symlink.logging_config import configure_logging

log = logging.getLogger(__name__)


def reload_config(config_path: str, current_config: dict[str, object]) -> dict[str, object]:
    """
    Reload the config file, so that changes take effect without restarting.
    If the config can't be loaded, continue with the last valid config.

    :param config_path: Path to the config file.
    :type config_path: str
    :param current_config: The last valid config.
    :type current_config: dict[str, object]
    :return: The reloaded config, or the current config if reloading failed.
    :rtype: dict[str, object]
    """
    try:
        config = load_config(config_path)
        log.info({"event_type": "config_loaded", "config_file": os.path.abspath(config_path)})
        return config
    except Exception as e:
        log.error({"event_type": "load_config_failed", "config_file": os.path.abspath(config_path), "error_type": type(e).__name__, "error": str(e)})
        return current_config


def symlink_run(config: dict[str, object], run: dict[str, object]):
    """
    Create symlinks for one run. A problem with one run shouldn't stop other runs from being symlinked.
    Any run symlinks dir that wasn't completed is re-built on the next scan.

    :param config: Application config.
    :type config: dict[str, object]
    :param run: Sequencing run info.
    :type run: dict[str, object]
    """
    try:
        core.symlink_run(config, run)
    except OSError as e:
        # Expected, environmental problems (missing files, permissions, network storage).
        # The message says what went wrong, so a traceback wouldn't add anything.
        log.error({"event_type": "symlink_run_failed", "sequencing_run_id": run['run_id'], "error_type": type(e).__name__, "error": str(e)})
    except Exception as e:
        # Unexpected problems are probably bugs, so include the traceback.
        log.error({"event_type": "symlink_run_failed", "sequencing_run_id": run['run_id'], "error_type": type(e).__name__, "error": str(e)}, exc_info=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', required=True)
    parser.add_argument('-i', '--scan-interval', type=float, default=10, help='Seconds between scans, if not set by "scan_interval_seconds" in the config (default: 10)')
    parser.add_argument('--log-level', default='info', choices=['debug', 'info', 'warning', 'error'])
    parser.add_argument('--log-file', help='Append logs to this file, rotated daily at midnight (default: log to stdout)')
    parser.add_argument('--log-retention-days', type=int, default=90, help='Number of rotated log files to keep. 0 keeps all of them (default: 90)')
    args = parser.parse_args()

    configure_logging(args.log_level, args.log_file, args.log_retention_days)
    log.debug({"event_type": "debug_logging_enabled"})

    # Ctrl-c (SIGINT) or SIGTERM requests a stop, which happens at a safe time:
    # in between runs, or while waiting for the next scan. A second signal stops immediately.
    stop_requested = threading.Event()

    def request_stop(signum, frame):
        if stop_requested.is_set():
            log.warning({"event_type": "quit_immediately", "signal": signal.Signals(signum).name})
            raise SystemExit(1)
        log.info({"event_type": "quit_when_safe_enabled", "signal": signal.Signals(signum).name})
        stop_requested.set()

    signal.signal(signal.SIGINT, request_stop)
    signal.signal(signal.SIGTERM, request_stop)

    # Fail fast if the config is invalid at startup.
    config = load_config(args.config)
    log.info({"event_type": "config_loaded", "config_file": os.path.abspath(args.config)})

    while not stop_requested.is_set():
        scan_start_timestamp = datetime.datetime.now()
        for run in core.scan(config):
            symlink_run(config, run)
            if stop_requested.is_set():
                break
        if stop_requested.is_set():
            break

        scan_interval_seconds = config.get('scan_interval_seconds', args.scan_interval)
        scan_duration_seconds = (datetime.datetime.now() - scan_start_timestamp).total_seconds()
        next_scan_timestamp = scan_start_timestamp + datetime.timedelta(seconds=scan_interval_seconds)
        log.info({
            "event_type": "scan_complete",
            "scan_duration_seconds": scan_duration_seconds,
            "scan_interval_seconds": scan_interval_seconds,
            "timestamp_next_scan_start": next_scan_timestamp.isoformat(),
        })
        if stop_requested.wait(max(0, (next_scan_timestamp - datetime.datetime.now()).total_seconds())):
            break
        config = reload_config(args.config, config)

    log.info({"event_type": "quit"})


if __name__ == '__main__':
    main()
