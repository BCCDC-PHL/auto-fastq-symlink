import logging

# Keep the application's logs out of the test output, without disabling them (so that tests can use `assertLogs`).
logging.getLogger('auto_fastq_symlink').addHandler(logging.NullHandler())
