import datetime
import json
import logging
import os
import re
import shutil
from typing import Iterable, Optional

import jsonschema

import auto_fastq_symlink.samplesheet as ss

log = logging.getLogger(__name__)


def _subdir_names(dir_path: str) -> list[str]:
    """
    :param dir_path: Path to a directory.
    :type dir_path: str
    :return: Names of the sub-directories of the directory.
    :rtype: list[str]
    """
    return [entry.name for entry in os.scandir(dir_path) if entry.is_dir()]


def _find_latest_analysis_dir(run_dir_path: str) -> Optional[str]:
    """
    NextSeq and MiSeq i100 runs have a numbered `Analysis/<N>` directory for each demultiplexing.
    Find the one with the highest number.

    :param run_dir_path: Path to the sequencing run directory.
    :type run_dir_path: str
    :return: Path to the latest analysis directory, or None if there isn't one.
    :rtype: str | None
    """
    analysis_dir_path = os.path.join(run_dir_path, "Analysis")
    if not os.path.isdir(analysis_dir_path):
        return None
    analysis_subdir_nums = [int(name) for name in _subdir_names(analysis_dir_path) if re.fullmatch("\\d+", name)]
    if len(analysis_subdir_nums) == 0:
        return None

    return os.path.join(analysis_dir_path, str(max(analysis_subdir_nums)))


def _find_latest_demultiplexing(run_dir_path: str, instrument_type: str) -> Optional[dict[str, Optional[str]]]:
    """
    Find the output of the latest demultiplexing for a run. A run may have been demultiplexed
    more than once, in which case each demultiplexing has its own numbered output directory,
    with its own SampleSheet.

    Only the SampleSheet for the demultiplexing is used. There's no fall-back to the SampleSheet at the
    top level of the run directory, because it may not match the latest demultiplexing.

    | Directory structure | Fastq directory                         | SampleSheet                                 |
    |---------------------|-----------------------------------------|---------------------------------------------|
    | MiSeq (old)         | `Data/Intensities/BaseCalls`            | `SampleSheet*.csv` at the top level         |
    | MiSeq (new)         | `Alignment_<N>/<timestamp>/Fastq`       | `Alignment_<N>/<timestamp>/SampleSheetUsed.csv` |
    | NextSeq             | `Analysis/<N>/Data/fastq`               | `Analysis/<N>/Data/SampleSheet*.csv`        |
    | MiSeq i100          | `Analysis/<N>/Data/BCLConvert/fastq`    | `Analysis/<N>/inputs/SampleSheet*.csv`      |

    The old MiSeq directory structure only supports one demultiplexing, so its SampleSheet is at the top level of the run directory.

    :param run_dir_path: Path to the sequencing run directory.
    :type run_dir_path: str
    :param instrument_type: Instrument type ('miseq', 'nextseq' or 'i100')
    :type instrument_type: str
    :return: Dict with keys 'output_directory', 'fastq_directory' and 'samplesheet_path' (which is None if the SampleSheet can't be found), or None if the demultiplexing output can't be found.
    :rtype: dict[str, str | None] | None
    """
    find_samplesheet = lambda dir_path: ss.choose_samplesheet_to_parse(ss.find_samplesheets(dir_path))
    latest_demultiplexing = None
    if instrument_type == 'miseq':
        alignment_subdir_nums = [int(name.split('_')[1]) for name in _subdir_names(run_dir_path) if re.fullmatch("Alignment_\\d+", name)]
        if len(alignment_subdir_nums) > 0:
            alignment_directory = os.path.join(run_dir_path, "Alignment_" + str(max(alignment_subdir_nums)))
            alignment_subdirs = sorted(_subdir_names(alignment_directory), reverse=True)
            if len(alignment_subdirs) > 0:
                output_directory = os.path.join(alignment_directory, alignment_subdirs[0])
                samplesheet_used_path = os.path.join(output_directory, "SampleSheetUsed.csv")
                latest_demultiplexing = {
                    'output_directory': output_directory,
                    'fastq_directory': os.path.join(output_directory, "Fastq"),
                    'samplesheet_path': samplesheet_used_path if os.path.isfile(samplesheet_used_path) else None,
                }
        else:
            data_intensities_basecalls_dir_path = os.path.join(run_dir_path, "Data", "Intensities", "BaseCalls")
            if os.path.exists(data_intensities_basecalls_dir_path):
                latest_demultiplexing = {
                    'output_directory': run_dir_path,
                    'fastq_directory': data_intensities_basecalls_dir_path,
                    'samplesheet_path': find_samplesheet(run_dir_path),
                }

    elif instrument_type == 'nextseq':
        analysis_dir_path = _find_latest_analysis_dir(run_dir_path)
        if analysis_dir_path is not None:
            latest_demultiplexing = {
                'output_directory': analysis_dir_path,
                'fastq_directory': os.path.join(analysis_dir_path, "Data", "fastq"),
                'samplesheet_path': find_samplesheet(os.path.join(analysis_dir_path, "Data")),
            }

    elif instrument_type == 'i100':
        analysis_dir_path = _find_latest_analysis_dir(run_dir_path)
        if analysis_dir_path is not None:
            latest_demultiplexing = {
                'output_directory': analysis_dir_path,
                'fastq_directory': os.path.join(analysis_dir_path, "Data", "BCLConvert", "fastq"),
                'samplesheet_path': find_samplesheet(os.path.join(analysis_dir_path, "inputs")),
            }

    return latest_demultiplexing


def _determine_libraries_section(samplesheet, instrument_type):
    """
    Determine which section of the SampleSheet contains the library info.

    :param samplesheet: The parsed SampleSheet.
    :type samplesheet: dict[str, object]
    :param instrument_type: Instrument type ('miseq', 'nextseq' or 'i100')
    :type instrument_type: str
    :return: The section of the SampleSheet that contains the library info.
    :rtype: str | None
    """
    if instrument_type == 'miseq':
        libraries_section = 'data'
    elif instrument_type in ('nextseq', 'i100'):
        libraries_section = 'cloud_data'
    else:
        libraries_section = None

    return libraries_section


def _determine_project_header(samplesheet, instrument_type):
    """
    For MiSeq runs, the SampleSheet has two fields that can either be used for the project ID.

    :param samplesheet: The parsed SampleSheet.
    :type samplesheet: dict[str, object]
    :param instrument_type: Instrument type ('miseq', 'nextseq' or 'i100')
    :type instrument_type: str
    :return: The header of the field that contains the project ID.
    :rtype: str | None
    """
    if instrument_type == 'miseq':
        project_header = 'sample_project'
    elif instrument_type in ('nextseq', 'i100'):
        project_header = 'project_name'
    else:
        project_header = None

    return project_header


def _determine_library_id_header(samplesheet, instrument_type):
    """
    For MiSeq runs, the SampleSheet has two fields that can either be used for the library ID.
    Either 'Sample_ID' or 'Sample_Name' may be used to input the library ID. The other field
    might be blank, or it might be filled with 'S1', 'S2', 'S3', etc.
    In order to find the correct field to use for the library ID, we need to look through those
    columns to check which was actually used for the librar ID on this run.

    :param samplesheet: The parsed SampleSheet.
    :type samplesheet: dict[str, object]
    :param instrument_type: Instrument type ('miseq', 'nextseq' or 'i100')
    :type instrument_type: str
    :return: The header of the field that contains the library ID.
    :rtype: str | None
    """
    libraries_section = _determine_libraries_section(samplesheet, instrument_type)
    if instrument_type == 'miseq':
        sample_ids = [library['sample_id'] for library in samplesheet[libraries_section]]
        all_sample_ids_blank = all([sample_id == "" for sample_id in sample_ids])
        all_sample_ids_s_plus_digits = all([re.match("S\\d+", sample_id) for sample_id in sample_ids])
        all_sample_ids_only_digits = all([re.match("\\d+", sample_id) for sample_id in sample_ids])
        if all_sample_ids_blank or all_sample_ids_s_plus_digits or all_sample_ids_only_digits:
            library_id_header = 'sample_name'
        else:
            library_id_header = 'sample_id'
    elif instrument_type in ('nextseq', 'i100'):
        library_id_header = 'sample_id'
    else:
        library_id_header = None

    return library_id_header


def _sanitize_library_id(library_id: str):
    """
    Sanitize library ID by replacing invalid characters with '-'.

    :param library_id: The library ID.
    :type library_id: str
    :return: The sanitized library ID.
    :rtype: str
    """
    sanitized_library_id = library_id.strip()
    sanitized_library_id = sanitized_library_id.replace('_', '-') # There shouldn't be underscores in library IDs, but if there are, they get swapped with '-' in fastq filenames.
    sanitized_library_id = sanitized_library_id.replace(' ', '-') # There shouldn't be spaces in library IDs, but if there are, they get swapped with '-' in fastq filenames.
    symbols = "\\`*{}[]()>#!$^"
    for s in symbols:
        sanitized_library_id = sanitized_library_id.replace(s, '-')
    sanitized_library_id = re.sub('-+', '-', sanitized_library_id) #

    return sanitized_library_id


def find_libraries(run: dict[str, object], samplesheet: dict[str, object], fastq_extensions: list[str]) -> list[dict[str, object]]:
    """
    Use parsed samplesheet and run info to find all libraries on the run, along with their project ID and fastq paths.

    :param run: The run info.
    :type run: dict[str, object]
    :param samplesheet: The parsed SampleSheet.
    :type samplesheet: dict[str, object]
    :param fastq_extensions: A list of valid fastq filename extensions (defined in config)
    :type fastq_extensions: list[str]
    :return: Libraries
    :rtype: list[dict[str, object]]
    """
    run_id = run['run_id']
    libraries = []
    # If we can't confirm that the fastq directory exists, no sense continuing.
    # Short-circuit and return an empty list.
    fastq_dir = None
    if 'fastq_directory' not in run or not run['fastq_directory']:
        log.error({"event_type": "find_libraries_failed", "sequencing_run_id": run_id})
        return libraries
    else:
        fastq_dir = run['fastq_directory']

    if not os.path.exists(fastq_dir):
        log.error({"event_type": "find_libraries_failed", "sequencing_run_id": run_id})
        return libraries

    libraries_section = _determine_libraries_section(samplesheet, run['instrument_type'])
    project_header = _determine_project_header(samplesheet, run['instrument_type'])

    has_correct_extension = lambda x: any([x.endswith(ext) for ext in fastq_extensions])
    run_fastq_dir_contents = os.listdir(fastq_dir)
    run_fastq_files = set(filter(lambda x: all([has_correct_extension(x), os.path.isfile(os.path.join(fastq_dir, x))]), run_fastq_dir_contents))

    library_id_header = _determine_library_id_header(samplesheet, run['instrument_type'])
    
    found_library_ids = set()
    for item in samplesheet[libraries_section]:
        required_item_keys = [project_header, library_id_header]
        if not all([k in item for k in required_item_keys]):
            log.warning({"event_type": "library_missing_required_fields", "required_fields": required_item_keys, "library": item})
            continue
        library = {}
        library_id = _sanitize_library_id(item[library_id_header])
        project_id = item[project_header]
        if library_id == "" or library_id in found_library_ids:
            continue
        library['library_id'] = library_id
        library['project_id'] = project_id
        log.debug({"event_type": "found_library", "library_id": library_id, "samplesheet_project_id": project_id})
        r1_fastq_filenames = list(filter(lambda x: re.match(library_id + '_S\\d+(_L\\d{3})?' + '_R1_' + '\\d{3}', x), run_fastq_files))
        if len(r1_fastq_filenames) > 0:
            r1_fastq_filename = r1_fastq_filenames[0]
        else:
            log.debug({"event_type": "failed_to_find_R1_fastq", "sequencing_run_id": run_id, "library_id": library_id, "samplesheet_project_id": project_id})
            r1_fastq_filename = None
        r2_fastq_filenames = list(filter(lambda x: re.match(library_id + '_S\\d+(_L\\d{3})?' + '_R2_' + '\\d{3}', x), run_fastq_files))
        if len(r2_fastq_filenames) > 0:
            r2_fastq_filename = r2_fastq_filenames[0]
        else:
            log.debug({"event_type": "failed_to_find_R2_fastq", "sequencing_run_id": run_id, "library_id": library_id, "samplesheet_project_id": project_id})
            r2_fastq_filename = None
        r1_fastq_path = os.path.join(fastq_dir, r1_fastq_filename) if r1_fastq_filename else None
        r2_fastq_path = os.path.join(fastq_dir, r2_fastq_filename) if r2_fastq_filename else None
        if r1_fastq_path and os.path.exists(r1_fastq_path):
            log.debug({"event_type": "found_library_fastq_file", "library_id": library_id, "read_type": "R1", "fastq_path": r1_fastq_path})
            library['fastq_path_r1'] = r1_fastq_path
        else:
            library['fastq_path_r1'] = None
        if r2_fastq_path and os.path.exists(r2_fastq_path):
            log.debug({"event_type": "found_library_fastq_file", "library_id": library_id, "read_type": "R2", "fastq_path": r2_fastq_path})
            library['fastq_path_r2'] = r2_fastq_path
        else:
            library['fastq_path_r2'] = None
        found_library_ids.add(library_id)
        libraries.append(library)
        

    return libraries
    

MISEQ_RUN_ID_REGEX = "\\d{6}_M\\d{5}_\\d+_\\d{9}-[A-Z0-9]{5}"
NEXTSEQ_RUN_ID_REGEX = "\\d{6}_VH\\d{5}_\\d+_[A-Z0-9]{9}"
I100_RUN_ID_REGEX = "\\d{8}_SH\\d{5}_\\d+_[A-Z0-9]{10}-[A-Z0-9]{3}"


def _load_run(config: dict[str, object], run_dir_path: str) -> Optional[dict[str, object]]:
    """
    Load a sequencing run, if it is ready to be symlinked. A run is ready once it has been uploaded and has passed its QC check.

    :param config: Application config.
    :type config: dict[str, object]
    :param run_dir_path: Path to the (possible) sequencing run directory.
    :type run_dir_path: str
    :return: Sequencing run info, including its libraries, or None if the directory isn't a run that is ready to be symlinked.
    :rtype: dict[str, object] | None
    """
    run_id = os.path.basename(run_dir_path)
    instrument_type = None
    if re.match(MISEQ_RUN_ID_REGEX, run_id):
        instrument_type = "miseq"
    elif re.match(NEXTSEQ_RUN_ID_REGEX, run_id):
        instrument_type = "nextseq"
    elif re.match(I100_RUN_ID_REGEX, run_id):
        instrument_type = "i100"

    qc_check_complete_path = os.path.join(run_dir_path, "qc_check_complete.json")
    qc_check_complete_file_exists = os.path.exists(qc_check_complete_path)
    passed_run_qc_check = False
    if qc_check_complete_file_exists:
        try:
            with open(qc_check_complete_path) as f:
                qc_check = json.load(f)
        except json.decoder.JSONDecodeError as e:
            log.error({"event_type": "qc_check_json_decode_error", "sequencing_run_id": run_id, "error": str(e)})
            qc_check = {}
        overall_pass_fail = qc_check.get('overall_pass_fail', None)
        if overall_pass_fail is not None and re.match("PASS", overall_pass_fail, re.IGNORECASE):
            passed_run_qc_check = True

    conditions_checked = {
        'subdir_is_directory': os.path.isdir(run_dir_path),
        'upload_complete': os.path.exists(os.path.join(run_dir_path, "upload_complete.json")),
        'qc_check_complete': qc_check_complete_file_exists,
        'passed_run_qc_check': passed_run_qc_check,
        'determined_instrument_type': instrument_type is not None,
    }
    if not all(conditions_checked.values()):
        log.info({"event_type": "skipped_run", "sequencing_run_id": run_id, "conditions_checked": conditions_checked})
        return None

    log.info({"event_type": "scan_run_start", "sequencing_run_id": run_id})
    latest_demultiplexing = _find_latest_demultiplexing(run_dir_path, instrument_type)
    if latest_demultiplexing is None:
        log.warning({"event_type": "demultiplexing_output_not_found", "sequencing_run_id": run_id})
        return None

    log.debug({"event_type": "sequencing_run_found", "sequencing_run_id": run_id})
    run = {
        "run_id": run_id,
        "instrument_type": instrument_type,
        "run_directory": run_dir_path,
        "fastq_directory": latest_demultiplexing['fastq_directory'],
    }
    samplesheet_to_parse = latest_demultiplexing['samplesheet_path']
    if samplesheet_to_parse:
        log.debug({"event_type": "samplesheet_found", "sequencing_run_id": run_id, "samplesheet_path": samplesheet_to_parse})
    else:
        log.error({"event_type": "samplesheet_not_found", "sequencing_run_id": run_id, "demultiplexing_output_directory": latest_demultiplexing['output_directory']})
        return None

    run['samplesheet'] = samplesheet_to_parse
    try:
        samplesheet = ss.parse_samplesheet(samplesheet_to_parse, run['instrument_type'])
    except jsonschema.ValidationError:
        # Already logged by the samplesheet parser.
        return None

    libraries = find_libraries(run, samplesheet, config['fastq_extensions'])
    for library in libraries:
        if library['project_id'] in config['project_id_translation']:
            samplesheet_project_id = library['project_id']
            symlinking_project_id = config['project_id_translation'][samplesheet_project_id]
            library['project_id'] = symlinking_project_id
        elif library['project_id'] == '':
            library['project_id'] = None
    run['libraries'] = libraries

    return run


def find_runs(config: dict[str, object]) -> Iterable[dict[str, object]]:
    """
    Find all sequencing runs under all of the `run_parent_dirs` from the config that are ready to be symlinked.
    Runs are found by matching sub-directory names against `MISEQ_RUN_ID_REGEX`, `NEXTSEQ_RUN_ID_REGEX` and `I100_RUN_ID_REGEX`.

    A problem with one run directory (or one run parent dir) is logged, and doesn't stop other runs from being found.
    The run will be tried again on the next scan.

    :param config: Application config.
    :type config: dict[str, object]
    :return: Sequencing runs that are ready to be symlinked.
    :rtype: Iterable[dict[str, object]]
    """
    for run_parent_dir in config['run_parent_dirs']:
        try:
            subdirs = list(os.scandir(run_parent_dir))
        except OSError as e:
            # Parent dirs may be on network storage that's temporarily unavailable.
            # Skip it for this scan, and try again on the next one.
            log.error({"event_type": "failed_to_scan_run_parent_dir", "run_parent_dir": run_parent_dir, "error": str(e)})
            continue

        for subdir in subdirs:
            try:
                run = _load_run(config, subdir.path)
            except OSError as e:
                # Expected, environmental problems (missing files, permissions, network storage).
                # The message says what went wrong, so a traceback wouldn't add anything.
                log.error({"event_type": "load_run_failed", "sequencing_run_id": subdir.name, "error_type": type(e).__name__, "error": str(e)})
                continue
            except Exception as e:
                # Unexpected problems are probably bugs, so include the traceback.
                log.error({"event_type": "load_run_failed", "sequencing_run_id": subdir.name, "error_type": type(e).__name__, "error": str(e)}, exc_info=True)
                continue

            if run is not None:
                yield run


def _symlink_filename(target: str, simplify_symlink_filenames: bool) -> str:
    """
    Determine the filename for a symlink to a fastq file. If `simplify_symlink_filenames` is set,
    the filename is reduced to the library ID followed by `_R1.fastq.gz` or `_R2.fastq.gz`.

    :param target: Path to the fastq file that the symlink will point to.
    :type target: str
    :param simplify_symlink_filenames: Whether to simplify the symlink filename.
    :type simplify_symlink_filenames: bool
    :return: Symlink filename
    :rtype: str
    """
    target_basename = os.path.basename(target)
    if not simplify_symlink_filenames:
        return target_basename

    r1_r2_match = re.search("_(R[12])_", target_basename)
    if r1_r2_match:
        r1_r2 = r1_r2_match.group(1)
    else:
        r1_r2 = ""
    symlink_filename = target_basename.split('_')[0] + '_' + r1_r2 + '.fastq.gz'

    return symlink_filename


def determine_symlinks_to_create_for_run(config: dict[str, object], run: dict[str, object]) -> dict[str, list[dict[str, str]]]:
    """
    Determine which symlinks should be created for a run, for each project.

    A project's symlinks directory for a run (`<fastq_symlinks_dir>/<run_id>`) is only ever created once.
    If it already exists, nothing is created for that project, even if the run has since been re-demultiplexed.
    To re-generate the symlinks for a run, delete (or move) that directory and they will be re-created on
    the next scan, using the latest demultiplexing of the run.

    :param config: Application config
    :type config: dict[str, object]
    :param run: Sequencing run info, including libraries.
    :type run: dict[str, object]
    :return: Symlinks to create, indexed by project ID. Each symlink has keys 'target' and 'path' (relative to the run's symlinks directory).
    :rtype: dict[str, list[dict[str, str]]]
    """
    run_id = run['run_id']
    log.debug({"event_type": "determine_symlinks_for_run_start", "sequencing_run_id": run_id})
    symlinks_to_create_by_project_id = {}

    for project_id, project in config['projects'].items():
        if run_id in project['excluded_runs']:
            continue

        project_libraries = [
            library for library in run['libraries']
            if library['project_id'] == project_id and library['library_id'] not in project['excluded_libraries']
        ]
        if len(project_libraries) == 0:
            continue

        run_symlinks_dir = os.path.join(project['fastq_symlinks_dir'], run_id)
        if os.path.lexists(run_symlinks_dir):
            if not os.path.exists(os.path.join(run_symlinks_dir, 'symlinks_complete.json')):
                # Most likely the application was interrupted between creating the directory and writing
                # symlinks_complete.json. Leave it alone, but make it visible. Delete the directory to re-generate it.
                log.warning({"event_type": "run_symlinks_dir_missing_symlinks_complete", "sequencing_run_id": run_id, "project_id": project_id, "run_symlinks_dir": run_symlinks_dir})
            else:
                log.debug({"event_type": "run_symlinks_dir_already_exists", "sequencing_run_id": run_id, "project_id": project_id, "run_symlinks_dir": run_symlinks_dir})
            continue

        symlinks = []
        for library in project_libraries:
            for read_type, fastq_path_key in [('R1', 'fastq_path_r1'), ('R2', 'fastq_path_r2')]:
                target = library[fastq_path_key]
                # Defend against cases where symlink target is None
                # See https://github.com/BCCDC-PHL/auto-fastq-symlink/issues/18
                if target is None:
                    log.warning({"event_type": "fastq_file_not_found", "sequencing_run_id": run_id, "project_id": project_id, "library_id": library['library_id'], "read_type": read_type})
                    continue
                symlinks.append({
                    'target': target,
                    'path': _symlink_filename(target, project['simplify_symlink_filenames']),
                })

        if len(symlinks) > 0:
            symlinks_to_create_by_project_id[project_id] = symlinks

    num_symlinks_to_create_by_project_id = {project_id: len(symlinks) for project_id, symlinks in symlinks_to_create_by_project_id.items()}
    log.debug({
        "event_type": "determine_symlinks_for_run_complete",
        "sequencing_run_id": run_id,
        "total_num_symlinks_to_create": sum(num_symlinks_to_create_by_project_id.values()),
        "num_symlinks_to_create_by_project_id": num_symlinks_to_create_by_project_id,
    })

    return symlinks_to_create_by_project_id


def create_symlinks(config: dict[str, object], symlinks_to_create_by_project_id: dict[str, list[dict[str, str]]], run_id: str) -> dict[str, list[dict[str, str]]]:
    """
    Create a symlinks directory for the run under each project's `fastq_symlinks_dir`.

    The symlinks are created under a temporary hidden directory (`.<run_id>.tmp`), which is then renamed into place.
    `symlinks_complete.json` (which signals to downstream systems that the directory is ready) is only written
    after the rename, and is itself written to a temporary file and renamed, so it never appears in the
    temporary directory and is never partially written. If the application is interrupted before the rename,
    the temporary directory is cleared out and re-built on the next scan.

    :param config: Application config.
    :type config: dict[str, object]
    :param symlinks_to_create_by_project_id: Symlinks to create, indexed by project ID (from `determine_symlinks_to_create_for_run`).
    :type symlinks_to_create_by_project_id: dict[str, list[dict[str, str]]]
    :param run_id: Sequencing run ID.
    :type run_id: str
    :return: Dictionary of symlinks created by project ID.
    :rtype: dict[str, list[dict[str, str]]]
    """
    log.debug({"event_type": "create_symlinks_start", "sequencing_run_id": run_id})
    symlinks_complete_by_project_id = {}
    for project_id, symlinks in symlinks_to_create_by_project_id.items():
        project_fastq_symlinks_dir = config['projects'][project_id]['fastq_symlinks_dir']
        run_symlinks_dir = os.path.join(project_fastq_symlinks_dir, run_id)
        tmp_run_symlinks_dir = os.path.join(project_fastq_symlinks_dir, '.' + run_id + '.tmp')

        if os.path.lexists(tmp_run_symlinks_dir):
            log.warning({"event_type": "removing_incomplete_run_symlinks_dir", "sequencing_run_id": run_id, "project_id": project_id, "path": tmp_run_symlinks_dir})
            shutil.rmtree(tmp_run_symlinks_dir)
        os.makedirs(tmp_run_symlinks_dir)

        symlinks_complete = []
        for symlink in symlinks:
            # The order and naming of the parameters to os.symlink are a bit confusing
            # src = 'target' = the original fastq file
            # dst = 'path' = the symlink
            os.symlink(src=symlink['target'], dst=os.path.join(tmp_run_symlinks_dir, symlink['path']))
            symlinks_complete.append({"target": symlink['target'], "path": os.path.join(run_symlinks_dir, symlink['path'])})

        try:
            os.rename(tmp_run_symlinks_dir, run_symlinks_dir)
        except OSError as e:
            # The run symlinks dir was created by something else since we checked for it. Leave it alone.
            log.error({"event_type": "rename_run_symlinks_dir_failed", "sequencing_run_id": run_id, "project_id": project_id, "run_symlinks_dir": run_symlinks_dir, "error": str(e)})
            shutil.rmtree(tmp_run_symlinks_dir)
            continue

        project_symlinks_complete = {
            'num_symlinks_created': len(symlinks_complete),
            'timestamp': datetime.datetime.now().isoformat(),
        }
        tmp_symlinks_complete_path = os.path.join(run_symlinks_dir, '.symlinks_complete.json.tmp')
        with open(tmp_symlinks_complete_path, 'w') as f:
            f.write(json.dumps(project_symlinks_complete, indent=2) + '\n')
        os.rename(tmp_symlinks_complete_path, os.path.join(run_symlinks_dir, 'symlinks_complete.json'))

        symlinks_complete_by_project_id[project_id] = symlinks_complete

    return symlinks_complete_by_project_id


def scan(config: dict[str, object]) -> Iterable[dict[str, object]]:
    """
    Scan for all sequencing runs. Each run is yielded as it is found, so that it can be symlinked
    before the scan continues.

    :param config: Application config.
    :type config: dict[str, object]
    :return: Sequencing runs that are ready to be symlinked.
    :rtype: Iterable[dict[str, object]]
    """
    log.info({"event_type": "scan_start"})
    num_runs_found = 0
    for run in find_runs(config):
        num_runs_found += 1
        yield run

    log.info({"event_type": "find_runs_complete", "num_runs_found": num_runs_found})


def symlink_run(config: dict[str, object], run: dict[str, object]):
    """
    Determine which symlinks need to be created for one run, based on the current state of the filesystem.
    Then create all symlinks that need to be created.

    :param config: Application config.
    :type config: dict[str, object]
    :param run: Sequencing run info.
    :type config: dict[str, object]
    :return: None
    :rtype: NoneType
    """
    run_id = run['run_id']
    log.debug({"event_type": "symlink_run_start", "sequencing_run_id": run_id})

    symlinks_to_create = determine_symlinks_to_create_for_run(config, run)
    symlinks_complete_by_project_id = create_symlinks(config, symlinks_to_create, run_id)
    log.debug({"event_type": "create_symlinks_complete", "sequencing_run_id": run_id})
    num_symlinks_created_by_project_id = {project_id: len(symlinks) for project_id, symlinks in symlinks_complete_by_project_id.items()}
    log.info({
        "event_type": "symlink_run_complete",
        "sequencing_run_id": run_id,
        "total_num_symlinks_created": sum(num_symlinks_created_by_project_id.values()),
        "num_symlinks_created_by_project_id": num_symlinks_created_by_project_id,
    })
