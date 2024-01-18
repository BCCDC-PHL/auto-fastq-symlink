import datetime
import json
import logging
import os
import re
from typing import Iterable, Optional

import jsonschema

import auto_fastq_symlink.samplesheet as ss
import auto_fastq_symlink.db as db


def collect_project_info(config: dict[str, object]) -> dict[str, str]:
    """
    Collect project info from config and return it as a dictionary.

    :param config: Application config.
    :type config: dict[str, object]
    :return: Project info, indexed by project ID.
    :rtype: dict[str, str]
    """
    projects = {}
    projects = config['projects']
    
    return projects


def _find_fastq_directory(run_dir_path, instrument_type):
    """
    Old (v1) MiSeq directory strucuture: Data/Intensities/BaseCalls
    New (v2) MiSeq directory structure: like: Alignment_1/20220619_120702/Fastq
    NextSeq directory structure: like: Analysis/1/Data/fastq

    :param run_dir_path: Path to the sequencing run directory.
    :type run_dir_path: str
    :param instrument_type: Instrument type ('miseq' or 'nextseq')
    :type instrument_type: str
    :return: Path to the fastq directory, or None if it can't be found.
    """
    fastq_directory = ""
    if instrument_type == 'miseq':
        run_subdirs = set([subdir.name for subdir in os.scandir(run_dir_path)])
        version_2_miseq_dir_structure = any([re.match("Alignment_\d+", subdir) for subdir in run_subdirs])
        if version_2_miseq_dir_structure:
            alignment_subdirs = list(filter(lambda x: re.match("Alignment_\d+", x), run_subdirs))
            alignment_subdir_nums = list(map(lambda x: int(x.split('_')[1]), alignment_subdirs))
            greatest_alignment_subdir_num = sorted(alignment_subdir_nums, reverse=True)[0]
            alignment_directory = os.path.join(run_dir_path, "Alignment_" + str(greatest_alignment_subdir_num))
            alignment_subdirs = sorted([subdir.name for subdir in os.scandir(alignment_directory)], reverse=True)
            greatest_alignment_subdir = alignment_subdirs[0]
            fastq_directory = os.path.join(alignment_directory, greatest_alignment_subdir, "Fastq")
        else:
            data_intensities_basecalls_dir_path = os.path.join(run_dir_path, "Data", "Intensities", "BaseCalls")
            if os.path.exists(data_intensities_basecalls_dir_path):
                fastq_directory = data_intensities_basecalls_dir_path
            else:
                fastq_directory = None

    elif instrument_type == 'nextseq':
        analysis_dir_path = os.path.join(run_dir_path, "Analysis")
        if os.path.exists(analysis_dir_path):
            analysis_subdirs = list([subdir.name for subdir in os.scandir(analysis_dir_path)])
            greatest_analysis_subdir_num = sorted(analysis_subdirs, reverse=True)[0]
            fastq_directory = os.path.join(run_dir_path, "Analysis", str(greatest_analysis_subdir_num), "Data", "fastq")
        else:
            fastq_directory = None

    return fastq_directory


def _determine_libraries_section(samplesheet, instrument_type):
    """
    Determine which section of the SampleSheet contains the library info.

    :param samplesheet: The parsed SampleSheet.
    :type samplesheet: dict[str, object]
    :param instrument_type: Instrument type ('miseq' or 'nextseq')
    :type instrument_type: str
    :return: The section of the SampleSheet that contains the library info.
    :rtype: str | None
    """
    if instrument_type == 'miseq':
        libraries_section = 'data'
    elif instrument_type == 'nextseq':
        libraries_section = 'cloud_data'
    else:
        libraries_section = None

    return libraries_section


def _determine_project_header(samplesheet, instrument_type):
    """
    For MiSeq runs, the SampleSheet has two fields that can either be used for the project ID.

    :param samplesheet: The parsed SampleSheet.
    :type samplesheet: dict[str, object]
    :param instrument_type: Instrument type ('miseq' or 'nextseq')
    :type instrument_type: str
    :return: The header of the field that contains the project ID.
    :rtype: str | None
    """
    if instrument_type == 'miseq':
        project_header = 'sample_project'
    elif instrument_type == 'nextseq':
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
    :param instrument_type: Instrument type ('miseq' or 'nextseq')
    :type instrument_type: str
    :return: The header of the field that contains the library ID.
    :rtype: str | None
    """
    libraries_section = _determine_libraries_section(samplesheet, instrument_type)
    if instrument_type == 'miseq':
        sample_ids = [library['sample_id'] for library in samplesheet[libraries_section]]
        all_sample_ids_blank = all([sample_id == "" for sample_id in sample_ids])
        all_sample_ids_s_plus_digits = all([re.match("S\d+", sample_id) for sample_id in sample_ids])
        all_sample_ids_only_digits = all([re.match("\d+", sample_id) for sample_id in sample_ids])
        if all_sample_ids_blank or all_sample_ids_s_plus_digits or all_sample_ids_only_digits:
            library_id_header = 'sample_name'
        else:
            library_id_header = 'sample_id'
    elif instrument_type == 'nextseq':
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
    if not run['fastq_directory'] or not os.path.exists(run['fastq_directory']):
        logging.error(json.dumps({"event_type": "find_libraries_failed", "sequencing_run_id": run_id, "run_fastq_directory": run['fastq_directory']}))
        return libraries
    libraries_section = _determine_libraries_section(samplesheet, run['instrument_type'])
    project_header = _determine_project_header(samplesheet, run['instrument_type'])

    has_correct_extension = lambda x: any([x.endswith(ext) for ext in fastq_extensions])
    run_fastq_dir_contents = os.listdir(run['fastq_directory'])
    run_fastq_files = set(filter(lambda x: all([has_correct_extension(x), os.path.isfile(os.path.join(run['fastq_directory'], x))]), run_fastq_dir_contents))

    library_id_header = _determine_library_id_header(samplesheet, run['instrument_type'])
    
    found_library_ids = set()
    for item in samplesheet[libraries_section]:
        required_item_keys = [project_header, library_id_header]
        if not all([k in item for k in required_item_keys]):
            logging.warning(json.dumps({"event_type": "library_missing_required_fields", "required_fields": required_item_keys, "library": item}))
            continue
        library = {}
        library_id = _sanitize_library_id(item[library_id_header])
        project_id = item[project_header]
        if library_id == "" or library_id in found_library_ids:
            continue
        library['library_id'] = library_id
        library['project_id'] = project_id
        logging.debug(json.dumps({"event_type": "found_library", "library_id": library_id, "samplesheet_project_id": project_id}))
        r1_fastq_filenames = list(filter(lambda x: re.match(library_id + '.*' + '_R1_' + '.*', x), run_fastq_files))
        if len(r1_fastq_filenames) > 0:
            r1_fastq_filename = r1_fastq_filenames[0]
        else:
            r1_fastq_filename = None
        r2_fastq_filenames = list(filter(lambda x: re.match(library_id + '.*' + '_R2_' + '.*', x), run_fastq_files))
        if len(r2_fastq_filenames) > 0:
            r2_fastq_filename = r2_fastq_filenames[0]
        else:
            r2_fastq_filename = None
        r1_fastq_path = os.path.join(run['fastq_directory'], r1_fastq_filename) if r1_fastq_filename else None
        r2_fastq_path = os.path.join(run['fastq_directory'], r2_fastq_filename) if r2_fastq_filename else None
        if r1_fastq_path and os.path.exists(r1_fastq_path):
            logging.debug(json.dumps({"event_type": "found_library_fastq_file", "library_id": library_id, "read_type": "R1", "fastq_path": r1_fastq_path}))
            library['fastq_path_r1'] = r1_fastq_path
        else:
            library['fastq_path_r1'] = None
        if r2_fastq_path and os.path.exists(r2_fastq_path):
            logging.debug(json.dumps({"event_type": "found_library_fastq_file", "library_id": library_id, "read_type": "R2", "fastq_path": r2_fastq_path}))
            library['fastq_path_r2'] = r2_fastq_path
        else:
            library['fastq_path_r2'] = None
        found_library_ids.add(library_id)
        libraries.append(library)
        

    return libraries
    

def find_runs(config: dict[str, object]) -> Iterable[Optional[dict[str, object]]]:
    """
    Find all sequencing runs under all of the `run_parent_dirs` from the config.
    Runs are found by matching sub-directory names against the following regexes: `"\d{6}_M\d{5}_\d+_\d{9}-[A-Z0-9]{5}"` (MiSeq) and `"\d{6}_VH\d{5}_\d+_[A-Z0-9]{9}"` (NextSeq)

    :param config: Application config.
    :type config: dict[str, object]
    :return: Dictionary of sequencin run info, indexed by sequencing run ID.
    :rtype: Iterable[dict[str, object]]
    """
    run = {}
    miseq_run_id_regex = "\d{6}_M\d{5}_\d+_\d{9}-[A-Z0-9]{5}"
    nextseq_run_id_regex = "\d{6}_VH\d{5}_\d+_[A-Z0-9]{9}"
    run_parent_dirs = config['run_parent_dirs']
    fastq_extensions = config['fastq_extensions']
    for run_parent_dir in run_parent_dirs:
        if not os.path.exists(run_parent_dir):
            logging.warning(json.dumps({"event_type": "run_parent_dir_does_not_exist", "run_parent_dir": run_parent_dir}))
            continue
        subdirs = os.scandir(run_parent_dir)
        for subdir in subdirs:
            run = {}
            instrument_type = None
            run_id = subdir.name
            if re.match(miseq_run_id_regex, run_id):
                instrument_type = "miseq"
            elif re.match(nextseq_run_id_regex, run_id):
                instrument_type = "nextseq"

            subdir_is_dir = os.path.isdir(subdir.path)
            upload_complete_file_exists = os.path.exists(os.path.join(subdir.path, "upload_complete.json"))
            qc_check_complete_file_exists = os.path.exists(os.path.join(subdir.path, "qc_check_complete.json"))
            determined_instrument_type = instrument_type != None

            passed_run_qc_check = False
            if qc_check_complete_file_exists:
                try:
                    qc_check = json.load(open(os.path.join(subdir.path, "qc_check_complete.json")))
                except json.decoder.JSONDecodeError as e:
                    logging.error(json.dumps({"event_type": "qc_check_json_decode_error", "sequencing_run_id": run_id, "error": str(e)}))
                    qc_check = {}
                overall_pass_fail = qc_check.get('overall_pass_fail', None)
                if overall_pass_fail is not None and re.match("PASS", overall_pass_fail, re.IGNORECASE):
                    passed_run_qc_check = True

            conditions_checked = {
                'subdir_is_directory': subdir_is_dir,
                'upload_complete': upload_complete_file_exists,
                'qc_check_complete': qc_check_complete_file_exists,
                'passed_run_qc_check': passed_run_qc_check,
                'determined_instrument_type': determined_instrument_type,
            }
            conditions_met = [v for k, v in conditions_checked.items()]
            if not all(conditions_met):
                logging.info(json.dumps({"event_type": "skipped_run", "sequencing_run_id": run_id, "conditions_checked": conditions_checked}))
                yield None

            if all(conditions_met):
                logging.info(json.dumps({"event_type": "scan_run_start", "sequencing_run_id": run_id}))
                samplesheet_paths = ss.find_samplesheets(subdir.path, instrument_type)
                fastq_directory = _find_fastq_directory(subdir.path, instrument_type)
                if fastq_directory != None:
                    logging.debug(json.dumps({"event_type": "sequencing_run_found", "sequencing_run_id": run_id}))
                    run = {
                        "run_id": run_id,
                        "instrument_type": instrument_type,
                        "samplesheet_files": samplesheet_paths,
                        "run_directory": subdir.path,
                        "fastq_directory": fastq_directory,
                    }
                    samplesheet_to_parse = ss.choose_samplesheet_to_parse(run['samplesheet_files'], run['instrument_type'], run_id)
                    if samplesheet_to_parse:
                        logging.debug(json.dumps({"event_type": "samplesheet_found", "sequencing_run_id": run_id, "samplesheet_path": samplesheet_to_parse}))
                    else:
                        logging.error(json.dumps({"event_type": "samplesheet_not_found", "sequencing_run_id": run_id, "samplesheet_path": samplesheet_to_parse}))
                        run['libraries'] = []
                        continue

                    run['parsed_samplesheet'] = samplesheet_to_parse
                    try:
                        samplesheet = ss.parse_samplesheet(samplesheet_to_parse, run['instrument_type'])
                    except jsonschema.ValidationError as e:
                        yield None
                    libraries = find_libraries(run, samplesheet, fastq_extensions)
                    for library in libraries:
                        if library['project_id'] in config['project_id_translation']:
                            samplesheet_project_id = library['project_id']
                            symlinking_project_id = config['project_id_translation'][samplesheet_project_id]
                            library['project_id'] = symlinking_project_id
                        elif library['project_id'] == '':
                            library['project_id'] = None
                    run['libraries'] = libraries

                    yield run


def find_symlinks(projects):
    """
    Find all existing symlinks under each project's `fastq_symlinks_dir`

    :param projects:
    :type projects: dict[str, object]
    :return: Dict of existing symlinks, indexed by project ID
    :rtype: dict[str, object]
    """
    symlinks_by_project = {}
    for project_id, project in projects.items():
        symlinks_by_project[project_id] = []
        fastq_symlinks_dir = project['fastq_symlinks_dir']
        if os.path.exists(fastq_symlinks_dir):
            project_symlinks_by_run_dirs = os.scandir(fastq_symlinks_dir)
            for project_symlinks_by_run_dir in project_symlinks_by_run_dirs:
                if not os.path.isdir(project_symlinks_by_run_dir):
                    continue
                project_symlinks_by_run_dir_contents = os.scandir(project_symlinks_by_run_dir)
                for dir_item in project_symlinks_by_run_dir_contents:
                    if os.path.islink(dir_item):
                        path = dir_item.path
                        target = os.path.realpath(dir_item.path)
                        symlink = {
                            "sequencing_run_id": os.path.basename(project_symlinks_by_run_dir),
                            "path": path,
                            "target": target,
                        }
                        symlinks_by_project[project_id].append(symlink)

    return symlinks_by_project


def determine_symlinks_to_create_for_run(config: dict[str, object], run_id: str) -> dict[str, dict[str, str]]:
    """
    :param config: Application config
    :type config: dict[str, object]
    :param run_id: Sequencing run identifier
    :type run_id: str
    :return: Dictionary of file paths to fastq files for which symlinks should be created, indexed by project ID
    :rtype: dict[str, dict[str, str]]
    """
    logging.debug(json.dumps({"event_type": "determine_symlinks_for_run_start", "sequencing_run_id": run_id}))
    symlinks_to_create_by_project_id = {}

    existing_symlinks = db.get_symlinks_by_run_id(config, run_id)
    existing_project_target_pairs = set()
    for symlink in existing_symlinks:
        project_target_pair = (symlink['project_id'], symlink['target'])
        existing_project_target_pairs.add(project_target_pair)

    for project_id in config['projects']:
        symlinks_to_create_by_project_id[project_id] = []
        project_libraries = db.get_libraries_by_project_id_and_run_id(config, project_id, run_id)
        project_excluded_runs = config['projects'][project_id]['excluded_runs']
        project_excluded_libraries = config['projects'][project_id]['excluded_libraries']

        for library in project_libraries:
            if (library['sequencing_run_id'] not in project_excluded_runs) and (library['library_id'] not in project_excluded_libraries):
                project_fastq_path_r1_pair = (library['project_id'], library['fastq_path_r1'])
                if project_fastq_path_r1_pair not in existing_project_target_pairs:
                    fastq_path_r1 = {
                        'project_id': library['project_id'],
                        'sequencing_run_id': library['sequencing_run_id'],
                        'target': library['fastq_path_r1'],
                    }
                    symlinks_to_create_by_project_id[project_id].append(fastq_path_r1)
                    fastq_path_r2 = {
                        'project_id': library['project_id'],
                        'sequencing_run_id': library['sequencing_run_id'],
                        'target': library['fastq_path_r2'],
                    }
                    symlinks_to_create_by_project_id[project_id].append(fastq_path_r2)

    total_num_symlinks_to_create = 0
    num_symlinks_to_create_by_project_id = {}
    for project_id, symlinks in symlinks_to_create_by_project_id.items():
        num_symlinks_to_create_by_project_id[project_id] = len(symlinks)
    logging.debug(json.dumps({
        "event_type": "determine_symlinks_for_run_complete",
        "sequencing_run_id": run_id,
        "total_num_symlinks_to_create": total_num_symlinks_to_create,
        "num_symlinks_to_create_by_project_id": num_symlinks_to_create_by_project_id,
    }))


    return symlinks_to_create_by_project_id


def create_symlinks(config: dict[str, object], symlinks_to_create_by_project_id: dict[str, list[dict[str, str]]], run_id: str):
    """
    :param config:
    :type config: dict[str, object]
    :param symlinks_to_create_by_project_id:
    :type symlinks_to_create_by_project_id: dict[str, object]
    :param run_id:
    :type run_id: str
    :return: Dictionary of symlinks created by project ID.
    :rtype: dict[str, list[dict[str, str]]]
    """
    logging.debug(json.dumps({"event_type": "create_symlinks_start", "sequencing_run_id": run_id}))
    symlinks_complete_by_project_id = {}
    for project_id, symlinks in symlinks_to_create_by_project_id.items():
        project_fastq_symlinks_dir = config['projects'][project_id]['fastq_symlinks_dir']

        symlinks_complete_by_project_id[project_id] = []
        for symlink in symlinks:
            # Defend against cases where symlink target is None
            # See https://github.com/BCCDC-PHL/auto-fastq-symlink/issues/18
            if symlink['target'] is None:
                continue

            symlink_parent_dir = os.path.join(project_fastq_symlinks_dir, symlink['sequencing_run_id'])

            if not os.path.exists(symlink_parent_dir):
                os.makedirs(symlink_parent_dir)

            if config['projects'][project_id]['simplify_symlink_filenames']:
                target_basename = os.path.basename(symlink['target'])
                r1_r2_match = re.search("_(R[12])_", target_basename)
                if r1_r2_match:
                    r1_r2 = r1_r2_match.group(1)
                else:
                    r1_r2 = ""
                symlink_filename = os.path.basename(symlink['target']).split('_')[0] + '_' + r1_r2 + '.fastq.gz'
            else:
                symlink_filename = os.path.basename(symlink['target'])

            symlink['path'] = os.path.join(symlink_parent_dir, symlink_filename)

            # The order and naming of the parameters to os.symlink are a bit confusing
            # src = 'target' = the original fastq file
            # dst = 'path' = the symlink
            try:
                os.symlink(src=symlink['target'], dst=symlink['path'])
                symlinks_complete_by_project_id[project_id].append({"target": symlink['target'], "path": symlink['path']})
            except FileExistsError as e:
                logging.warning(json.dumps({
                    "event_type": "attempted_to_create_existing_symlink",
                    "sequencing_run_id": run_id,
                    "symlink_target": symlink['target'],
                    "symlink_path": symlink['path'],
                }))

            project_symlinks_complete = {}
            project_symlinks_complete['num_symlinks_created'] = len(symlinks_complete_by_project_id[project_id])
            timestamp = datetime.datetime.now().isoformat()
            project_symlinks_complete['timestamp'] = timestamp

            with open(os.path.join(symlink_parent_dir, 'symlinks_complete.json'), 'w') as f:
                f.write(json.dumps(project_symlinks_complete, indent=2) + '\n')
                
    return symlinks_complete_by_project_id


def scan(config: dict[str, object]) -> Iterable[Optional[dict[str, object]]]:
    """
    Scanning involves looking for all existing runs and storing them to the database,
    then looking for all existing symlinks and storing them to the database.
    At the end of a scan, we should be able to determine which (if any) symlinks need to be created.

    :param config: Application config.
    :type config: dict[str, object]
    :return: None
    :rtype: NoneType
    """
    logging.info(json.dumps({"event_type": "scan_start"}))
    logging.debug(json.dumps({"event_type": "collect_projects_start"}))
    projects = collect_project_info(config)
    num_projects = len(projects)
    logging.debug(json.dumps({"event_type": "collect_projects_complete", "num_projects": num_projects}))

    logging.debug(json.dumps({"event_type": "store_projects_start"}))
    db.store_projects(config, projects)
    logging.debug(json.dumps({"event_type": "store_projects_complete"}))

    logging.debug(json.dumps({"event_type": "find_symlinks_start"}))
    num_symlinks_found = 0
    symlinks_by_destination_dir = find_symlinks(config['projects'])
    for destination_dir, symlinks in symlinks_by_destination_dir.items():
        num_symlinks_found += len(symlinks)
    logging.debug(json.dumps({"event_type": "find_symlinks_complete", "num_symlinks_found": num_symlinks_found}))

    logging.debug(json.dumps({"event_type": "store_symlinks_start"}))
    db.store_symlinks(config, symlinks_by_destination_dir)
    logging.debug(json.dumps({"event_type": "store_symlinks_complete"}))

    logging.debug(json.dumps({"event_type": "delete_nonexistent_symlinks_start"}))
    db.delete_nonexistent_symlinks(config)
    logging.debug(json.dumps({"event_type": "delete_nonexistent_symlinks_complete"}))

    logging.debug(json.dumps({"event_type": "find_and_store_runs_start"}))
    num_runs_found = 0
    for run in find_runs(config):
        if run is not None:
            db.store_run(config, run)
            num_runs_found += 1
        yield run

    logging.info(json.dumps({"event_type": "find_and_store_runs_complete", "num_runs_found": num_runs_found}))


def symlink_run(config: dict[str, object], run: dict[str, object]):
    """
    Determine which symlinks need to be created for one run, based on the current state of the database.
    Then create all symlinks that need to be created.

    :param config: Application config.
    :type config: dict[str, object]
    :param run: Sequencing run info.
    :type config: dict[str, object]
    :return: None
    :rtype: NoneType
    """
    run_id = run['run_id']
    logging.debug(json.dumps({"event_type": "symlink_run_start", "sequencing_run_id": run_id}))

    symlinks_to_create = determine_symlinks_to_create_for_run(config, run_id)
    symlinks_complete_by_project_id = create_symlinks(config, symlinks_to_create, run_id)
    total_num_symlinks_created = 0
    for project_id, symlinks_complete in symlinks_complete_by_project_id.items():
        total_num_symlinks_created += len(symlinks_complete)
    logging.debug(json.dumps({"event_type": "create_symlinks_complete", "sequencing_run_id": run_id}))
    num_symlinks_created_by_project_id = {project_id: len(symlinks) for project_id, symlinks in symlinks_complete_by_project_id.items()}
    logging.info(json.dumps({
        "event_type": "symlink_run_complete",
        "sequencing_run_id": run_id,
        "total_num_symlinks_created": total_num_symlinks_created,
        "num_symlinks_created_by_project_id": num_symlinks_created_by_project_id,
    }))
