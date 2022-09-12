import datetime
import json
import logging
import os
import re

import auto_fastq_symlink.samplesheet as ss
import auto_fastq_symlink.db as db


def collect_project_info(config):
    projects = {}
    projects = config['projects']
    # projects_json_safe = auto_fastq_symlink.config.make_config_json_serializable(config)
    # print(json.dumps(projects_json_safe, indent=2))
    # exit(0)
    
    return projects


def _find_fastq_directory(run_dir_path, instrument_type):
    """
    Old (v1) MiSeq directory strucuture: Data/Intensities/BaseCalls
    New (v2) MiSeq directory structure: like: Alignment_1/20220619_120702/Fastq
    NextSeq directory structure: like: Analysis/1/Data/fastq
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


def sanitize_library_id(library_id: str):
    """
    """
    sanitized_library_id = library_id.strip()
    sanitized_library_id = sanitized_library_id.replace('_', '-') # There shouldn't be underscores in library IDs, but if there are, they get swapped with '-' in fastq filenames.
    sanitized_library_id = sanitized_library_id.replace(' ', '-') # There shouldn't be spaces in library IDs, but if there are, they get swapped with '-' in fastq filenames.
    symbols = "\\`*{}[]()>#!$^"
    for s in symbols:
        sanitized_library_id = sanitized_library_id.replace(s, '-')
    sanitized_library_id = re.sub('-+', '-', sanitized_library_id) #
    return sanitized_library_id


def find_libraries(run, samplesheet, fastq_extensions):
    """
    """
    libraries = []
    # If we can't confirm that the fastq directory exists, no sense continuing.
    # Short-circuit and return an empty list.
    if not run['fastq_directory'] or not os.path.exists(run['fastq_directory']):
        logging.error(json.dumps({"event_type": "find_libraries_failed", "run_fastq_directory": run['fastq_directory']}))
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
            continue
        library = {}
        library_id = sanitize_library_id(item[library_id_header])
        if library_id == "" or library_id in found_library_ids:
            continue
        logging.debug(json.dumps({"event_type": "found_library", "library_id": library_id}))
        library['library_id'] = library_id
        library['project_id'] = item[project_header]
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
            logging.debug(json.dumps({"event_type": "found_library_fastq_file", "library_id": library_id, "fastq_path": r1_fastq_path}))
            library['fastq_path_r1'] = r1_fastq_path
        else:
            library['fastq_path_r1'] = None
        if r2_fastq_path and os.path.exists(r2_fastq_path):
            logging.debug(json.dumps({"event_type": "found_library_fastq_file", "library_id": library_id, "fastq_path": r2_fastq_path}))
            library['fastq_path_r2'] = r2_fastq_path
        else:
            library['fastq_path_r2'] = None
        found_library_ids.add(library_id)
        libraries.append(library)

    return libraries
    

def find_runs(run_parent_dirs: list[str], fastq_extensions: list[str]):
    """
    """
    run = {}
    miseq_run_id_regex = "\d{6}_M\d{5}_\d+_\d{9}-[A-Z0-9]{5}"
    nextseq_run_id_regex = "\d{6}_VH\d{5}_\d+_[A-Z0-9]{9}"
    for run_parent_dir in run_parent_dirs:
        subdirs = os.scandir(run_parent_dir)
        for subdir in subdirs:
            run = {}
            instrument_type = None
            run_id = subdir.name
            if re.match(miseq_run_id_regex, run_id):
                instrument_type = "miseq"
            elif re.match(nextseq_run_id_regex, run_id):
                instrument_type = "nextseq"
            if subdir.is_dir() and instrument_type != None:
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
                    samplesheet_to_parse = ss.choose_samplesheet_to_parse(run['samplesheet_files'], run['instrument_type'])
                    if samplesheet_to_parse:
                        logging.debug(json.dumps({"event_type": "samplesheet_found", "sequencing_run_id": run_id, "samplesheet_path": samplesheet_to_parse}))
                    else:
                        logging.error(json.dumps({"event_type": "samplesheet_not_found", "sequencing_run_id": run_id, "samplesheet_path": samplesheet_to_parse}))
                        run['libraries'] = []
                        continue

                    run['parsed_samplesheet'] = samplesheet_to_parse
                    samplesheet = ss.parse_samplesheet(samplesheet_to_parse, run['instrument_type'])
                    libraries = find_libraries(run, samplesheet, fastq_extensions)
                    run['libraries'] = libraries

                    yield run


def find_symlinks(projects):
    """
    """
    symlinks_by_project = {}
    for project_id, project in projects.items():
        symlinks_by_project[project_id] = []
        fastq_symlinks_dir = project['fastq_symlinks_dir']
        if os.path.exists(fastq_symlinks_dir):
            project_symlinks_by_run_dirs = os.scandir(fastq_symlinks_dir)
            for project_symlinks_by_run_dir in project_symlinks_by_run_dirs:
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


def determine_symlinks_to_create(config):
    """
    """
    symlinks_to_create_by_project_id = {}

    existing_symlinks = db.get_symlinks(config)
    existing_project_target_pairs = set()
    for symlink in existing_symlinks:
        project_target_pair = (symlink['project_id'], symlink['target'])
        existing_project_target_pairs.add(project_target_pair)

    for project_id in config['projects']:
        symlinks_to_create_by_project_id[project_id] = []
        project_libraries = db.get_libraries_by_project_id(config, project_id)
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

    return symlinks_to_create_by_project_id


def create_symlinks(config: dict[str, object], symlinks_to_create_by_project_id: dict[str, object]):
    """
    """
    symlinks_complete = {'num_symlinks_created': 0}
    for project_id, symlinks in symlinks_to_create_by_project_id.items():
        project_fastq_symlinks_dir = config['projects'][project_id]['fastq_symlinks_dir']

        symlinks_complete = {'num_symlinks_created': 0}
        for symlink in symlinks:
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
            except FileExistsError as e:
                logging.warning(json.dumps({
                    "event_type": "attempted_to_create_existing_symlink",
                    "symlink_target": symlink['target'],
                    "symlink_path": symlink['path'],
                }))

            symlinks_complete['num_symlinks_created'] += 1

            timestamp = datetime.datetime.now().isoformat()
            symlinks_complete['timestamp'] = timestamp

            with open(os.path.join(symlink_parent_dir, 'symlinks_complete.json'), 'w') as f:
                f.write(json.dumps(symlinks_complete, indent=2) + '\n')
                
    return symlinks_complete


def scan(config: dict[str, object]):
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

    logging.debug(json.dumps({"event_type": "find_and_store_runs_start"}))
    num_runs_found = 0
    for run in find_runs(config['run_parent_dirs'], config['fastq_extensions']):
        db.store_run(config, run)
        num_runs_found += 1
    logging.debug(json.dumps({"event_type": "find_and_store_runs_complete", "num_runs_found": num_runs_found}))

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

    logging.info(json.dumps({"event_type": "scan_complete", "num_runs_found": num_runs_found, "num_symlinks_found": num_symlinks_found}))


def symlink(config: dict[str, object]):
    """
    Determine which symlinks need to be created, based on the current state of the database.
    Then create all symlinks that need to be created.

    :param config: Application config.
    :type config: dict[str, object]
    :return: None
    :rtype: NoneType
    """
    logging.info(json.dumps({"event_type": "symlink_start"}))
    logging.debug(json.dumps({"event_type": "determine_symlinks_start"}))
    symlinks_to_create = determine_symlinks_to_create(config)
    num_symlinks_to_create = 0
    for project_id, symlinks in symlinks_to_create.items():
        num_symlinks_to_create += len(symlinks)
    logging.debug(json.dumps({"event_type": "determine_symlinks_complete", "num_symlinks_to_create": num_symlinks_to_create}))

    logging.debug(json.dumps({"event_type": "create_symlinks_start"}))
    symlinks_complete = create_symlinks(config, symlinks_to_create)
    num_symlinks_created = symlinks_complete['num_symlinks_created']
    logging.debug(json.dumps({"event_type": "create_symlinks_complete"}))
    logging.info(json.dumps({"event_type": "symlink_complete", "num_symlinks_created": num_symlinks_created}))

