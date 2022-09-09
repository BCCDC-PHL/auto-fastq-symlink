import datetime
import json
import logging
import os
import re
import time

import auto_fastq_symlink.config
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
            fastq_directory = os.path.join(run_dir_path, "Data", "Intensities", "BaseCalls")

    elif instrument_type == 'nextseq':
        analysis_subdirs = list([subdir.name for subdir in os.scandir(os.path.join(run_dir_path, "Analysis"))])
        greatest_analysis_subdir_num = sorted(analysis_subdirs)[0]
        fastq_directory = os.path.join(run_dir_path, "Analysis", str(greatest_analysis_subdir_num), "Data", "fastq")

    return fastq_directory    


def _find_libraries(run, samplesheet, fastq_extensions):
    """
    """
    libraries = []
    libraries_section = None
    if run['instrument_type'] == 'miseq':
        libraries_section = 'data'
    elif run['instrument_type'] == 'nextseq':
        libraries_section = 'cloud_data'

    project_header = None
    if run['instrument_type'] == 'miseq':
        project_header = 'sample_project'
    elif run['instrument_type'] == 'nextseq':
        project_header = 'project_name'
        
    has_correct_extension = lambda x: any([x.endswith(ext) for ext in fastq_extensions])
    run_fastq_dir_contents = os.listdir(run['fastq_directory'])
    run_fastq_files = set(filter(lambda x: all([has_correct_extension(x), os.path.isfile(os.path.join(run['fastq_directory'], x))]), run_fastq_dir_contents))

    for item in samplesheet[libraries_section]:
        library = {}
        library['library_id'] = item['sample_id']
        library['project_id'] = item[project_header]
        r1_fastq_filename = list(filter(lambda x: re.match(item['sample_id'] + '.*' + '_R1_' + '.*', x), run_fastq_files))[0]
        r2_fastq_filename = list(filter(lambda x: re.match(item['sample_id'] + '.*' + '_R2_' + '.*', x), run_fastq_files))[0]
        r1_fastq_path = os.path.join(run['fastq_directory'], r1_fastq_filename)
        r2_fastq_path = os.path.join(run['fastq_directory'], r2_fastq_filename)
        library['fastq_path_r1'] = r1_fastq_path
        library['fastq_path_r2'] = r2_fastq_path
        libraries.append(library)

    return libraries
    

def find_runs(run_parent_dirs, fastq_extensions):
    """
    """
    runs = {}
    miseq_run_id_regex = "\d{6}_M\d{5}_\d+_\d{9}-[A-Z0-9]{5}"
    nextseq_run_id_regex = "\d{6}_VH\d{5}_\d+_[A-Z0-9]{9}"
    for run_parent_dir in run_parent_dirs:
        subdirs = os.scandir(run_parent_dir)
        for subdir in subdirs:
            instrument_type = None
            if re.match(miseq_run_id_regex, subdir.name):
                instrument_type = "miseq"
            elif re.match(nextseq_run_id_regex, subdir.name):
                instrument_type = "nextseq"
            if subdir.is_dir() and instrument_type != None:
                samplesheet_paths = ss.find_samplesheets(subdir.path, instrument_type)
                fastq_directory = _find_fastq_directory(subdir.path, instrument_type)
                
                runs[subdir.name] = {
                    "run_id": subdir.name,
                    "instrument_type": instrument_type,
                    "samplesheet_files": samplesheet_paths,
                    "run_directory": subdir.path,
                    "fastq_directory": fastq_directory,
                }

    for run_id, run in runs.items():
        samplesheet_to_parse = ss.choose_samplesheet_to_parse(run['samplesheet_files'], run['instrument_type'])
        run['parsed_samplesheet'] = samplesheet_to_parse

        samplesheet = ss.parse_samplesheet(samplesheet_to_parse, run['instrument_type'])
        libraries = _find_libraries(run, samplesheet, fastq_extensions)
        run['libraries'] = libraries

    return runs


def find_symlinks(projects):
    """
    """
    symlinks_by_project = {}
    for project_id, project in projects.items():
        symlinks_by_project[project_id] = []
        fastq_symlinks_dir = project['fastq_symlinks_dir']
        if os.path.exists(fastq_symlinks_dir):
            symlinks = []
            project_symlinks_by_run_dirs = os.scandir(fastq_symlinks_dir)
            for project_symlinks_by_run_dir in project_symlinks_by_run_dirs:
                project_symlinks_by_run_dir_contents = os.scandir(project_symlinks_by_run_dir)
                for dir_item in project_symlinks_by_run_dir_contents:
                    if os.path.islink(dir_item):
                        path = dir_item.path
                        target = os.path.realpath(dir_item.path)
                        symlink = {
                            "path": path,
                            "target": target,
                        }
                        symlinks_by_project[project_id].append(symlink)
            
    return symlinks_by_project


def determine_symlinks_to_create(config):
    """
    """
    symlinks_to_create_by_project_id = {}

    libraries_by_project_id = {}

    existing_symlinks = db.get_symlinks(config)
    existing_project_target_pairs = set()
    for symlink in existing_symlinks:
        project_target_pair = (symlink['project_id'], symlink['target'])
        existing_project_target_pairs.add(project_target_pair)

    for project_id in config['projects']:
        symlinks_to_create_by_project_id[project_id] = []
        project_libraries = db.get_libraries_by_project_id(config, project_id)
        for library in project_libraries:
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


def create_symlinks(config, symlinks_to_create_by_project_id):
    """
    """
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


def scan(config):
    """
    Scanning involves looking for all existing runs and storing them to the database,
    then looking for all existing symlinks and storing them to the database.
    At the end of a scan, we should be able to determine which (if any) symlinks need to be created.
    """
    logging.info(json.dumps({"event_type": "scan_start"}))
    logging.debug(json.dumps({"event_type": "collect_projects_start"}))
    projects = collect_project_info(config)
    num_projects = len(projects)
    logging.debug(json.dumps({"event_type": "collect_projects_complete", "num_projects": num_projects}))

    logging.debug(json.dumps({"event_type": "store_projects_start"}))
    db.store_projects(config, projects)
    logging.debug(json.dumps({"event_type": "store_projects_complete"}))

    logging.debug(json.dumps({"event_type": "find_runs_start"}))
    runs = find_runs(config['run_parent_dirs'], config['fastq_extensions'])
    num_runs = len(runs)
    logging.debug(json.dumps({"event_type": "find_runs_complete", "num_runs_found": num_runs}))

    logging.debug(json.dumps({"event_type": "store_runs_start"}))
    db.store_runs(config, runs)
    logging.debug(json.dumps({"event_type": "store_runs_complete"}))

    logging.debug(json.dumps({"event_type": "find_symlinks_start"}))
    symlinks_by_destination_dir = find_symlinks(config['projects'])
    num_symlinks = len(symlinks_by_destination_dir)
    logging.debug(json.dumps({"event_type": "find_symlinks_complete", "num_symlinks_found": num_symlinks}))
    logging.debug(json.dumps({"event_type": "store_symlinks_start"}))
    db.store_symlinks(config, symlinks_by_destination_dir)
    logging.debug(json.dumps({"event_type": "store_symlinks_complete"}))
    logging.debug(json.dumps({"event_type": "delete_nonexistent_symlinks_start"}))
    db.delete_nonexistent_symlinks(config)
    logging.debug(json.dumps({"event_type": "delete_nonexistent_symlinks_complete"}))
    logging.info(json.dumps({"event_type": "scan_complete", "num_runs_found": num_runs, "num_symlinks_found": num_symlinks}))


def symlink(config):
    """
    Determine which symlinks need to be created, based on the current state of the database.
    Then create all symlinks that need to be created.
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

