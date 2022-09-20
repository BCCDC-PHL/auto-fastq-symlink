import csv
import json
import logging

def _parse_list(list_path: str) -> list[str]:
    """
    Parse a 'list' file. List files are simple single-column text files, with one entry per line.

    :param list_path: Path to the list file.
    :type list_path: str
    :return: Contents of list file.
    :rtype: list[str]
    """
    items = []
    with open(list_path, 'r') as f:
        for line in f:
            item = line.strip()
            items.append(item)

    return items


def parse_projects_definition_file(projects_definition_file_path: str) -> dict[str, object]:
    """
    Parse a 'projects definition file'. These files are .csv format, and should include the following fields:

    `project_id`
    `fastq_symlinks_dir`
    `excluded_runs_list`
    `excluded_libraries_list`
    `simplify_symlink_filenames`

    :param projects_definition_file_path: Path to the projects definition file.
    :type projects_definition_file_path: str
    :return: Dictionary containing project-related info, indexed by project ID.
    :rtype: dict[str, str]
    """
    projects = {}
    with open(projects_definition_file_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['simplify_symlink_filenames'].lower() in ['true', '1', 't']:
                row['simplify_symlink_filenames'] = True
            else:
                row['simplify_symlink_filenames'] = False
            projects[row['project_id']] = row

    return projects


def parse_project_id_translation_file(project_id_translation_file_path, projects):
    """
    :param project_id_translation_file_path:
    :type project_id_translation_file_path: str
    :param projects:
    :type projects: dict[str, object]
    :return: 
    :rtype: dict[str, str]
    """

    project_translation = {}
    with open(project_id_translation_file_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            samplesheet_project_id = row['samplesheet_project_id']
            symlinking_project_id = row['symlinking_project_id']
            if symlinking_project_id in projects:
                project_translation[samplesheet_project_id] = symlinking_project_id

    return project_translation


def load_config(config_path: str) -> dict[str, object]:
    """
    :param config_path: Path to application config file (json format)
    :type config_path: str
    :return: Application config dictionary.
    :rtype: dict[str, object]
    """
    config = None
    with open(config_path, 'r') as f:
        config = json.load(f)

    if 'projects_definition_file' in config:
        projects = parse_projects_definition_file(config['projects_definition_file'])
        config['projects'] = projects
        if 'project_id_translation_file' in config:
            project_id_translation = parse_project_id_translation_file(config['project_id_translation_file'], config['projects'])
            config['project_id_translation'] = project_id_translation
        else:
            config['project_id_translation'] = {}
            
        for project_id, project in projects.items():
            if 'excluded_runs_list' in project:
                try:
                    excluded_runs_list = _parse_list(project['excluded_runs_list'])
                except Exception as e:
                    # If we can't load the excluded runs, we shouldn't proceed.
                    # If we did, we might create symlinks that shouldn't be created,
                    # which could have further effects downstream
                    logging.error("Error loading excluded runs list: " + project['excluded_runs_list'])
                    exit(-1)
                excluded_runs = set(excluded_runs_list)
                config['projects'][project_id]['excluded_runs'] = excluded_runs
            else:
                config['projects'][project_id]['excluded_runs'] = set()

            if 'excluded_libraries_list' in project:
                try:
                    excluded_libraries_list = _parse_list(project['excluded_libraries_list'])
                except Exception as e:
                    # If we can't load the excluded libraries, we shouldn't proceed.
                    # If we did, we might create symlinks that shouldn't be created,
                    # which could have further effects downstream
                    logging.error("Error loading excluded libraries list: " + project['excluded_libraries_list'])
                    exit(-1)
                excluded_libraries = set(excluded_libraries_list)
                config['projects'][project_id]['excluded_libraries'] = excluded_libraries
            else:
                config['projects'][project_id]['excluded_libraries'] = set()

    return config


def make_config_json_serializable(original_config):
    """
    """
    config = original_config.copy()
    for project_id, project in config['projects'].items():
        excluded_runs = config['projects'][project_id]['excluded_runs']
        config['projects'][project_id]['excluded_runs'] = list(excluded_runs)
        excluded_libraries = config['projects'][project_id]['excluded_libraries']
        config['projects'][project_id]['excluded_libraries'] = list(excluded_libraries)

    return config
