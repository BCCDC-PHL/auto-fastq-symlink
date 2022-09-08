import csv
import json
import logging

def parse_list(list_path):
    items = []
    with open(list_path, 'r') as f:
        for line in f:
            item = line.strip()
            items.append(item)

    return items


def parse_projects_definition_file(projects_definition_file_path):
    projects = {}
    with open(projects_definition_file_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            projects[row['project_id']] = row

    return projects


def load_config(config_path):
    config = None
    with open(config_path, 'r') as f:
        config = json.load(f)

    if 'projects_definition_file' in config:
        projects = parse_projects_definition_file(config['projects_definition_file'])
        config['projects'] = projects
        for project_id, project in projects.items():
            if 'excluded_runs_list' in project:
                try:
                    excluded_runs_list = parse_list(project['excluded_runs_list'])
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
                    excluded_libraries_list = parse_list(project['excluded_libraries_list'])
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
    config = original_config.copy()
    for project_id, project in config['projects'].items():
        excluded_runs = config['projects'][project_id]['excluded_runs']
        config['projects'][project_id]['excluded_runs'] = list(excluded_runs)
        excluded_libraries = config['projects'][project_id]['excluded_libraries']
        config['projects'][project_id]['excluded_libraries'] = list(excluded_libraries)

    return config
