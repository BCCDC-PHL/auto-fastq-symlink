import datetime
import json
import logging
import os

from sqlalchemy.orm import Session, sessionmaker
from sqlalchemy import create_engine

import auto_fastq_symlink.util as util

from auto_fastq_symlink.model import *


def store_projects(config, projects):
    """
    """
    connection_uri = config['database_connection_uri']
    engine = create_engine(connection_uri)
    Session = sessionmaker(bind=engine)
    session = Session()

    existing_projects = session.query(Project).all()
    existing_project_ids = set([project.project_id for project in existing_projects])

    projects_to_store = []
    for project_id, project in projects.items():
        if project_id not in existing_project_ids:
            p = Project(
                project_id = project_id,
                fastq_symlinks_directory = project['fastq_symlinks_dir']
            )
            projects_to_store.append(p)

    for project in projects_to_store:
        logging.debug(json.dumps({"event_type": "project_stored", "project_id": project.project_id}))

    session.add_all(projects_to_store)
    session.commit()


def store_libraries(session: Session, run: dict[str, object]):
    """
    """
    run_id = run['run_id']
    
    libraries_to_store = []
    for library in run['libraries']:
        if library['project_id'] == "":
            project_id = None
        else:
            project_id = library['project_id']

        logging.debug(json.dumps({
            "event_type": "preparing_to_queue_library_for_storage",
            "sequencing_run_id": run_id,
            "library_id": library['library_id'],
            "project_id": library['project_id'],
            "fastq_path_r1": library['fastq_path_r1'],
            "fastq_path_r2": library['fastq_path_r2'],
        }))
        l = Library(
            library_id = library['library_id'],
            sequencing_run_id = run_id,
            project_id = project_id,
            fastq_path_r1 = library['fastq_path_r1'],
            fastq_path_r2 = library['fastq_path_r2'],
        )
        libraries_to_store.append(l)
        logging.debug(json.dumps({
            "event_type": "queued_library_for_storage",
            "sequencing_run_id": run_id,
            "library_id": library['library_id'],
            "project_id": library['project_id'],
            "fastq_path_r1": library['fastq_path_r1'],
            "fastq_path_r2": library['fastq_path_r2']
        }))

    session.add_all(libraries_to_store)
    session.commit()
    logging.debug(json.dumps({"event_type": "run_libraries_stored", "sequencing_run_id": run_id}))


def update_libraries(session: Session, run: dict[str, object]):
    """
    """
    run_id = run['run_id']

    for library in run['libraries']:
        existing_library = session.query(Library).filter(
            Library.sequencing_run_id == run_id,
            Library.library_id == library['library_id']).first()

        if existing_library:
            existing_library.library_id = library['library_id']
            existing_library.project_id = library['project_id']
            existing_library.fastq_path_r1 = library['fastq_path_r1']
            existing_library.fastq_path_r2 = library['fastq_path_r2']
            session.commit()
            logging.debug(json.dumps({"event_type": "library_updated", "sequencing_run_id": run_id, "library_id": library['library_id']}))
            

def store_run(config: dict[str, object], run: dict[str, object]):
    """
    """

    connection_uri = str(config['database_connection_uri'])
    engine = create_engine(connection_uri)
    Session = sessionmaker(bind=engine)
    session = Session()

    run_id = run['run_id']
    existing_run = session.query(SequencingRun).filter(SequencingRun.sequencing_run_id == run_id).first()

    six_digit_date = run_id.split('_')[0]
    year = int("20" + six_digit_date[0:2])
    month = int(six_digit_date[2:4])
    day = int(six_digit_date[4:6])
    run_date = datetime.date(year, month, day)

    instrument_id = run_id.split('_')[1]

    if existing_run:
        existing_run.samplesheet = run['parsed_samplesheet']
        existing_run.fastq_directory = run['fastq_directory']
        session.commit()
        logging.debug(json.dumps({"event_type": "run_updated", "run_id": run_id}))
        update_libraries(session, run)
    else:
        r = SequencingRun(
            sequencing_run_id = run_id,
            run_date = run_date,
            instrument_type = run['instrument_type'],
            instrument_id = instrument_id,
            samplesheet = run['parsed_samplesheet'],
            run_directory = run['run_directory'],
            fastq_directory = run['fastq_directory'],
        )
        session.add(r)
        store_libraries(session, run)
        logging.debug(json.dumps({"event_type": "run_stored", "run_id": run_id}))


def store_symlinks(config: dict[str, object], symlinks_by_project_id: dict[str, object]):
    """
    """
    connection_uri = str(config['database_connection_uri'])
    engine = create_engine(connection_uri)
    Session = sessionmaker(bind=engine)
    session = Session()

    existing_symlinks = session.query(Symlink).all()
    existing_symlink_path_target_tuples = set([(symlink.path, symlink.target) for symlink in existing_symlinks])

    symlinks_to_store = []
    for project_id, symlinks in symlinks_by_project_id.items():
        for symlink in symlinks:
            path_target_tuple = (symlink['path'], symlink['target'])
            if path_target_tuple not in existing_symlink_path_target_tuples:
                library_id = os.path.basename(symlink['target']).split('_')[0]
                s = Symlink(
                    project_id = project_id,
                    sequencing_run_id = symlink['sequencing_run_id'],
                    library_id = library_id,
                    path = symlink['path'],
                    target = symlink['target'],
                )
                symlinks_to_store.append(s)

    session.add_all(symlinks_to_store)
    session.commit()


def delete_nonexistent_symlinks(config):
    """
    """
    connection_uri = config['database_connection_uri']
    engine = create_engine(connection_uri)
    Session = sessionmaker(bind=engine)
    session = Session()

    existing_symlinks = session.query(Symlink).all()

    for symlink in existing_symlinks:
        if not os.path.exists(symlink.path):
            session.delete(symlink)
            session.commit()
    

def get_symlinks(config):
    """
    """
    connection_uri = config['database_connection_uri']
    engine = create_engine(connection_uri)
    Session = sessionmaker(bind=engine)
    session = Session()

    query_result = session.query(Symlink).all()

    existing_symlinks = []
    for row in query_result:
        existing_symlinks.append(util.row2dict(row))
        

    return existing_symlinks


def get_symlinks_by_run_id(config, run_id):
    """
    """
    connection_uri = config['database_connection_uri']
    engine = create_engine(connection_uri)
    Session = sessionmaker(bind=engine)
    session = Session()

    query_result = session.query(Symlink).filter(Symlink.sequencing_run_id == run_id)

    existing_symlinks_for_run = []
    for row in query_result:
        existing_symlinks_for_run.append(util.row2dict(row))
        

    return existing_symlinks_for_run


def get_libraries(config):
    """
    """
    connection_uri = config['database_connection_uri']
    engine = create_engine(connection_uri)
    Session = sessionmaker(bind=engine)
    session = Session()

    query_result = session.query(Library).all()

    all_libraries = []
    for row in query_result:
        all_libraries.append(util.row2dict(row))

    return all_libraries


def get_libraries_by_project_id(config, project_id):
    """
    """
    connection_uri = config['database_connection_uri']
    engine = create_engine(connection_uri)
    Session = sessionmaker(bind=engine)
    session = Session()

    query_result = session.query(Library).filter(Library.project_id == project_id)

    project_libraries = []
    for row in query_result:
        project_libraries.append(util.row2dict(row))

    return project_libraries


def get_libraries_by_project_id_and_run_id(config, project_id, run_id):
    """
    """
    connection_uri = config['database_connection_uri']
    engine = create_engine(connection_uri)
    Session = sessionmaker(bind=engine)
    session = Session()

    query_result = session.query(Library).filter(
        Library.project_id == project_id,
        Library.sequencing_run_id == run_id,
    )

    project_libraries = []
    for row in query_result:
        project_libraries.append(util.row2dict(row))

    return project_libraries

    
    
