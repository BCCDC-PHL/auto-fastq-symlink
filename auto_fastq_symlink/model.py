import datetime

from sqlalchemy import Column
from sqlalchemy import ForeignKey
from sqlalchemy import Boolean
from sqlalchemy import Integer
from sqlalchemy import Float
from sqlalchemy import String
from sqlalchemy import Date
from sqlalchemy import DateTime
from sqlalchemy import JSON

from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()


class SequencingRun(Base):
    __tablename__ = 'sequencing_run'

    sequencing_run_id = Column(String, primary_key=True)
    instrument_type = Column(String)
    instrument_id = Column(String)
    run_date = Column(Date)
    run_directory = Column(String)
    fastq_directory = Column(String)
    samplesheet = Column(String)
    timestamp_updated = Column(DateTime, default=datetime.datetime.now, onupdate=datetime.datetime.now)

    libraries = relationship(
        "Library", back_populates="sequencing_run", cascade="all, delete-orphan"
    )


class Library(Base):
    __tablename__ = 'library'

    library_id = Column(String, primary_key=True)
    sequencing_run_id = Column(String, ForeignKey("sequencing_run.sequencing_run_id"), primary_key=True)
    project_id = Column(String, ForeignKey("project.project_id"))
    fastq_path_r1 = Column(String)
    fastq_path_r2 = Column(String)
    timestamp_updated = Column(DateTime, default=datetime.datetime.now, onupdate=datetime.datetime.now)

    sequencing_run = relationship("SequencingRun", back_populates="libraries")
    project = relationship("Project", back_populates="libraries")


class Project(Base):
    __tablename__ = 'project'

    project_id = Column(String, primary_key=True)
    fastq_symlinks_directory = Column(String)
    timestamp_updated = Column(DateTime, default=datetime.datetime.now, onupdate=datetime.datetime.now)

    libraries = relationship(
        "Library", back_populates="project", cascade="all, delete-orphan"
    )


class Symlink(Base):
    __tablename__ = 'symlink'

    project_id = Column(Integer, ForeignKey("project.project_id"))
    library_id = Column(Integer, ForeignKey("library.library_id"))
    path = Column(String, primary_key=True)
    target = Column(String, primary_key=True)
    timestamp_updated = Column(DateTime, default=datetime.datetime.now, onupdate=datetime.datetime.now)
