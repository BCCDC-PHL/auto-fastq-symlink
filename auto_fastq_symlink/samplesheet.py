import csv
import glob
import json
import logging
import os
import re
from typing import Optional

import jsonschema

import auto_fastq_symlink.util as util

log = logging.getLogger(__name__)


def _parse_header_section_miseq_v1(samplesheet_path):
    header_lines = []
    header = {}
    header['instrument_type'] = 'MiSeq'

    with open(samplesheet_path, 'r') as f:
        for line in f:
            if line.strip().startswith('[Header]'):
                continue
            if line.strip().startswith('[Reads]'):
                break
            else:
                header_lines.append(line.strip().rstrip(','))

    for line in header_lines:
        header_key = line.split(',')[0].lower().replace(" ", "_")

        if len(line.split(',')) > 1:
            header_value = line.split(',')[1]
        else:
            header_value = ""

        if header_key != "":
            header[header_key] = header_value
              
    return header


def _parse_reads_section_miseq_v1(samplesheet_path):
    """
    """
    reads_lines = []
    reads = []
    with open(samplesheet_path, 'r') as f:
        for line in f:
            if line.strip().startswith('[Reads]'):
                break
        for line in f:
            if line.strip().startswith('[Settings]'):
                break
            reads_lines.append(line.strip().rstrip(','))

    for line in reads_lines:
        if line != "":
            read_len = int(line.split(',')[0])
            reads.append(read_len)

    return reads


def _parse_settings_section_miseq_v1(samplesheet_path):
    """
    """
    settings_lines = []
    settings = {}
    with open(samplesheet_path, 'r') as f:
        for line in f:
            if line.strip().startswith('[Settings]'):
                break
        for line in f:
            if line.strip().startswith('[Data]'):
                break
            settings_lines.append(line.strip().rstrip(','))

    for line in settings_lines:
        settings_key = line.split(',')[0].lower().replace(" ", "_")

        if len(line.split(',')) > 1:
            settings_value = line.split(',')[1]
        else:
            settings_value = ""

        if settings_key != "":
            settings[settings_key] = settings_value
              
    return settings


def _parse_data_section_miseq_v1(samplesheet_path):
    """
    """
    data = []
    with open(samplesheet_path, 'r') as f:
        for line in f:
            if not line.strip().startswith('[Data]'):
                continue
            else:
                break
          
        data_header = [x.lower() for x in next(f).strip().split(',')]
      
        for line in f:
            if not all([x == '' for x in line.strip().split(',')]):
                data_line = {}
                for idx, data_element in enumerate(line.strip().split(',')):
                    try:
                        data_line[data_header[idx]] = data_element
                    except IndexError as e:
                        pass
                data.append(data_line)

    return data


def _determine_samplesheet_version(samplesheet_path, instrument_type):
    """
    """
    samplesheet_version = None
    if instrument_type == 'miseq':
        samplesheet_version = 1
    elif instrument_type == 'nextseq':
        samplesheet_version = 1

    return samplesheet_version


def _parse_samplesheet_miseq_v1(samplesheet_path: str) -> dict[str, object]:
    """
    :param samplesheet_path: Path to SampleSheet to be parsed
    :type samplesheet_path: str
    :return: The parsed SampleSheet
    :rtype: dict[str, object]
    :raises jsonschema.ValidationError: If the parsed samplesheed doesn't conform to the schema.
    """
    samplesheet = {}
    samplesheet['header'] = _parse_header_section_miseq_v1(samplesheet_path)
    samplesheet['reads'] = _parse_reads_section_miseq_v1(samplesheet_path)
    samplesheet['settings'] = _parse_settings_section_miseq_v1(samplesheet_path)
    samplesheet['data'] = _parse_data_section_miseq_v1(samplesheet_path)
    schema_path = os.path.join(os.path.dirname(__file__), "resources", "samplesheet_miseq_v1.schema.json")
    schema = None
    with open(schema_path, 'r') as f:
        schema = json.load(f)

    try:
        jsonschema.validate(instance=samplesheet, schema=schema)
    except jsonschema.ValidationError as e:
        log.error({"event_type": "samplesheet_validation_failed", "samplesheet_path": samplesheet_path, "schema_path": schema_path})
        raise e

    return samplesheet


def _read_sections_nextseq_v1(samplesheet_path: str) -> dict[str, list[str]]:
    """
    Split a NextSeq (or MiSeq i100) SampleSheet into its sections. Each section starts with a
    header line like `[Cloud_Data]`, and continues until the next section header. Blank lines
    and trailing commas are dropped.

    :param samplesheet_path: Path to the SampleSheet.
    :type samplesheet_path: str
    :return: Lines of each section, indexed by section name (eg. 'Cloud_Data').
    :rtype: dict[str, list[str]]
    """
    sections = {}
    section_name = None
    with open(samplesheet_path, 'r', encoding='utf-8-sig') as f:
        for line in f:
            line = line.strip().rstrip(',')
            section_header_match = re.fullmatch("\\[(.+)\\]", line)
            if section_header_match:
                section_name = section_header_match.group(1)
                sections.setdefault(section_name, [])
            elif section_name is not None and line != "":
                sections[section_name].append(line)

    return sections


def _parse_key_value_section_nextseq_v1(lines: list[str]) -> dict[str, str]:
    """
    Parse a section made up of `Key,Value` lines (eg. `[Header]` or `[BCLConvert_Settings]`).

    :param lines: Lines of the section.
    :type lines: list[str]
    :return: Values, indexed by key (converted to snake_case).
    :rtype: dict[str, str]
    """
    parsed = {}
    for line in lines:
        fields = line.split(',')
        key = util.camel_to_snake(fields[0])
        if key != "":
            parsed[key] = fields[1] if len(fields) > 1 else ""

    return parsed


def _parse_data_section_nextseq_v1(lines: list[str]) -> list[dict[str, str]]:
    """
    Parse a section made up of a line of column headers, followed by one line per record (eg. `[Cloud_Data]`).
    Records with fewer fields than there are column headers (which happens when trailing fields are empty) are
    padded with empty strings.

    :param lines: Lines of the section.
    :type lines: list[str]
    :return: Records, with column headers converted to snake_case.
    :rtype: list[dict[str, str]]
    """
    if len(lines) == 0:
        return []

    keys = [util.camel_to_snake(x) for x in lines[0].split(',')]
    records = []
    for line in lines[1:]:
        values = line.split(',')
        if all(value == "" for value in values):
            continue
        records.append({key: values[idx] if idx < len(values) else "" for idx, key in enumerate(keys)})

    return records


def _parse_samplesheet_nextseq_v1(samplesheet_path):
    """
    Parse a NextSeq SampleSheet. MiSeq i100 SampleSheets use the same format.
    Only the `[Cloud_Data]` section is required.

    :param samplesheet_path: Path to SampleSheet to be parsed
    :type samplesheet_path: str
    :return: The parsed SampleSheet
    :rtype: dict[str, object]
    :raises jsonschema.ValidationError: If the parsed samplesheet doesn't conform to the schema.
    """
    sections = _read_sections_nextseq_v1(samplesheet_path)

    samplesheet = {}
    samplesheet['header'] = {'instrument_type': 'NextSeq2000'}
    samplesheet['header'].update(_parse_key_value_section_nextseq_v1(sections.get('Header', [])))
    samplesheet['reads'] = {}
    for key, value in _parse_key_value_section_nextseq_v1(sections.get('Reads', [])).items():
        samplesheet['reads'][key] = int(value) if value.isdigit() else value
    samplesheet['sequencing_settings'] = _parse_key_value_section_nextseq_v1(sections.get('Sequencing_Settings', []))
    samplesheet['bclconvert_settings'] = _parse_key_value_section_nextseq_v1(sections.get('BCLConvert_Settings', []))
    samplesheet['bclconvert_data'] = _parse_data_section_nextseq_v1(sections.get('BCLConvert_Data', []))
    samplesheet['cloud_settings'] = _parse_key_value_section_nextseq_v1(sections.get('Cloud_Settings', []))
    samplesheet['cloud_data'] = _parse_data_section_nextseq_v1(sections.get('Cloud_Data', []))

    schema_path = os.path.join(os.path.dirname(__file__), "resources", "samplesheet_nextseq_v1.schema.json")
    schema = None
    with open(schema_path, 'r') as f:
        schema = json.load(f)

    try:
        jsonschema.validate(instance=samplesheet, schema=schema)
    except jsonschema.ValidationError as e:
        log.error({"event_type": "samplesheet_validation_failed", "samplesheet_path": samplesheet_path, "schema_path": schema_path})
        raise e

    return samplesheet


def parse_samplesheet_ont(samplesheet_path: str) -> dict[str, object]:
    """
    Parse an Oxford Nanopore (GridION or PromethION) sample sheet. It's a plain `.csv` file with one record
    per sample. The sample ID and project ID are combined in the `alias` field, as `<sample_id>_<project_id>`.

    :param samplesheet_path: Path to the sample sheet.
    :type samplesheet_path: str
    :return: The parsed sample sheet, with one record per sample under 'data'.
    :rtype: dict[str, object]
    :raises jsonschema.ValidationError: If the parsed sample sheet doesn't conform to the schema (eg. it doesn't have an `alias` field).
    """
    with open(samplesheet_path, 'r', encoding='utf-8-sig', newline='') as f:
        samplesheet = {'data': list(csv.DictReader(f))}

    schema_path = os.path.join(os.path.dirname(__file__), "resources", "samplesheet_ont_v1.schema.json")
    with open(schema_path, 'r') as f:
        schema = json.load(f)

    try:
        jsonschema.validate(instance=samplesheet, schema=schema)
    except jsonschema.ValidationError as e:
        log.error({"event_type": "samplesheet_validation_failed", "samplesheet_path": samplesheet_path, "schema_path": schema_path})
        raise e

    return samplesheet


def find_samplesheets(directory: str) -> list[str]:
    """
    Find all SampleSheet files (named like `SampleSheet*.csv`) directly inside a directory.

    :param directory: Directory to search.
    :type directory: str
    :return: List of paths to SampleSheet files.
    :rtype: list[str]
    """
    samplesheet_paths = sorted(glob.glob(os.path.join(directory, "SampleSheet*.csv")))

    return samplesheet_paths


def choose_samplesheet_to_parse(samplesheet_paths: list[str]) -> Optional[str]:
    """
    A directory may have multiple SampleSheet files in it. Choose only one to parse.
    A file named exactly `SampleSheet.csv` is preferred. Otherwise, if there is only one
    SampleSheet, it is chosen.

    :param samplesheet_paths: List of paths to SampleSheet files, all in the same directory.
    :type samplesheet_paths: list[str]
    :return: Path to the SampleSheet to parse, or None if one can't be chosen.
    :rtype: str | None
    """
    samplesheet_to_parse = None
    for samplesheet_path in samplesheet_paths:
        if os.path.basename(samplesheet_path) == "SampleSheet.csv":
            samplesheet_to_parse = samplesheet_path
    if not samplesheet_to_parse:
        if len(samplesheet_paths) == 1:
            samplesheet_to_parse = samplesheet_paths[0]
        else:
            # If there isn't a "SampleSheet.csv", and there are more than
            # one SampleSheet, then we have no other way of deciding which is preferable.
            pass

    return samplesheet_to_parse


def parse_samplesheet_miseq(samplesheet_path):
    """
    """
    samplesheet = None
    version = _determine_samplesheet_version(samplesheet_path, 'miseq')
    if version == 1:
        samplesheet = _parse_samplesheet_miseq_v1(samplesheet_path)

    return samplesheet


def parse_samplesheet_nextseq(samplesheet_path):
    """
    """
    samplesheet = None
    version = _determine_samplesheet_version(samplesheet_path, 'nextseq')
    if version == 1:
        samplesheet = _parse_samplesheet_nextseq_v1(samplesheet_path)

    return samplesheet


def parse_samplesheet(samplesheet_path: str, instrument_type: str) -> dict[str, object]:
    """
    :param samplesheet_path:
    :type samplesheet_path: str
    :param instrument_type: One of `miseq`, `nextseq`, `i100`, `gridion` or `promethion`
    :type instrument_type: str
    """
    samplesheet = None
    if instrument_type == 'miseq':
        samplesheet = parse_samplesheet_miseq(samplesheet_path)
    elif instrument_type in ('nextseq', 'i100'):
        # MiSeq i100 SampleSheets use the same format as NextSeq SampleSheets.
        samplesheet = parse_samplesheet_nextseq(samplesheet_path)
    elif instrument_type in ('gridion', 'promethion'):
        samplesheet = parse_samplesheet_ont(samplesheet_path)

    return samplesheet
    


