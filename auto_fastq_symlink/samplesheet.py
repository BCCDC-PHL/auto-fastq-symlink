import csv
import glob
import json
import logging
import os
import re

import jsonschema

import auto_fastq_symlink.util as util


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


def _parse_samplesheet_miseq_v1(samplesheet_path):
    """
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
        logging.error({"event_type": "samplesheet_validation_failed", "samplesheet_path": samplesheet_path, "schema_path": schema_path})

    return samplesheet


def _parse_header_section_nextseq_v1(samplesheet_path):
    """
    """
    header_lines = []
    header = {}
    header['instrument_type'] = 'NextSeq2000'

    with open(samplesheet_path, 'r') as f:
        for line in f:
            if line.strip().startswith('[Header]'):
                continue
            if line.strip().startswith('[Reads]') or line.strip().rstrip(',') == "":
                break
            else:
                header_lines.append(line.strip().rstrip(','))

    for line in header_lines:
        header_key = util.camel_to_snake(line.split(',')[0])

        if len(line.split(',')) > 1:
            header_value = line.split(',')[1]
        else:
            header_value = ""

        if header_key != "":
            header[header_key] = header_value
              
    return header


def _parse_reads_section_nextseq_v1(samplesheet_path):
    """
    """
    reads_lines = []
    reads = {}
    with open(samplesheet_path, 'r') as f:
        for line in f:
            if line.strip().startswith('[Reads]'):
                break
        for line in f:
            if line.strip().startswith('[Sequencing_Settings]') or line.strip().rstrip(',') == "":
                break
            reads_lines.append(line.strip().rstrip(','))

    for line in reads_lines:
        reads_key = util.camel_to_snake(line.split(',')[0])
        if len(line.split(',')) > 1:
            reads_value = int(line.split(',')[1])
        else:
            reads_value = ""

        if reads_key != "":
            reads[reads_key] = reads_value

    return reads


def _parse_sequencing_settings_section_nextseq_v1(samplesheet_path):
    """
    """
    sequencing_settings_lines = []
    sequencing_settings = {}
    with open(samplesheet_path, 'r') as f:
        for line in f:
            if line.strip().startswith('[Sequencing_Settings]'):
                break
        for line in f:
            if line.strip().startswith('[BCLConvert_Settings]') or line.strip().rstrip(',') == "":
                break
            sequencing_settings_lines.append(line.strip().rstrip(','))

    for line in sequencing_settings_lines:
        sequencing_settings_key = util.camel_to_snake(line.split(',')[0])

        if len(line.split(',')) > 1:
            sequencing_settings_value = line.split(',')[1]
        else:
            sequencing_settings_value = ""

        if sequencing_settings_key != "":
            sequencing_settings[sequencing_settings_key] = sequencing_settings_value
              
    return sequencing_settings


def _parse_bclconvert_settings_section_nextseq_v1(samplesheet_path):
    """
    """
    bclconvert_settings_lines = []
    bclconvert_settings = {}
    with open(samplesheet_path, 'r') as f:
        for line in f:
            if line.strip().startswith('[BCLConvert_Settings]'):
                break
        for line in f:
            if line.strip().startswith('[BCLConvert_Data]') or line.strip().rstrip(',') == "":
                break
            bclconvert_settings_lines.append(line.strip().rstrip(','))

    for line in bclconvert_settings_lines:
        bclconvert_settings_key = util.camel_to_snake(line.split(',')[0])

    if len(line.split(',')) > 1:
        bclconvert_settings_value = line.split(',')[1]
    else:
        bclconvert_settings_value = ""

    if bclconvert_settings_key != "":
        bclconvert_settings[bclconvert_settings_key] = bclconvert_settings_value
              
    return bclconvert_settings


def _parse_bclconvert_data_section_nextseq_v1(samplesheet_path):
    """
    """
    bclconvert_data_lines = []
    bclconvert_data = []
    with open(samplesheet_path, 'r') as f:
        for line in f:
            if line.strip().startswith('[BCLConvert_Data]'):
                break
        for line in f:
            if line.strip().startswith('[Cloud_Settings]') or line.strip().rstrip(',') == "":
                break
            bclconvert_data_lines.append(line.strip().rstrip(','))

    bclconvert_data_keys = [util.camel_to_snake(x) for x in bclconvert_data_lines[0].split(',')]
    for line in bclconvert_data_lines[1:]:
        d = {}
        values = line.split(',')
        for idx, key in enumerate(bclconvert_data_keys):
            d[key] = values[idx]
        bclconvert_data.append(d)

    return bclconvert_data


def _parse_cloud_settings_section_nextseq_v1(samplesheet_path):
    """
    """
    cloud_settings_lines = []
    cloud_settings = {}
    with open(samplesheet_path, 'r') as f:
        for line in f:
            if line.strip().startswith('[Cloud_Settings]'):
                break
        for line in f:
            if line.strip().startswith('[BCLConvert_Settings]') or line.strip().rstrip(',') == "":
                break
            cloud_settings_lines.append(line.strip().rstrip(','))

    for line in cloud_settings_lines:
        cloud_settings_key = util.camel_to_snake(line.split(',')[0])

        if len(line.split(',')) > 1:
            cloud_settings_value = line.split(',')[1]
        else:
            cloud_settings_value = ""

        if cloud_settings_key != "":
            cloud_settings[cloud_settings_key] = cloud_settings_value
              
    return cloud_settings


def _parse_cloud_data_section_nextseq_v1(samplesheet_path):
    cloud_data_lines = []
    cloud_data = []
    with open(samplesheet_path, 'r') as f:
        for line in f:
            if line.strip().startswith('[Cloud_Data]'):
                break
        for line in f:
            if line.strip().startswith('[Cloud_Settings]') or line.strip().rstrip(',') == "":
                break
            cloud_data_lines.append(line.strip().rstrip(','))

    if cloud_data_lines:
        cloud_data_keys = [util.camel_to_snake(x) for x in cloud_data_lines[0].split(',')]
        for line in cloud_data_lines[1:]:
            d = {}
            values = line.strip().split(',')
            if not all([x == '' for x in values]):
                for idx, key in enumerate(cloud_data_keys):
                    try:
                        d[key] = values[idx]
                    except IndexError as e:
                        d[key] = ""
                cloud_data.append(d)
  
    return cloud_data



def _parse_samplesheet_nextseq_v1(samplesheet_path):
    """
    """
    samplesheet = {}
    samplesheet['header'] = _parse_header_section_nextseq_v1(samplesheet_path)
    samplesheet['reads'] = _parse_reads_section_nextseq_v1(samplesheet_path)
    samplesheet['sequencing_settings'] = _parse_sequencing_settings_section_nextseq_v1(samplesheet_path)
    samplesheet['bclconvert_settings'] = _parse_bclconvert_settings_section_nextseq_v1(samplesheet_path)
    samplesheet['bclconvert_data'] = _parse_bclconvert_data_section_nextseq_v1(samplesheet_path)
    samplesheet['cloud_settings'] = _parse_cloud_settings_section_nextseq_v1(samplesheet_path)
    samplesheet['cloud_data'] = _parse_cloud_data_section_nextseq_v1(samplesheet_path)

    return samplesheet


def find_samplesheets(run_dir, instrument_type):
    """
    """
    samplesheet_paths = None
    if instrument_type == 'miseq':
        samplesheet_paths = glob.glob(os.path.join(run_dir, "SampleSheet*.csv"))
    elif instrument_type == 'nextseq':
        top_level_samplesheets = glob.glob(os.path.join(run_dir, "SampleSheet*.csv"))
        analysis_dir_samplesheets = glob.glob(os.path.join(run_dir, "Analysis", "*", "Data", "SampleSheet*.csv"))
        samplesheet_paths = top_level_samplesheets + analysis_dir_samplesheets                                    
    return samplesheet_paths


def choose_samplesheet_to_parse(samplesheet_paths: list[str], instrument_type: str, run_id: str):
    """
    A run directory may have multiple SampleSheet.csv files in it. Choose only one to parse.

    :param samplesheet_paths: List of paths to SampleSheet.csv files
    :type samplesheet_paths: list[str]
    :param instrument_type: Instrument type, should be one of: "miseq", "nextseq"
    :type instrument_type: str
    :param run_id: Sequencing run ID
    :type run_id: str
    """
    samplesheet_to_parse = None
    if instrument_type == 'miseq':
        for samplesheet_path in samplesheet_paths:
            if re.match("SampleSheet\.csv", os.path.basename(samplesheet_path)):
                samplesheet_to_parse = samplesheet_path
        if not samplesheet_to_parse:
            if len(samplesheet_paths) == 1:
                samplesheet_to_parse = samplesheet_paths[0]
            else:
                # If there isn't a top-level "SampleSheet.csv", and there are more than
                # one SampleSheet, then we have no other way of deciding which is preferable.
                pass
    elif instrument_type == 'nextseq':
        samplesheets_by_analysis_num = {}
        for samplesheet_path in samplesheet_paths:
            match = re.search("Analysis/(\d+)/Data", samplesheet_path)
            if match:
                analysis_num = int(match.group(1))
                samplesheets_by_analysis_num[analysis_num] = samplesheet_path
        largest_analysis_num = 0
        for analysis_num, samplesheet_path in samplesheets_by_analysis_num.items():
            if analysis_num > largest_analysis_num:
                largest_analysis_num = analysis_num
        if largest_analysis_num > 0:
            samplesheet_to_parse = samplesheets_by_analysis_num[largest_analysis_num]
        if not samplesheet_to_parse:
            if len(samplesheet_paths) == 1:
                samplesheet_to_parse = samplesheet_paths[0]
            else:
                # If there isn't a top-level "SampleSheet.csv", and there are more than
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


def parse_samplesheet(samplesheet_path, instrument_type):
    """
    """
    samplesheet = None
    if instrument_type == 'miseq':
        samplesheet = parse_samplesheet_miseq(samplesheet_path)
    elif instrument_type == 'nextseq':
        samplesheet = parse_samplesheet_nextseq(samplesheet_path)

    return samplesheet
    


