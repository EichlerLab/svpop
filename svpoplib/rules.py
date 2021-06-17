"""
Functions designed to directly support rules, such as input functions for Snakemake rules.
"""

import os
import pandas as pd
import re

import svpoplib

SAMPLE_TABLE_COL_TYPES = {
    'NAME': object,
    'SAMPLE': object,
    'TYPE': object,
    'DATA': object,
    'VERSION': object,
    'PARAMS': object
}


def get_sample_table(sample_table_file_name):
    """
    Read the sample table as a DataFrame, check format, and apply default values.

    :param sample_table_file_name: Sample table file name. If missing, a default empty table is returned.

    :return: Sample table DataFrame.
    """
    if os.path.isfile(sample_table_file_name):
        sample_table = pd.read_csv(sample_table_file_name, sep='\t', header=0, dtype=SAMPLE_TABLE_COL_TYPES)

        # Error on missing columns
        missing_columns = [col for col in ('NAME', 'TYPE', 'DATA') if col not in sample_table.columns]

        if missing_columns:
            raise RuntimeError('Missing sample table columns: {}'.format(', '.join(missing_columns)))

        if 'SAMPLE' not in sample_table.columns:
            sample_table['SAMPLE'] = 'DEFAULT'

        sample_table['SAMPLE'] = sample_table['SAMPLE'].fillna('DEFAULT')

        if 'VERSION' not in sample_table.columns:
            sample_table['VERSION'] = np.nan

        if 'PARAMS' not in sample_table.columns:
            sample_table['PARAMS'] = np.nan

        sample_table.set_index(['NAME', 'SAMPLE'], inplace=True, drop=False)
    else:
        sample_table = pd.DataFrame(
            [], columns=['NAME', 'SAMPLE', 'TYPE', 'DATA', 'VERSION', 'PARAMS']
        ).set_index(
            ['NAME', 'SAMPLE'], drop=False
        )

    return sample_table


def sample_table_entry(name, sample_table, sample=None, wildcards=None, type=None):
    """
    Get an entry from sample_table (SAMPLE_TABLE read by `Snakefile`).

    :param name: Sample table record name.
    :param sample_table: Sample table DataFrame.
    :param sample: Sample name. If not present, extract from "sample" entry in `wildcards`.
    :param wildcards: If present, parse the table entry with this value. Otherwise, parse with sample=sample.
    :param type: Entry must be declared as this type (points to parser used to process the input data).

    :return: A configuration entry (Pandas.Series).
    """

    # Check name
    if name is None:
        raise RuntimeError('Cannot get input for entry name: None')

    # Get sample name
    if sample is None:
        if wildcards is None:
            raise RuntimeError('Cannot get global sample table entry: sample and wildcards are both None')

        if 'sample' not in wildcards.keys():
            raise RuntimeError('Cannot get sample from wildcards: "sample" is not a wildcards key')

        sample = wildcards.sample

    # Find table entry

    if (name, sample) in sample_table.index:
        entry_sample = sample

    elif (name, 'DEFAULT') in sample_table.index:
        entry_sample = 'DEFAULT'

    else:
        raise RuntimeError('Cannot find sample table entry for name "{}" and sample "{}" (or sample "DEFAULT")'.format(
            name, sample
        ))

    sample_entry = sample_table.loc[(name, entry_sample)].copy()

    # Check type
    if type is not None:
        if sample_entry['TYPE'] != type:
            raise RuntimeError('Source type mismatch for sample entry ({}, {}): Expected type "{}", entry has type "{}"'.format(
                name, sample, type, sample_entry['TYPE']
            ))

    # Save wildcards (some rules need to know which wildcards were replaced)
    # Track sample name (actual sample name or DEFAULT) an the original DATA pattern.
    sample_entry['WILDCARDS'] = set(re.findall(r'\{([^\{\}]+)\}', sample_entry['DATA']))
    sample_entry['ENTRY_SAMPLE'] = entry_sample
    sample_entry['DATA_PATTERN'] = sample_entry['DATA']

    # Replace wildcards
    if wildcards is not None:

        missing_wildcards = sample_entry['WILDCARDS'] - set(wildcards.keys())

        if missing_wildcards:
            raise RuntimeError(
                'Cannot get sample entry from sample table: Wildcards object is missing keys to fill format patterns in DATA: {}: Entry ({}, {})'.format(
                    ', '.join(missing_wildcards),
                    name, entry_sample
                )
            )

        sample_entry['DATA'] = sample_entry['DATA'].format(**wildcards)
    else:

        missing_wildcards = sample_entry['WILDCARDS'] - {'sample'}

        if missing_wildcards:
            raise RuntimeError(
                'Cannot get sample entry from sample table: No wildcards to fill format patterns in DATA: {}: Entry ({}, {})'.format(
                    ', '.join(missing_wildcards),
                    name, entry_sample
                )
            )

        sample_entry['DATA'] = sample_entry['DATA'].format(sample=sample)

    # Set params
    sample_entry['PARAM_STRING'] = sample_entry['PARAMS']
    sample_entry['PARAMS'] = svpoplib.util.parse_param_string(sample_entry['PARAM_STRING'])

    return sample_entry


def get_bed_fa_input(sample_entry, wildcards, default=None):
    """
    Locate FASTA sequence input file.

    :param sample_entry: Entry from the sample table.

    :return: FASTA input file location
    """

    # Check for explicit location
    if 'fa_pattern' in sample_entry['PARAMS']:
        if sample_entry['PARAMS'] is None or not sample_entry['PARAMS'].strip():
            return default

        return sample_entry['PARAMS'].format(**wildcards)

    fa_file_name = os.path.join(
        os.path.dirname(sample_entry['DATA']),
        'fa',
        '{vartype}_{svtype}.fa.gz'.format(**wildcards)
    )

    if os.path.isfile(fa_file_name):
        return fa_file_name

    return default


def parse_wildcards(file_pattern, name, sample_table, sample=None, wildcards=None, type=None):
    """
    Parse a file pattern with wildcards derived from rule wildcards and the sample config.

    :param file_pattern: Pattern to parse. Format patterns are filled by wildcards with "sourcename"
        and "source_type" derived from the sample entry.
    :param name: Sample table record name.
    :param sample_table: Sample table DataFrame.
    :param sample: Sample name. If not present, extract from "sample" entry in `wildcards`.
    :param wildcards: If present, parse the table entry with this value. Otherwise, parse with sample=sample.
    :param type: Entry must be declared as this type (points to parser used to process the input data).

    :return: `file_pattern` with patterns filled in.
    """

    # Get sample entry
    sample_entry = sample_table_entry(name=name, sample_table=sample_table, sample=sample, wildcards=wildcards, type=type)

    # Set wildcards
    wildcards = dict(wildcards)

    wildcards['sourcename'] = sample_entry['NAME']
    wildcards['source_type'] = sample_entry['TYPE']

    return file_pattern.format(**wildcards)


# def variant_global_sample_info_entry(sample):
#     """
#     Get an entry from SAMPLE_INFO_TABLE.
#
#     :param sample: Sample name.
#
#     :return: A configuration entry (Pandas.Series) or None if there is no sample information.
#     """
#
#     if sample not in SAMPLE_INFO_TABLE.index:
#         return None
#
#     return SAMPLE_INFO_TABLE.loc[sample]
#
# def variant_global_sample_info_entry_element(sample, element, default=None, null_default=True):
#     """
#     Get an entry from SAMPLE_INFO_TABLE.
#
#     :param sample: Sample name (sample info table row).
#     :param element: Element name (sample info table column).
#     :param default: Default value if sample or element is not defined.
#     :param null_default: If the element is present but is null (`np.nan`), then return `default` instead of `np.nan`.
#
#     :return: A configuration entry (Pandas.Series) or None if there is no sample information.
#     """
#
#     # Get entry
#     sample_info = variant_global_sample_info_entry(sample)
#
#     if sample_info is None:
#         return default
#
#     # Get element value
#     if element not in sample_info:
#         return default
#
#     value = sample_info[element]
#
#     # Translate null to default (if null_default)
#     if null_default and pd.isnull(value):
#         return default
#
#     # Return value
#     return value
