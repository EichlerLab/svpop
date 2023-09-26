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

    return sample_table.sort_index().copy()


def sample_table_entry(name, sample_table, sample=None, wildcards=None, caller_type=None, format_data=True):
    """
    Get an entry from sample_table (SAMPLE_TABLE read by `Snakefile`).

    :param name: Sample table record name.
    :param sample_table: Sample table DataFrame.
    :param sample: Sample name. If not present, extract from "sample" entry in `wildcards`.
    :param wildcards: If present, parse the table entry with this value. Otherwise, parse with sample=sample.
    :param caller_type: Entry must be declared as this type (points to parser used to process the input data).
    :param format_data: If True, parse wildcards into 'DATA' entry. If False, leave 'DATA' with wildcards in it.

    :return: A configuration entry (Pandas.Series).
    """

    # Get wildcards
    if wildcards is None:
        wildcards = dict()
        wildcards_is_none = True

    else:
        wildcards = dict(wildcards)
        wildcards_is_none = False

    # Expand wildcard aliases
    if 'varsvtype' in wildcards.keys():
        varsvtype_tok = wildcards['varsvtype'].split('_')

        if len(varsvtype_tok) != 2:
            raise RuntimeError('Found wildcard "varsvtype" with value "{}", but expected two "_" separated fields'.format(wildcards['varsvtype']))

        if 'vartype' not in wildcards.keys():
            wildcards['vartype'] = varsvtype_tok[0]

        if 'svtype' not in wildcards.keys():
            wildcards['svtype'] = varsvtype_tok[1]

    # Check name
    if name is None:
        raise RuntimeError('Cannot get input for entry name: None')

    # Get sample name
    if sample is None:
        sample = wildcards.get('sample', None)

    if sample is None:
        if wildcards_is_none:
            raise RuntimeError('Cannot get global sample table entry: sample and wildcards are both None')

        raise RuntimeError('Cannot get sample from wildcards: "sample" is not a wildcards key')

    # Find table entry
    if (name, sample) in sample_table.index:
        entry_sample = sample

    elif (name, 'DEFAULT') in sample_table.index:
        entry_sample = 'DEFAULT'

    else:
        raise RuntimeError('Cannot find sample table entry for name "{}" and sample "{}" (or sample "DEFAULT")'.format(
            name, sample
        ))

    sample_entry = sample_table.loc[[(name, entry_sample)]]

    if sample_entry.shape[0] == 0:
        raise RuntimeError(f'No entries in sample TSV: NAME="{name}", SAMPLE="{sample}"')

    if sample_entry.shape[0] > 1:
        raise RuntimeError(f'Multiple entries in sample TSV: NAME="{name}", SAMPLE="{sample}": rows={sample_entry.shape[0]}')

    sample_entry = sample_entry.iloc[0].copy()

    # Set type field
    sample_entry_type = sample_entry['TYPE']

    if not pd.isnull(sample_entry_type):
        sample_entry_type = sample_entry_type.strip()

    if pd.isnull(sample_entry_type) or sample_entry_type == '':
        raise RuntimeError(f'Emtpy TYPE for table entry: NAME="{name}", SAMPLE="{sample}"')

    sample_entry['TYPE'] = sample_entry['TYPE'].strip().lower()

    # Check type
    if caller_type is not None:
        if sample_entry['TYPE'] != caller_type:
            raise RuntimeError('Source type mismatch for sample entry ({}, {}): Expected type "{}", entry has type "{}"'.format(
                name, sample, caller_type, sample_entry['TYPE']
            ))

    # Save wildcards (some rules need to know which wildcards were replaced)
    # Track sample name (actual sample name or DEFAULT) an the original DATA pattern.
    sample_entry['WILDCARDS'] = set(re.findall(r'\{([^\{\}]+)\}', sample_entry['DATA']))
    sample_entry['ENTRY_SAMPLE'] = entry_sample
    sample_entry['DATA_PATTERN'] = sample_entry['DATA']

    # Replace wildcards in DATA
    if format_data:
        missing_wildcards = sample_entry['WILDCARDS'] - set(wildcards.keys())

        if missing_wildcards:
            raise RuntimeError(
                'Cannot get sample entry from sample table: Wildcards object is missing keys to fill format patterns in DATA: {}: Entry ({}, {})'.format(
                    ', '.join(missing_wildcards),
                    name, entry_sample
                )
            )

        sample_entry['DATA'] = sample_entry['DATA'].format(**wildcards)

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
        if sample_entry['PARAMS'] is None or not sample_entry['PARAMS']['fa_pattern'].strip():
            return default

        return sample_entry['PARAMS']['fa_pattern'].format(**wildcards)

    fa_file_name = os.path.join(
        os.path.dirname(sample_entry['DATA']),
        'fa',
        '{vartype}_{svtype}.fa.gz'.format(**wildcards)
    )

    if os.path.isfile(fa_file_name):
        return fa_file_name

    return default


def parse_wildcards(file_pattern, name, sample_table, sample=None, wildcards=None, caller_type=None):
    """
    Parse a file pattern with wildcards derived from rule wildcards and the sample config.

    :param file_pattern: Pattern to parse. Format patterns are filled by wildcards with "sourcename"
        and "callertype" derived from the sample entry.
    :param name: Sample table record name.
    :param sample_table: Sample table DataFrame.
    :param sample: Sample name. If not present, extract from "sample" entry in `wildcards`.
    :param wildcards: If present, parse the table entry with this value. Otherwise, parse with sample=sample.
    :param caller_type: Entry must be declared as this type (points to parser used to process the input data).

    :return: `file_pattern` with patterns filled in.
    """

    # Get sample entry
    sample_entry = sample_table_entry(name=name, sample_table=sample_table, sample=sample, wildcards=wildcards, caller_type=caller_type)

    # Set wildcards
    wildcards = dict(wildcards)

    wildcards['sourcename'] = sample_entry['NAME']
    wildcards['callertype'] = sample_entry['TYPE']

    return file_pattern.format(**wildcards)


def get_sample_list(sample_list_name, config):
    """
    Get a named list of samples ("samplelist" in config) or a single sample name if the sample-list is not defined.

    :param sample_list_name: Name of the sample list to retrieve.
    :param config: Config object.

    :return: List of sample names or None if the sample list was not found.
    """

    sample_list_name = sample_list_name.strip()

    if not sample_list_name:
        raise RuntimeError('Sample list name is empty')

    # Get attributes
    hap_expand = False
    expand_match = False

    if ':' in sample_list_name:
        tok = sample_list_name.split(':')

        list_name = tok[0].strip()

        for attr in tok[1:]:
            attr = attr.strip()

            if not attr:
                continue

            if attr == 'hap':
                hap_expand = True

            elif attr == 'hapmatch':
                expand_match = True

            else:
                raise RuntimeError(f'Unrecognized attribute "{attr}" in sample list name {sample_list_name}')

    else:
            list_name = sample_list_name

    # Check attributes
    if hap_expand and expand_match:
        raise RuntimeError(f'Conflicting sample list attributes: Cannot set "hapexp" and "expmatch" for the same list')

    # Get list
    sample_list = config.get('samplelist', dict()).get(list_name, [list_name])

    # Process attributes
    if hap_expand and sample_list is not None:
        exp_sample_list = list()

        for sample_name in sample_list:
            exp_sample_list.append(f'{sample_name}-h1')
            exp_sample_list.append(f'{sample_name}-h2')

        sample_list = exp_sample_list

    if expand_match and sample_list is not None:
        exp_sample_list = list()

        for sample_name in sample_list:
            exp_sample_list.append(f'{sample_name}')
            exp_sample_list.append(f'{sample_name}')

        sample_list = exp_sample_list

    # Return list
    return sample_list
