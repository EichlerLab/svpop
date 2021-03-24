"""
Functions designed to directly support rules, such as input functions for Snakemake rules.
"""

import re

def sample_table_entry(name, sample_table, sample=None, wildcards=None):
    """
    Get an entry from sample_table (SAMPLE_TABLE read by `Snakefile`).

    :param name: Sample table record name.
    :param sample: Sample name. If not present, extract from "sample" entry in `wildcards`.
    :param wildcards: If present, parse the table entry with this value. Otherwise, parse with sample=sample.

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

    # Get seq_set
    if wildcards is not None and 'seq_set' in wildcards.keys():
        seq_set = wildcards.seq_set
    else:
        seq_set = 'DEFAULT'

    # Find table entry

    if (name, seq_set, sample) in sample_table.index:
        entry_sample = sample

    elif (name, seq_set, 'DEFAULT') in sample_table.index:
        entry_sample = 'DEFAULT'

    else:
        raise RuntimeError('Cannot find sample table entry for name "{}" and sample "{}" (or sample "DEFAULT") with seq-set {}'.format(
            name, sample, ('"' + seq_set + '"') if not pd.isnull(seq_set) else '<Undefined>'
        ))

    sample_entry = sample_table.loc[(name, seq_set, entry_sample)].copy()

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
                'Cannot get sample entry from sample table: Wildcards object is missing keys to fill format patterns in DATA: {}: Entry ({}, {}, {})'.format(
                    ', '.join(missing_wildcards),
                    name, seq_set, entry_sample
                )
            )

        sample_entry['DATA'] = sample_entry['DATA'].format(**wildcards)
    else:

        missing_wildcards = sample_entry['WILDCARDS'] - {'sample'}

        if missing_wildcards:
            raise RuntimeError(
                'Cannot get sample entry from sample table: No wildcards to fill format patterns in DATA: {}: Entry ({}, {}, {})'.format(
                    ', '.join(missing_wildcards),
                    name, seq_set, entry_sample
                )
            )

        sample_entry['DATA'] = sample_entry['DATA'].format(sample=sample)

    return sample_entry


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
