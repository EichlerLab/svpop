"""
Intersect support routines.
"""

import numpy as np
import pandas as pd
import svpoplib

#
# Definitions
#

MERGE_INFO_FIELD_LIST = [
    'MERGE_OFFSET', 'MERGE_RO', 'MERGE_SZRO', 'MERGE_OFFSZ', 'MERGE_MATCH'
]


#
# Functions
#

def intersect_is_read_seq(wildcards, config):
    """
    Determine if merge requires input sequence.

    :param wildcards: Rule wildcards.
    :param config: Configuration.

    :return: `True` if sequences should be read.
    """

    config_def = svpoplib.svmerge.get_merge_def(wildcards.merge_def, config)

    if config_def is None:
        config_def = wildcards.merge_def

    return svpoplib.svmergeconfig.params.get_merge_config(config_def).read_seq


def intersect_get_input(pattern, wildcards, config):
    """
    Format a pattern according to input parameters. Used to find the correct files for intersect input.

    :param pattern: A list of file two input patterns [a, b].
    :param wildcards: Rule wildcards.
    :param config: Pipeline config.

    :return: `pattern` with wildcards filled in.
    """

    # Pattern must be a list of two elements (a and b)
    if type(pattern) != list or len(pattern) != 2:
        raise RuntimeError(f'Pattern must be a list of two strings')

    # Determine type
    if all([file_name.endswith('.bed.gz') for file_name in pattern]):
        file_type = 'bed'
    elif all([file_name.endswith('.fa.gz') for file_name in pattern]):
        file_type = 'fa'
    else:
        raise RuntimeError(f'Pattern does not contain a list of BED or FASTA files')

    # Stop if FASTA and no sequence
    if file_type == 'fa' and not intersect_is_read_seq(wildcards, config):
        return []

    # Dict of wildcards
    wc_dict = dict(wildcards)

    # Split filter
    if '_vs_' in wc_dict['filter']:
        tok = wc_dict['filter'].split('_vs_', 1)

        wc_dict['filter_a'] = tok[0]
        wc_dict['filter_b'] = tok[1]
    else:
        wc_dict['filter_a'] = wc_dict['filter']
        wc_dict['filter_b'] = wc_dict['filter']

    # Split svset
    if '_vs_' in wc_dict['svset']:
        svset_tok = wc_dict['svset'].split('_vs_', 1)

        wc_dict['svset_a'] = svset_tok[0]
        wc_dict['svset_b'] = svset_tok[1]
    else:
        wc_dict['svset_a'] = wc_dict['svset']
        wc_dict['svset_b'] = wc_dict['svset']

    # Format
    return [
        pattern_str.format(**wc_dict) for pattern_str in pattern
    ]


def run_intersect(bed_list, strategy, fa_list=None, threads=1):
    """
    Intersect two callsets and generate a table of matching and non-matching variants.

    :param bed_list: List of input BED files.
    :param strategy: Merge strategy.
    :param fa_list: List of sequence FASTA files. If sequences are not needed for the match, this should be `None`.
    :param threads: Number of threads to use for the match.

    :return: Table of matching and non-matching variants.
    """

    # Check arguments
    if len(bed_list) != 2:
        raise RuntimeError(f'Intersect only supports 2 samples, found {len(bed_list)}')

    if fa_list is not None and len(fa_list) != len(bed_list):
        raise RuntimeError(f'FASTA list length does not match the BED list length: bed_list={len(bed_list)}, fa_list={len(fa_list)} (may be None for no sequence match)')

    # Merge
    df = svpoplib.svmerge.merge_variants(
        bed_list=bed_list,
        sample_names=['A', 'B'],
        strategy=strategy,
        threads=threads,
        fa_list=fa_list
    )

    support_col_list = [col for col in df.columns if col in MERGE_INFO_FIELD_LIST]

    # Subset columns
    df = df.loc[:, ['ID', 'MERGE_SAMPLES', 'MERGE_SRC', 'MERGE_VARIANTS'] + support_col_list]

    # Create an ID column for sample (empty string if the variant was not in that sample)
    df['MERGE_SAMPLES'] = df['MERGE_SAMPLES'].apply(lambda val:
        val + ',' if val == 'A' else (',' + val if val == 'B' else val)
    )

    df['MERGE_VARIANTS'] = df.apply(lambda row:
        row['MERGE_VARIANTS'] + ',' if row['MERGE_SAMPLES'] == 'A,' else (',' + row['MERGE_VARIANTS'] if row['MERGE_SAMPLES'] == ',B' else row['MERGE_VARIANTS']),
        axis=1
    )

    df['ID_A'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[0])
    df['ID_B'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[1])

    # Set support columns
    new_col_list = list()

    for col in support_col_list:
        new_col = col[len('MERGE_'):]
        new_col_list.append(new_col)

        split_list = df[col].apply(lambda val: val.split(',') if not pd.isnull(val) else '')

        df[new_col] = split_list.apply(lambda val: val[1] if len(val) > 1 else np.nan)

    # Subset
    df['SOURCE_SET'] = df['MERGE_SAMPLES']
    df = df.loc[:, ['ID_A', 'ID_B', 'SOURCE_SET'] + new_col_list]

    return df
