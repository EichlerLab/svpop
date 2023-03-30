"""
Functions for making tracks.
"""

# Packages
import pandas as pd
import numpy as np
import re

import subprocess

import os

import svpoplib.util

### Definitions ###

FAI_FILE_NAME = '/net/eichler/vol27/projects/autism_genome_assembly/nobackups/sv/reference/hg38.no_alt.fa.fai'

TRACKS_SVTYPE_COLOR_DICT = {
    'INS': '64,64,255',
    'DEL': '255,64,64',
    'INV': '96,255,96',
    'DUP': '255,64,255',
    'SNV': '0,0,0',
    'RGN': '0,0,0',
    'SUB': '0,0,0',
    'TD-DONOR': '72,72,72'  # Not used, but can be for INS donor sites
}


COL_DTYPE = {
    'POS': np.int32,
    'END': np.int32,
    'POS_THICK': np.int32,
    'END_THICK': np.int32,
    'SVLEN': np.int32,
    'DP': np.int32,
    'GQ': np.int32,
    'QA': np.int32,
    'CALLERSET_N': np.int32,
    'MERGE_AC': np.int32
}

COL_NA_FILL = {
    'POS': 0,
    'END': 0,
    'POS_THICK': 0,
    'END_THICK': 0,
    'SVLEN': -1,
    'DP': -1,
    'GQ': -1,
    'QA': -1,
    'CALLERSET_N': -1,
    'MERGE_AC': -1,
}

TYPE_DICT = {
    'string': str,
    'lstring': str,
    'uint': int,
    'int': int,
    'float': float,
}


def report_status(message, verbose):
    """
    Report a message if verbose is true.
    """
    if verbose:
        print('Status: ' + message)


def make_bb_track(df, df_fai, bed_file_name, as_file_name, track_name, track_description, field_table_file_name=None, verbose=True):

    # Trim large INS records that extend past chromosome end
    df['END'] = df['POS'] + df['SVLEN']
    df['END'] = df.apply(lambda row: np.min([row['END'], df_fai[row['#CHROM']]]), axis=1)

    # Add BED fields
    report_status('Formatting BED', verbose)

    df['SCORE'] = 1000
    df['STRAND'] = '.'
    df['POS_THICK'] = df['POS']
    df['END_THICK'] = df['END']

    df['COL'] = df.apply(lambda row: TRACKS_SVTYPE_COLOR_DICT[row['SVTYPE']], axis=1)

    # Arrange columns
    head_cols = ['#CHROM', 'POS', 'END', 'ID', 'SCORE', 'STRAND', 'POS_THICK', 'END_THICK', 'COL']
    tail_cols = [col for col in df.columns if col not in head_cols]

    df = df.loc[:, head_cols + tail_cols]

    df.sort_values(['#CHROM', 'POS', 'ID'], inplace=True)

    # Set types and default values
    for col in set(COL_NA_FILL.keys()) & set(df.columns):
        df[col] = df[col].fillna(COL_NA_FILL[col])

    for col in [col for col in COL_DTYPE.keys() if col in df.columns]:
        df[col] = df[col].astype(COL_DTYPE[col])

    ### Define AS columns (AutoSQL, needed to make a BigBed) ###

    if field_table_file_name is None:
        field_table_file_name = os.path.join(svpoplib.util.get_install_dir(), 'files/tracks/ucsc_track_fields.tsv')

    df_as = pd.read_csv(
        field_table_file_name,
        sep='\t', header=0,
        dtype={'DEFAULT': object},
        na_values=[''], keep_default_na=False
    )

    df_as.set_index('FIELD', inplace=True, drop=False)

    missing_list = [col for col in tail_cols if col not in set(df_as['FIELD'])]

    if missing_list:
        raise RuntimeError('Missing AS definitions for columns: {}'.format(', '.join(missing_list)))

    # Reformat columns
    for col in df.columns:

        if 'DEFAULT' in df_as.columns and not pd.isnull(df_as.loc[col, 'DEFAULT']):
            default_val = df_as.loc[col, 'DEFAULT']
        else:
            default_val = '.'

        try:
            df[col] = format_column(df[col], TYPE_DICT.get(df_as.loc[col, 'TYPE'], str), default_val=default_val)
        except Exception as ex:
            raise RuntimeError('Error formatting {} as {}: {}'.format(col, df_as.loc[col, 'TYPE'], ex))

    # Write
    report_status('Writing', verbose)

    with open(as_file_name, 'w') as out_file:
        # Heading
        out_file.write('table {}\n"{}"\n(\n'.format(track_name, track_description))

        # Column definitions
        for col in head_cols + tail_cols:
            out_file.write('{TYPE} {NAME}; "{DESC}"\n'.format(**df_as.loc[col]))

        # Closing
        out_file.write(')\n')

    # Set column defaults based on AS field type

    ### Make BED ###

    df.to_csv(bed_file_name, sep='\t', na_rep='.', index=False)

    report_status('Done', verbose)


def format_column(row, type, default_val='.'):

    non_missing = (~ pd.isnull(row)) & (row != '.')

    # Process non-nan values (to string)
    #row_nonan = row.loc[non_missing].astype(type).astype(str)
    row_nonan = row.loc[non_missing].apply(lambda val_list:
        ','.join(
            [str(type(val)) for val in (val_list.split(',') if isinstance(val_list, str) else [val_list])]
        )
    )

    # Process nan as empty strings
    nan_vals = row.loc[~ non_missing]

    row_nan = pd.Series([default_val] * nan_vals.shape[0], index=nan_vals.index)

    # Concat, order, return
    return pd.concat([row_nonan, row_nan]).loc[row.index]

