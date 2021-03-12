"""
Reference genome utilities.
"""

import hashlib
import numpy as np
import os
import pandas as pd
import pysam
import re

from Bio import SeqIO

import svpoplib.util
import svpoplib.seq


def grc_to_hg_chrom(chroms, grc_build, rev=False, missing=False, grc_tab='files/ref/{grc_build}/grc_report.tab'):
    """
    Translate chromosome names between GRC (GRCh) and UCSC (hg) chromosome names. If not `rev` (default), then translate
    from GRC to UCSC. If `rev`, then UCSC to GRC.

    :param chroms: List of Series of chromosome names.
    :param grc_build: GRC build name (e.g. "GRCh38"). This is parsed into `grc_tab`.
    :param rev: UCSC to GRC if set, GRC to UCSC otherwise.
    :param missing: Return `np.nan` for values with no mapping if `True`.
    :param grc_tab: Pattern to locate the table containing "UCSC-style-name" and "Sequence-Name". String contains format
        element "{grc_build}" for the GRC reference name (e.g. "GRCh38").

    :return: List or Series of translated chromosome names.
    """

    # Get column names
    if not rev:
        from_col = 'Sequence-Name'
        to_col = 'UCSC-style-name'
    else:
        from_col = 'UCSC-style-name'
        to_col = 'Sequence-Name'

    # Read table
    grc_tab_path = os.path.join(
        svpoplib.util.get_install_dir(),
        grc_tab.format(**{'grc_build': grc_build})
    )

    df_grc_report = pd.read_csv(grc_tab_path, sep='\t', header=0)

    # Get primary contig names. Also get chromosomes keyed on GenBank and RefSeq accessions.
    grc_table = df_grc_report.loc[:, [from_col, to_col]].set_index(from_col).squeeze().append(
        df_grc_report.loc[:, ['GenBank-Accn', to_col]].set_index('GenBank-Accn').squeeze()
    ).append(
        df_grc_report.loc[:, ['RefSeq-Accn', to_col]].set_index('RefSeq-Accn').squeeze()
    )

    # Get a series of chromosome names (if not already)
    if not isinstance(chroms, pd.Series):
        chroms_as_list = True
    else:
        chroms_as_list = False
        chroms = pd.Series(chroms)


    # Translate, set NaN for missing values
    key_set = set(grc_table.index)

    chroms_trans = chroms.apply(str).apply(
        lambda val: grc_table[val] if val in key_set else np.nan
    )

    # Check for missing values
    if not missing:
        missing_list = list(chroms.loc[pd.isnull(chroms_trans)])

        if len(missing_list) > 0:
            missing_set = set(missing_list)

            missing_str = ', '.join(sorted(missing_set)[:3])

            if len(missing_set) > 3:
                missing_str += ', ...'

            raise RuntimeError(
                'Cannot translate chromosome for {0} record(s): Chromosome(s) = {1}'.format(
                    len(missing_list),
                    missing_str
                )
            )

    # Return translated values
    return chroms_trans


def hg_chr_scaffold(chroms):
    """
    Determine if each chromosome name is a proper chromosome (primary assembly, not unplaced or unlocalized).

    :param chroms: List of Series of chromosome names.

    :return: A list or Series of boolean values (`True` if proper chromosome).
    """

    if isinstance(chroms, pd.Series):
        return chroms.apply(lambda val: bool(re.match('^chr[1-9XY][0-9]?$', val)))
    else:
        return [bool(re.match('^chr[1-9XY][0-9]?$', val)) for val in chroms]


def get_df_fai(fai_file_name, usecols=('CHROM', 'LEN'), index_col='CHROM', squeeze=True):
    """
    Read an FAI File name. By default, return a Series of chromosome lengths keyed by the chromosome (or contig) name.

    :param fai_file_name: File to read.
    :param usecols: Use these columns. Must include index column if set.
    :param index_col: Index column name.
    :param squeeze: Squeeze to Series if only one column is selected (not counting the index column).

    :return: A DataFrame or Series of the FAI file.
    """

    return pd.read_csv(
        fai_file_name,
        sep='\t',
        names=['CHROM', 'LEN', 'POS', 'LINE_BP', 'LINE_BYTES'],
        usecols=usecols,
        index_col=index_col,
        squeeze=squeeze,
        dtype={'CHROM': str, 'LEN': np.int, 'POS': np.int, 'LINE_BP': np.int, 'LINE_BYTES': np.int}
    )


def get_ref_region(df, ref_fa):
    """
    Get a reference region.

    :param df: Dataframe with "#CHROM", "POS", and "END" fields.

    :return: A Series of sequences corresponding to the rows in df.
    """

    with pysam.FastaFile(ref_fa) as fa_file:
        return df.apply(lambda row: fa_file.fetch(row['#CHROM'], row['POS'], row['END']), axis=1)

def get_ref_info(ref_fa):
    """
    Get a table of reference information from a FASTA file. FASTA must have a ".fai".

    :param ref_fa: Reference FASTA.

    :return: Pandas DataFrame with index "CHROM" and columns "LEN", "POS", "LINE_BP", "LINE_BYTES" derived from the FAI
        file. Also has "MD5" (Sequence MD5) and "ORDER" (order of records in the FASTA).
    """

    # Read FAI
    df_fai = svpoplib.ref.get_df_fai(ref_fa + '.fai', usecols=None)

    record_list = list()
    record_count = 0

    # Read sequneces and get MD5
    with svpoplib.seq.PlainOrGzReader(ref_fa, 'rt') as in_file:
        for record in SeqIO.parse(in_file, 'fasta'):
            record_list.append(pd.Series(
                [
                    record.id,
                    len(record),
                    hashlib.md5(str(record.seq).upper().encode()).hexdigest(),
                    record_count
                ],
                index=['CHROM', 'LEN', 'MD5', 'ORDER']
            ))

            record_count += 1

    df = pd.concat(record_list, axis=1).T
    df.set_index('CHROM', inplace=True)

    # Check completeness and sizes against the FAI
    if set(df.index) != set(df_fai.index):
        missing_df = set(df_fai.index) - set(df.index)
        missing_fai = set(df.index) - set(df_fai.index)

        raise RuntimeError('Reference and FAI mismatch: {} missing in FASTA ({}{}), {} missing in FAI ({}{})'.format(
            len(missing_df), ', '.join(sorted(missing_df)[:3]) if missing_df else 'OK',
            '...' if len(missing_df) > 3 else '',
            len(missing_fai), ', '.join(sorted(missing_fai)[:3]) if missing_fai else 'OK',
            '...' if len(missing_fai) > 3 else '',
        ))

    if np.any(df['LEN'] != df_fai['LEN']):
        df_mismatch = df.loc[df['LEN'] != df_fai['LEN']]

        raise RuntimeError('Reference record lengths do not match for {} records: {}{}'.format(
            df_mismatch.shape[0],
            ', '.join(
                ['{} (FASTA={:,d}, FAI={:,d})'.format(chrom, df_mismatch.loc[chrom, 'LEN'], df_fai[chrom]['LEN']) for
                 chrom in df_mismatch.index[:3]]),
            '...' if df_mismatch.shape[0] > 3 else ''
        ))

    # Add fields to df_fai
    df_fai['MD5'] = df['MD5']
    df_fai['ORDER'] = df['ORDER']

    return df_fai
