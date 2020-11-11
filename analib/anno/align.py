"""
Routines for alignments to annotations.
"""

import numpy as np
import pandas as pd
import pysam


def get_align_path(alignset, sample, config):
    """
    Get path to alignment file (BAM/CRAM) for an alignment set (e.g. "illumina" or "clr").

    :param alignset: Alignment set name.
    :param sample: Aligned sample.

    :return: Path to alignment file.
    """

    align_path = config.get('align', {}).get(alignset, None)

    if align_path is None:
        raise RuntimeError('Cannot find alignment file for alignment set {}'.format('alignset'))

    if '{sample}' not in align_path:
        raise RuntimeError('Wildcard "{{sample}}" not in alignment path for set {}: {}'.format(alignset, align_path))

    return align_path.format(sample=sample)

def get_depth(df, alignset, sample, config):
    """
    Get alignment depth over the variant.

    :param df: Dataframe of variant calls.
    :param alignset: Alignment set (in "align" config section).
    :param sample: Sample name.
    :param config: Configuration dict.

    :return: Dataframe with ID and DEPTH columns.
    """

    # Get path to alignment file
    align_path = get_align_path(alignset, sample, config)

    # Empty table if no variants
    if df is None or df.shape[0] == 0:
        return pd.DataFrame([], columns=['ID', 'DEPTH'])

    df = df.set_index('ID')

    # Get depth
    with pysam.AlignmentFile(align_path) as in_file:
        df_depth = df.apply(lambda row: len(list(in_file.fetch(row['#CHROM'], row['POS'], row['END']))), axis=1)

    df_depth.index.name = 'ID'
    df_depth.name = 'DEPTH'

    return pd.DataFrame(df_depth).reset_index()
