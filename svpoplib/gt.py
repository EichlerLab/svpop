"""
Genotyping routines.
"""

import numpy as np
import pandas as pd


def fst_wc(row):
    """
    Weir and Cockerham's FST (Code from Benson).

    :param row: Row of genotype data with multi-index where the first level is "POP_A" or "POP_B" and the second level
        is "N" (number of samples with a genotype call) or "AF" (allele frequency of genotyped calls).
    """

    # Get population sizes
    list_popsizes = (
        row[('POP_A', 'N')] * 2,
        row[('POP_B', 'N')] * 2
    )

    # Get allele frequencies
    list_freqs = (
        row[('POP_A', 'AF')],
        row[('POP_B', 'AF')]
    )

    # No calculation if AF is NA (no alleles, all no-calls)
    if pd.isnull(list_freqs[0]) or pd.isnull(list_freqs[1]):
        return np.nan

    # No calculation if AF is 1 for both populations
    if all(((1 - val) < np.float64(0.0000001) for val in list_freqs)):
        return np.float64(0.0)

    # calcuate Ht = 2 * p_bar * q_bar; p refer to the sum of ingroup allele freqs. q is for the same quality for outgroup.
    # calculate Hs = (1/n) * sum of 2 * pi * qi, for i  in range(1,n+1) and n =numPop
    list_popsizes = [np.float64(x) for x in list_popsizes]
    list_freqs = [np.float64(x) for x in list_freqs]

    # Get number of unique alleles (r) and total number of alleles (n)
    r = float(len(list_popsizes))
    n = float(sum(list_popsizes))

    p_bar = 0.0

    for i, v in enumerate(list_freqs):
        p_bar += (list_popsizes[i] / n) * list_freqs[i]

    n_bar = sum(list_popsizes) / r

    S_sq_part1 = 1.0 / ((r - 1) * n_bar)
    S_sq_part2 = 0.0

    sum_n_i = 0.0
    sum_n_i_sq = 0.0

    for i, v in enumerate(list_freqs):
        S_sq_part2 += list_popsizes[i] * ((list_freqs[i] - p_bar)**2)
        sum_n_i += list_popsizes[i]
        sum_n_i_sq += (list_popsizes[i])**2

    S_sq = S_sq_part1 * S_sq_part2

    T1 = S_sq - (1.0 / (2*n_bar - 1)) * (p_bar * (1-p_bar) - (r-1)/r * S_sq)

    nc_part1 = 1.0 / (r - 1)
    nc_part2 = (sum_n_i - (sum_n_i_sq / sum_n_i))

    nc = nc_part1 * nc_part2

    T2_part1 = (2 * nc - 1) / (2 * n_bar - 1)
    T2_part2 = p_bar * (1- p_bar)
    T2_part3 = (1 + (2*(r-1)*(n_bar - nc)) / (2 * n_bar - 1))
    T2_part4 = S_sq / r

    T2 = T2_part1 * T2_part2 + T2_part3 * T2_part4

    locus_fst = T1 / T2

    # Return value
    locus_fst = np.float64(locus_fst)

    if locus_fst >= 0.0:
        return locus_fst
    else:
        return np.float64(0.0)


def get_gt_summary_table(ins_bed, del_bed):
    """
    Get summary statistics for variants by GT call.

    Output table:
    * rows: HOM-REF, HOM-ALT, HET, NO-CALL
    * cols: INS_N, INS_SUM, INS_MED, DEL_N, DEL_SUM, DEL_MED, ALL_N, ALL_SUM, ALL_MED

    Input BED files must have a 'SVLEN' field.

    :param ins_bed: Insertions BED file.
    :param del_bed: Deletions BED file.

    :return: A dataframe of summary statistics for each genotype call.
    """

    # Read
    df_ins = pd.read_csv(ins_bed, sep='\t', header=0)
    df_del = pd.read_csv(del_bed, sep='\t', header=0)

    # ABS SVLEN
    df_del['SVLEN'] = df_del['SVLEN'].apply(abs)

    # Create data frame
    df = pd.DataFrame(data=0,
                      columns=(
                          'INS_N', 'INS_SUM', 'INS_MED',
                          'DEL_N', 'DEL_SUM', 'DEL_MED',
                          'ALL_N', 'ALL_SUM', 'ALL_MED'
                      ),
                      index=('HOM-REF', 'HOM-ALT', 'HET', 'NO-CALL'), dtype=np.int32
                      )


    # Summarize insertions
    group = df_ins.groupby('GT')

    df['INS_N'] = group.count()['SVLEN']
    df['INS_SUM'] = group.sum()['SVLEN']
    df['INS_MED'] = group.median()['SVLEN'].astype(np.int32)

    # Summarize deletions
    group = df_del.groupby('GT')

    df['DEL_N'] = group.count()['SVLEN']
    df['DEL_SUM'] = group.sum()['SVLEN']
    df['DEL_MED'] = group.median()['SVLEN'].astype(np.int32)

    # Summarize merged ins & del
    df_all = pd.concat([df_ins, df_del], axis=0)

    group = df_all.groupby('GT')

    df['ALL_N'] = group.count()['SVLEN']
    df['ALL_SUM'] = group.sum()['SVLEN']
    df['ALL_MED'] = group.median()['SVLEN'].astype(np.int32)

    # Return
    return df


def get_gt_db_summary_table(summary_bed_file):
    """
    Generate summary statistics for a variant database genotyping is run on.

    Table:
    * rows: One row per samples
    * cols: INS_N, INS_SUM, INS_MED, DEL_N, DEL_SUM, DEL_MED, ALL_N, ALL_SUM, ALL_MED

    :param summary_bed_file: Summary BED file where the samples name is prepended to the CONTIG field (with a vertical
        bar). Must contain SVLEN and SVTYPE fields.

    :return: A dataframe of summary statistics.
    """

    # Read table
    df = pd.read_csv(summary_bed_file, sep='\t', header=0)

    # Get samples name variant was taken from (prepended to the contig)
    df['SAMPLE'] = df['CONTIG'].apply(lambda contig_name: contig_name.split('|', 2)[0])

    # Translate SVLEN
    df['SVLEN'] = df['SVLEN'].apply(abs)

    # Rearrange and drop columns
    df = df.ix[:, ('SAMPLE', 'SVTYPE', 'SVLEN')]

    # Split
    df_ins = df.ix[df['SVTYPE'] == 'INS', ('SAMPLE', 'SVLEN')]
    df_del = df.ix[df['SVTYPE'] == 'DEL', ('SAMPLE', 'SVLEN')]

    # Group by source samples
    group_ins = df_ins.groupby('SAMPLE')
    group_del = df_del.groupby('SAMPLE')
    group_all = df.groupby('SAMPLE')

    return pd.concat(
        [
            group_ins.count()['SVLEN'],
            group_ins.sum()['SVLEN'],
            group_ins.median()['SVLEN'],
            group_del.count()['SVLEN'],
            group_del.sum()['SVLEN'],
            group_del.median()['SVLEN'],
            group_all.count()['SVLEN'],
            group_all.sum()['SVLEN'],
            group_all.median()['SVLEN']
        ],
        keys=('INS_N', 'INS_SUM', 'INS_MED', 'DEL_N', 'DEL_SUM', 'DEL_MED', 'ALL_N', 'ALL_SUM', 'ALL_MED'),
        axis=1
    )
