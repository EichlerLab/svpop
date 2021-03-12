"""
Variant processing and comparison functions.
"""

import intervaltree
import multiprocessing
import numpy as np
import os
import pandas as pd
import re


def reciprocal_overlap(begin_a, end_a, begin_b, end_b):
    """
    Get reciprocal overlap of two intervals. Intervals are expected to be half-open coordinates (length is end - start).

    :param begin_a: Begin of interval a.
    :param end_a: End of interval a.
    :param begin_b: Begin of interval b.
    :param end_b: End of interval b.

    :return: A value between 0 and 1 if intervals a and b overlap and some negative number if they do not.
    """

    overlap = min(end_a, end_b) - max(begin_a, begin_b)

    if overlap < 0:
        return 0.0

    return min([
        overlap / (end_a - begin_a),
        overlap / (end_b - begin_b)
    ])


def var_nearest(df_a, df_b, ref_alt=False, verbose=False):
    """
    For each variant in `df_a`, get the nearest variant in `df_b`. All `df_a` variants are in the output except those
    where `df_b` has no variant call on the same chromosome.

    :param df_a: Variants to match.
    :param df_b: Variants to match against.
    :param ref_alt: Find distance to nearest variant where "REF" and "ALT" columns match (for SNVs of the same type).
    :param verbose: Print status information.

    :return: A dataframe with columns "ID_A", "ID_B", and "DISTANCE". If a variant from `df_b` is downstream, the
        distance is positive.
    """

    # Check
    if ref_alt:
        if 'REF' not in df_a.columns or 'ALT' not in df_a.columns:
            raise RuntimeError('Missing required column(s) in dataframe A for ref-alt comparisons: REF, ALT')

        if 'REF' not in df_b.columns or 'ALT' not in df_b.columns:
            raise RuntimeError('Missing required column(s) in dataframe B for ref-alt comparisons: REF, ALT')

    # Process each chromosome
    match_list = list()

    for chrom in sorted(set(df_a['#CHROM'])):
        if verbose:
            print('Chrom: ' + chrom)

        # Subset by chromosome
        df_b_chrom = df_b.loc[df_b['#CHROM'] == chrom]

        if df_b_chrom.shape[0] == 0:
            continue

        df_a_chrom = df_a.loc[df_a['#CHROM'] == chrom]

        # Get a set of REF-ALT tuples from df_a (if ref_alt)
        if ref_alt:
            ref_alt_set = set(df_a_chrom.loc[:, ['REF', 'ALT']].apply(tuple, axis=1))
        else:
            ref_alt_set = {(None, None)}

        # Process each subset for this chromosome
        for ref, alt in ref_alt_set:

            # Separate on REF/ALT (if ref_alt is True)
            if ref is not None and alt is not None:
                df_a_sub = df_a_chrom.loc[(df_a_chrom['REF'] == ref) & (df_a_chrom['ALT'] == alt)]
                df_b_sub = df_b_chrom.loc[(df_b_chrom['REF'] == ref) & (df_b_chrom['ALT'] == alt)]
            else:
                df_a_sub = df_a_chrom
                df_b_sub = df_b_chrom

            if df_a_sub.shape[0] == 0 or df_b_sub.shape[0] == 0:
                continue

            # Get arrays from b for comparisons
            pos_array = np.array(df_b_sub['POS'])
            end_array = np.array(df_b_sub['END'])
            id_array = np.array(df_b_sub['ID'])

            # Process each record on chromosome
            for index, row in df_a_sub.iterrows():

                min_pos_index = np.argmin(np.abs(pos_array - row['POS']))
                min_end_index = np.argmin(np.abs(end_array - row['END']))

                min_pos = row['POS'] - pos_array[min_pos_index]
                min_end = row['END'] - end_array[min_end_index]

                # Make intersect record
                if np.abs(min_pos) < np.abs(min_end):
                    match_list.append(pd.Series(
                        [
                            row['ID'],
                            id_array[min_pos_index],
                            min_pos,
                        ],
                        index=['ID_A', 'ID_B', 'DISTANCE']
                    ))

                else:
                    match_list.append(pd.Series(
                        [
                            row['ID'],
                            id_array[min_end_index],
                            min_end,
                        ],
                        index=['ID_A', 'ID_B', 'DISTANCE']
                    ))

    # Return merged dataframe
    return pd.concat(match_list, axis=1).T


def nr_interval_merge(df_chr, overlap=0.5):
    """
    Reduce a dataframe to non-redundant intervals based on reciprocal overlap. All records in the dataframe must be
    on the same chromosome.

    :param df_chr: DataFrame of one chromosome.
    :param overlap: Reciprocal overlap (0, 1].

    :return: Dataframe subset using the first record in a unique interval.
    """

    index_list = list()  # Dataframe indices to return

    interval_tree = intervaltree.IntervalTree()  # Tree of intervals

    # Iterate rows
    for index, row in df_chr.iterrows():
        ri_match = False

        # Find matches
        for interval in interval_tree[row['POS']:row['END']]:
            if reciprocal_overlap(row['POS'], row['END'], interval.begin, interval.end) >= 0.50:
                ri_match = True
                break

        # Append to non-redundant records if no match
        if not ri_match:
            index_list.append(index)

        # All records are added to the tree
        interval_tree[row['POS']:row['END']] = True

    return df_chr.loc[index_list]


def order_variant_columns(
        df, head_cols=('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN'), tail_cols=None, allow_missing=False
):
    """
    Rearrange columns with a set list first (in defined order of `head_cols`) and leave the remaining columns
    in the order they were found.

    :param df: Data frame.
    :param head_cols: Columns to move to the first columns. Set variant BED order by default.
    :param tail_cols: Columns to move to the end. May be set to `None`.
    :param flex_col: All

    :return: Data frame with rearranged columns.
    """

    # Check head columns
    if head_cols is not None:
        head_cols = list(head_cols)
    else:
        head_cols = list()

    if not allow_missing:
        for col in head_cols:
            if col not in df.columns:
                raise RuntimeError('Missing head column in variant file: {}'.format(col))
    else:
        head_cols = [col for col in head_cols if col in df.columns]

    # Check tail columns
    if tail_cols is not None:
        tail_cols = list(tail_cols)
    else:
        tail_cols = list()

    if not allow_missing:
        for col in tail_cols:
            if col not in df.columns:
                raise RuntimeError('Missing tail column in variant file: {}'.format(col))
    else:
        tail_cols = [col for col in tail_cols if col in df.columns]

    # Give precedence for head columns if a col is head and tail
    tail_cols = [col for col in tail_cols if col not in head_cols]

    # Check for empty column sets
    if len(head_cols) == 0 and len(tail_cols) == 0:
        raise RuntimeError('No head or tail columns to sort (after filtering for missing columns if allow_missing=True)')

    # Define middle columns
    head_tail_set = set(head_cols).union(set(tail_cols))

    mid_cols = [col for col in df.columns if col not in head_tail_set]

    # Arrange with head columns first. Leave remaining columns in order
    return df.loc[:, head_cols + mid_cols + tail_cols]


def get_variant_id(df):
    """
    Get variant IDs using '#CHROM', 'POS', 'SVTYPE', and 'SVLEN' columns.

    :param df: Dataframe.

    :return: A Series of variant IDs for `df`.
    """

    # Set ID
    return df.apply(get_variant_id_from_row, axis=1)


def get_variant_id_from_row(row):
    """
    Get variant ID for one row.

    :param row: Variant row.

    :return: Variant ID.
    """

    if row['SVTYPE'] != 'SNV':
        return '{}-{}-{}-{}'.format(
            row['#CHROM'], row['POS'] + 1, row['SVTYPE'], row['SVLEN']
        )

    else:
        return '{}-{}-{}-{}-{}'.format(
            row['#CHROM'], row['POS'] + 1, row['SVTYPE'], row['REF'].upper(), row['ALT'].upper()
        )

def vcf_fields_to_seq(row, pos_row='POS', ref_row='REF', alt_row='ALT'):
    """
    Get call for one VCF record and one sample.

    Example command-line to generate input table from a VCF:

    bcftools query -H -f"%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT][\t%GQ][\t%DP][\t%AD]\n" INPUT.vcf.gz | gzip > OUTPUT.vcf.tab.gz

    :param row: Table row.
    :param pos_row: Row name for the variant position (POS).
    :param ref_row: Row name for the reference sequence (REF).
    :param alt_row: Row name for the alternate sequence (ALT).

    :return: Standard BED format of variant calls with "POS", "END", "VARTYPE", "SVTYPE", "SVLEN", "SEQ", "REF".
    """

    pos = row[pos_row]
    ref = row[ref_row].upper()
    alt = row[alt_row].upper()

    # Handle symbolic SV variants
    if alt[0] == '<' and alt[-1] == '>':

        # Get type
        svtype = alt[1:-1].split(':', 1)[0]

        if svtype not in {'INS', 'DEL', 'INV', 'DUP', 'CNV'}:
            raise RuntimeError('Unrecognized symbolic variant type: {}: Row {}'.format(svtype, row.name))

        # Set pos
        pos = row['POS']

        # Get length
        if 'SVLEN' not in row or pd.isnull(row['SVLEN']) or row['SVLEN'] == 0:

            if 'END' not in row:
                raise RuntimeError('Missing or 0-length SVLEN and no END for symbolic SV: Row {}'.format(row.name))

            if svtype == 'INS':
                raise RuntimeError('Cannot calculate SVLEN for insertion: Row {}'.format(row.name))

            svlen = row['END'] - row['POS']

        else:
            svlen = row['SVLEN']


        # Set variant type
        if svlen < 50:
            raise RuntimeError('Symbolic ALT for non SV: Row {}'.format(row.name))

        vartype = 'SV'

        # Set end
        if svtype == 'INS':
            end = pos + 1
        else:
            end = pos + svlen

        # Sequence
        seq = row['SEQ'] if 'SEQ' in row else np.nan

    else:

        min_len = min(len(ref), len(alt))

        trim_left = 0

        # Trim left
        while min_len and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            trim_left += 1
            min_len -= 1

        # Trim right
        while min_len and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
            min_len -= 1

        # Check variant type
        if ref == '' and alt != '':
            svtype = 'INS'
            seq = alt
            svlen = len(seq)
            vartype = 'INDEL' if svlen < 50 else 'SV'
            pos = pos + trim_left - 1
            end = pos + 1

        elif ref != '' and alt == '':
            svtype = 'DEL'
            seq = ref
            svlen = len(seq)
            vartype = 'INDEL' if svlen < 50 else 'SV'
            pos = pos + trim_left - 1
            end = pos + svlen
            ref = ''

        elif len(ref) == 1 and len(alt) == 1:
            vartype = 'SNV'
            svtype = 'SNV'
            seq = alt
            svlen = 1
            pos = pos + trim_left - 1
            end = pos + 1

        else:
            vartype = 'SUB'
            svtype = 'SUB'
            seq = alt
            svlen = len(seq)
            pos = pos + trim_left - 1
            end = pos + svlen

    # Return with AC
    if 'GT' in row.index:
        ac = gt_to_ac(row['GT'], no_call=-1, no_call_strict=False)

        return pd.Series(
            [pos, end, vartype, svtype, svlen, ac, seq, ref],
            index=['POS', 'END', 'VARTYPE', 'SVTYPE', 'SVLEN', 'AC', 'SEQ', 'REF']
        )

    # Return without AC
    return pd.Series(
        [pos, end, vartype, svtype, svlen, seq, ref],
        index=['POS', 'END', 'VARTYPE', 'SVTYPE', 'SVLEN', 'SEQ', 'REF']
    )


def gt_to_ac(gt, no_call=-1, no_call_strict=False):
    """
    Convert a genotype string to an allele count.

    :param gt: Genotype string. Examples are "0/1", "0|1", "./.", ".|1".
    :param no_call: If all alleles are no-call ("."), then return this value.
    :param no_call_strict: If any alleles are no-call, return `no_call`.

    :return: Genotype count.
    """
    ac = 0

    gt_list = re.split('[/|]', gt)

    # Handle all no-call
    if '.' in gt_list:
        if set(gt_list) == {'.'} or no_call_strict:
            return no_call

    # Get AC
    for gt_allele in gt_list:
        if gt_allele != '.' and int(gt_allele) > 0:
            ac += 1

    return ac


def get_filter_bed(filter_name, ucsc_ref_name, config, svpop_dir):
    """
    Get a BED file defining a filter. Searches config['filter'] for the filter name (key) and path (value). If not
    found, the default is "files/filter/{ucsc_ref_name}/{filter}.bed" within SV-Pop pipeline directory.

    If wildcard "ucsc_ref_name" is in the filter path, then "ucsc_ref_name" is parsed into it.

    :param filter_name: Fliter name.
    :param ucsc_ref_name: Name of the UCSC reference (e.g. "hg38").
    :param config: SV-Pop config.
    :param svpop_dir: SV-Pop pipeline directory.

    :return: Filter path.
    """

    # Get path
    filter_path = config.get('filter', dict()).get(filter_name, None)

    if filter_path is None:
        filter_path = os.path.join(
            svpop_dir, 'files/filter/{ucsc_ref_name}/{filter}.bed'.format(
                ucsc_ref_name=ucsc_ref_name,
                filter=filter_name
            )
        )

    elif '{ucsc_ref_name}' in filter_path:
        filter_path = filter_path.format(
            ucsc_ref_name=ucsc_ref_name,
        )

    # Check path
    if not os.path.isfile(os.path.join(svpop_dir, filter_path)):
        raise RuntimeError('Cannot find filter {}: {}'.format(filter_name, filter_path))

    # Return
    return filter_path


def qual_to_filter(row, min_qv=30.0):
    """
    Use VCF "QUAL" field to fill "FILTER".

    :return: "PASS" if "QUAL" is numeric and greater than or equal to `min_qv`, "FAIL" if "PASS" if "QUAL" is numeric
        and less than `min_qv`, and "." if "QUAL" is not numeric.
    """

    if row['FILTER'] == '.':
        try:
            return 'PASS' if float(row['QUAL']) >= min_qv else 'FAIL'
        except ValueError:
            return '.'
    else:
        return row['FILTER']
