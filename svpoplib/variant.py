"""
Variant processing and comparison functions.
"""

import collections
import intervaltree
import multiprocessing
import numpy as np
import os
import pandas as pd
import re

import svpoplib


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
        df, head_cols=None, tail_cols=None, allow_missing=False, subset=False
):
    """
    Rearrange columns with a set list first (in defined order of `head_cols`) and leave the remaining columns
    in the order they were found.

    :param df: Data frame.
    :param head_cols: Columns to move to the first columns. If None, defaults to ['#CHROM', 'POS', 'END', 'ID',
        'SVTYPE', 'SVLEN'] and includes 'REF' and 'ALT' if they are in the df.
    :param tail_cols: Columns to move to the end. May be set to `None`.
    :param allow_missing: Do not throw an error if the dataframe is missing one or more columns.
    :param subset: If True, subset to defined columns and drop all others.

    :return: Data frame with rearranged columns.
    """

    # Check head columns
    if head_cols is not None:
        head_cols = list(head_cols)
    else:
        head_cols = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN']

        if 'REF' in df.columns:
            head_cols += ['REF']

        if 'ALT' in df.columns:
            head_cols += ['ALT']

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

    if subset:
        mid_cols = []
    else:
        mid_cols = [col for col in df.columns if col not in head_tail_set]

    # Arrange with head columns first. Leave remaining columns in order
    return df.loc[:, head_cols + mid_cols + tail_cols]


def get_variant_id(df, apply_version=True, existing_id_set=None):
    """
    Get variant IDs using '#CHROM', 'POS', 'SVTYPE', and 'SVLEN' columns.

    :param df: Dataframe.
    :param apply_version: Version ID (add "." and a number for duplicated IDs). SV-Pop does not allow duplicate IDs, so
        this should be explicitly turned on unless duplicates are checked and handled explicitly. If there are no
        duplicate IDs before versioning, this option has no effect on the output ("." in only added if necessary).
    :param existing_id_set: A set of existing variant IDs that must also be avoided. If any IDs in id_col match these
        IDs, they are altered as if the variant intersects another ID in id_col.

    :return: A Series of variant IDs for `df`.
    """

    id_col = df.apply(get_variant_id_from_row, axis=1)

    if apply_version:
        id_col = version_id(id_col, existing_id_set=existing_id_set)

    return id_col


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
    ref = row[ref_row].upper().strip()
    alt = row[alt_row].upper().strip()

    # This function does not handle multiple alleles or missing ALTs (one variant per record)
    if ',' in alt:
        raise RuntimeError(f'Multiple alleles in ALT, separate before calling vcf_fields_to_seq(): "{alt}"')

    if alt == '.':
        raise RuntimeError('Missing ALT in record')

    # Fix common malformed VCFs (Sniffles2)
    if ref in {'N', 'n'} and re.match('^(?![Nn])[ACGTNacgtn]+(?<![Nn])$', alt) is not None:
        alt = ref + alt

    if alt in {'N', 'n'} and re.match('^(?![Nn])[ACGTNacgtn]+(?<![Nn])$', ref) is not None:
        ref = alt + ref

    # Fix common malformed VCFs (SVIM-asm)
    if ref == '.':
        ref = 'N'
        alt = 'N' + alt

    if alt == '.':
        alt = 'N'
        ref = 'N' + alt

    # Handle symbolic SV variants
    if alt[0] == '<' and alt[-1] == '>':

        # Get type
        svtype = alt[1:-1].split(':', 1)[0]

        if svtype not in {'INS', 'DEL', 'INV', 'DUP', 'CNV'}:
            raise RuntimeError('Unrecognized symbolic variant type: {}: Row {}'.format(svtype, row.name))

        # Get length
        svlen = None

        if 'SVLEN' in row:
            try:
                svlen = abs(int(row['SVLEN']))
            except:
                svlen = None
        else:
            svlen = None

        if svlen is None:
            if 'END' not in row:
                raise RuntimeError('Missing or 0-length SVLEN and no END for symbolic SV: Row {}'.format(row.name))

            try:
                svlen = abs(int(row['END'])) - pos
            except:
                raise RuntimeError('Variant has no SVLEN and END is not an integer: {}: Row {}'.format(row['END'], row.name))

        # Set variant type
        vartype = 'INDEL' if svlen < 50 else 'SV'

        # Set end
        if svtype == 'INS':
            end = pos + 1
        else:
            end = pos + svlen

        # Sequence
        seq = row['SEQ'] if 'SEQ' in row else np.nan

    elif alt == '.':
        vartype = 'NONE'
        svtype = 'NONE'
        seq = np.nan
        ref = np.nan
        end = pos + 1
        svlen = 0

    elif '[' in alt or ']' in alt or '.' in alt:
        vartype = 'BND'
        svtype = 'BND'
        seq = np.nan
        ref = np.nan
        end = pos + 1
        svlen = 0

    elif re.match('^[a-zA-Z]+$', alt) and re.match('^[a-zA-Z]+$', ref):

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

    else:
        raise RuntimeError(f'Unknown variant type: REF="{ref}", ALT="{alt}')

    # Return with AC
    if 'GT' in row.index:
        ac = gt_to_ac(row['GT'], no_call=-1, no_call_strict=False)
    else:
        ac = np.nan

    return pd.Series(
        [pos, end, vartype, svtype, svlen, ac, seq, ref],
        index=['POS', 'END', 'VARTYPE', 'SVTYPE', 'SVLEN', 'AC', 'SEQ', 'REF']
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

    gt_list = re.split('[/|]', str(gt))

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


def left_homology(pos_tig, seq_tig, seq_sv):
    """
    Duplicated from from PAV (https://github.com/EichlerLab/pav).

    Determine the number of perfect-homology bp upstream of an SV/indel using the SV/indel sequence (seq_sv), a contig
    or reference sequence (seq_tig) and the position of the first base upstream of the SV/indel (pos_tig) in 0-based
    coordinates. Both the contig and SV/indel sequence must be in the same orientation (reverse-complement if needed).
    Generally, the SV/indel sequence is in reference orientation and the contig sequence is the reference or an
    aligned contig in reference orientation (reverse-complemented if needed to get to the + strand).

    This function traverses from `pos_tig` to upstream bases in `seq_tig` using bases from the end of `seq_sv` until
    a mismatch between `seq_sv` and `seq_tig` is found. Search will wrap through `seq_sv` if homology is longer than
    the SV/indel.

    WARNING: This function assumes upper-case for the sequences. Differing case will break the homology search. If any
    sequence is None, 0 is returned.

    :param pos_tig: Contig/reference position (0-based) in reference orientation (may have been reverse-complemented by an
        alignment) where the homology search begins.
    :param seq_tig: Contig sequence as an upper-case string and in reference orientation (may have been reverse-
        complemented by the alignment).
    :param seq_sv: SV/indel sequence as an upper-case string.

    :return: Number of perfect-homology bases between `seq_sv` and `seq_tig` immediately upstream of `pos_tig`. If any
        of the sequneces are None, 0 is returned.
    """

    if seq_sv is None or seq_tig is None:
        return 0

    svlen = len(seq_sv)

    hom_len = 0

    while hom_len <= pos_tig:  # Do not shift off the edge of a contig.
        seq_tig_base = seq_tig[pos_tig - hom_len]

        # Do not match ambiguous bases
        if seq_tig_base not in {'A', 'C', 'G', 'T'}:
            break

        # Match the SV sequence (dowstream SV sequence with upstream reference/contig)
        if seq_sv[-((hom_len + 1) % svlen)] != seq_tig_base:
            # Circular index through seq in reverse from last base to the first, then back to the first
            # if it wraps around. If the downstream end of the SV/indel matches the reference upstream of
            # the SV/indel, shift left. For tandem repeats where the SV was placed in the middle of a
            # repeat array, shift through multiple perfect copies (% oplen loops through seq).
            break

        hom_len += 1

    # Return shifted amount
    return hom_len


def right_homology(pos_tig, seq_tig, seq_sv):
    """
    Duplicated from from PAV (https://github.com/EichlerLab/pav).

    Determine the number of perfect-homology bp downstream of an SV/indel using the SV/indel sequence (seq_sv), a contig
    or reference sequence (seq_tig) and the position of the first base downstream of the SV/indel (pos_tig) in 0-based
    coordinates. Both the contig and SV/indel sequence must be in the same orientation (reverse-complement if needed).
    Generally, the SV/indel sequence is in reference orientation and the contig sequence is the reference or an
    aligned contig in reference orientation (reverse-complemented if needed to get to the + strand).

    This function traverses from `pos_tig` to downstream bases in `seq_tig` using bases from the beginning of `seq_sv` until
    a mismatch between `seq_sv` and `seq_tig` is found. Search will wrap through `seq_sv` if homology is longer than
    the SV/indel.

    WARNING: This function assumes upper-case for the sequences. Differing case will break the homology search. If any
    sequence is None, 0 is returned.

    :param pos_tig: Contig/reference position (0-based) in reference orientation (may have been reverse-complemented by an
        alignment) where the homology search begins.
    :param seq_tig: Contig sequence as an upper-case string and in reference orientation (may have been reverse-
        complemented by the alignment).
    :param seq_sv: SV/indel sequence as an upper-case string.

    :return: Number of perfect-homology bases between `seq_sv` and `seq_tig` immediately downstream of `pos_tig`. If any
        of the sequences are None, 0 is returned.
    """

    if seq_sv is None or seq_tig is None:
        return 0

    svlen = len(seq_sv)
    tig_len = len(seq_tig)

    hom_len = 0
    pos_tig_limit = tig_len - pos_tig

    while hom_len < pos_tig_limit:  # Do not shift off the edge of a contig.
        seq_tig_base = seq_tig[pos_tig + hom_len]

        # Do not match ambiguous bases
        if seq_tig_base not in {'A', 'C', 'G', 'T'}:
            break

        # Match the SV sequence (dowstream SV sequence with upstream reference/contig)
        if seq_sv[hom_len % svlen] != seq_tig_base:
            # Circular index through seq in reverse from last base to the first, then back to the first
            # if it wraps around. If the downstream end of the SV/indel matches the reference upstream of
            # the SV/indel, shift left. For tandem repeats where the SV was placed in the middle of a
            # repeat array, shift through multiple perfect copies (% oplen loops through seq).
            break

        hom_len += 1

    # Return shifted amount
    return hom_len


def version_id(id_col, existing_id_set=None):
    """
    Take a column of IDs (Pandas Series object, `id_col`) and transform all duplicate IDs by appending "." and an
    integer so that no duplicate IDs remain.

    Example: If "chr1-1000-INS-10" has no duplicates, it will remain "chr1-1000-INS-10" in the output column. If it
    appears 3 times, then they will be named "chr1-1000-INS-10.1", "chr1-1000-INS-10.2", and "chr1-1000-INS-10.3".

    If an ID is duplicated and already versioned, the first appearance of the versioned ID will remain unchanged and the
    version will be incremented for subsequent appearances. If upon incrementing the name conflicts with another variant
    that was already versioned (e.g. a ".2" version already exists in the callset), the version will be incremented
    until it does not collide with any variant IDs.

    :param id_col: ID column as a Pandas Series object.
    :param existing_id_set: A set of existing variant IDs that must also be avoided. If any IDs in id_col match these
        IDs, they are altered as if the variant intersects another ID in id_col.

    :return: `id_col` unchanged if there are no duplicate IDs, or a new copy of `id_col` with IDs de-duplicated and
        versioned.
    """

    # Get counts
    id_count = collections.Counter()

    if existing_id_set is not None:
        id_count.update(existing_id_set)

    id_count.update(id_col)

    dup_set = {val for val, count in id_count.items() if count > 1}

    if len(dup_set) == 0:
        return id_col

    # Create a map: old name to new name
    id_col = id_col.copy()
    id_set = set(id_col) - dup_set

    if existing_id_set is not None:
        id_set |= existing_id_set

    for index in range(id_col.shape[0]):
        name = id_col.iloc[index]

        if name in dup_set:

            # Get current variant version (everything after "." if present, 1 by default)
            if not re.match(r'.*\.\d+$', name):
                name_version = 1
                tok = [name]

            else:
                try:
                    tok = name.rsplit('.', 1)
                    name_version = int(tok[1]) + 1

                except ValueError:
                    raise RuntimeError(f'Error de-duplicating variant ID field: Split "{name}" on "." and expected to find an integer at the end')

            # Find unique name
            new_name = '.'.join([tok[0], str(name_version)])

            while new_name in id_set:
                name_version += 1
                new_name = '.'.join([tok[0], str(name_version)])

            # Add to map
            id_col.iloc[index] = new_name
            id_set.add(new_name)

    # Append new variants
    return id_col


def check_unique_ids(df, message=''):
    """
    Check for unique IDs in a dataframe or an ID row. If IDs are not unique, throw a runtime error.

    :param df: Dataframe with an 'ID' column (pandas.DataFrame) or an ID column (pandas.Series).
    :param message: Prefix the exception message with this value if defined. (e.g. "{message}: Found X duplicate...").

    :return: No return value, throws an exception on failure, no effect if IDs are unique.
    """

    # Get ID row
    if issubclass(df.__class__, pd.DataFrame):
        if 'ID' not in df.columns:
            raise RuntimeError('Cannot check for unique IDs is DataFrame: No ID column')

        id_row = df['ID']

    elif issubclass(df.__class__, pd.Series):
        id_row = df

    else:
        raise RuntimeError(f'Unrecognized data type for check_unique_ids(): Expected Pandas DataFrame or Series: {df.__class__}')

    # Check for unique IDs
    if len(set(id_row)) < id_row.shape[0]:
        dup_id_set = [val for val, count in collections.Counter(id_row).items() if count > 1]

        if message is None:
            message = ''
        else:
            message = message.strip()

        raise RuntimeError(
            '{}Found {} duplicate variant IDs: {}{}'.format(
                message + ': ' if message else '',
                len(dup_id_set),
                ', '.join(dup_id_set[:3]),
                '...' if len(dup_id_set) > 3 else ''
            )
        )


def vcf_tsv_to_bed(tsv_in, sample=None, bed_file=None, filt_file=None, chunk_size=6000, threads=1, callback_pre_bed=None):
    """
    Convert a TSV generated by "bcftools query" to a formatted variant BED file.

    Function reads the TSV in chunks. If an open output file is given, records are written to it.

    The function call returns an iterator for processed chunks. Empty chunks are skipped, and if no record pass, then
    the iterator returns nothing.

    :param tsv_in: Input file name it TSV format. Typically generated by "bcftools query".
    :param bed_file: Open output file to write passed records or None to skip writing.
    :param filt_file: Open output file to write dropped records or None to skip writing.
    :param chunk_size: Number of lines in `tsv_in` to read in each iteration. Larger values use more memory, excessively
        low values will affect performance.
    :param threads: Number of threads to use while converting TSV fields to BED format.
    :param callback_pre_bed: Call this function before converting TSV fields to BED format. Each chunk is read and pre-
        processed before this function is called (i.e. column names are formatted, basic filtering has been done, and
        multiple ALTs are separated). This function takes one argument, the chunk dataframe, and returns a filtered
        dataframe. All records dropped by the function call are written to `filt_file`. If None is returned, then all
        records in the chunk are dropped and written to `filt_file`.

    :return: `None` if an open output file is given, otherwise, an iterator of processed chunks.
    """

    # Check arguments
    if type(chunk_size) == str:
        try:
            chunk_size = int(chunk_size)
        except ValueError as ex:
            raise ValueError(f'chunk_size: {ex}')

    elif type(chunk_size) == float:
        chunk_size = int(chunk_size)

    elif type(chunk_size) != int:
        raise ValueError(f'chunk_size is not an integer: {chunk_size}')

    if type(threads) == str:
        try:
            threads = int(threads)
        except ValueError as ex:
            raise ValueError(f'chunk_size: {ex}')

    elif type(threads) == float:
        threads = int(threads)

    elif type(threads) != int:
        raise ValueError(f'chunk_size is not an integer: {threads}')

    # Definitions
    COL_RENAME = {
        'CHROM': '#CHROM',
        'REF': 'VCF_REF',
        'POS': 'VCF_POS',
        'ALT': 'VCF_ALT'
    }

    REQUIRED_COLS = {'CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'GT'}

    # Read (as iterator)
    df_iter = pd.read_csv(tsv_in, sep='\t', header=0, low_memory=False, iterator=True, chunksize=chunk_size)

    # Process input in chunks
    col_names = None  # Input column names (after first regex transform)
    out_col_order = None  # Output column order

    id_set = set()  # Set of IDs already seen (used for ensuring IDs are unique)

    write_header = True  # Write header for this chunk (only true for first chunk)
    filt_header = True  # Write header for this chunk's dropped SV records (only true for first chunk)

    df = None

    for df in df_iter:

        # Remove prefixes added by bcftools (e.g. "# [1]CHROM" to "CHROM", "[2]POS" to "POS"
        df.columns = [re.sub('^#?\s*\[[^\]]+\]', '', col) for col in df.columns]

        if col_names is None:

            # Check samples. If single-sample VCF, accept it regardless of the sample name. If multi-sample VCF, must
            # contain sample (if not None) as a sample in the VCF.
            samples = {col.split(':', 1)[0] for col in df.columns if re.match('.+:.+', col)}

            if len(samples) > 1 and sample is not None:

                # Multi-samples and none match wildcards.sample
                if sample not in samples:
                    raise RuntimeError(
                        'Detected multiple samples in VCF with no samples "{sample}": {}'.format(
                            ', '.join(sample, sorted(samples))
                        )
                    )

                # Subset to non-sample columns (CHROM, POS, REF, ALT, etc) and sample columns (GT, etc)
                col_names = [
                    col for col in df.columns if not re.match('.+:.+', col) or col.startswith(sample + ':')
                ]

            else:
                col_names = df.columns

        # Subset to required columns
        df = df[col_names].copy()

        # Remove sample name from columns
        df.columns = [(col if ':' not in col else col.split(':', 1)[1]) for col in df.columns]

        if 'GT' not in df.columns:
            if len(samples) > 1:
                raise RuntimeError('Multi-sample VCF missing GT column')

            df['GT'] = '1/.'

        # Check columns and rename
        missing_cols = REQUIRED_COLS - set(df.columns)

        if missing_cols:
            raise RuntimeError('Missing VCF columns: {}'.format(', '.join(sorted(missing_cols))))

        df.columns = [COL_RENAME[col] if col in COL_RENAME else col for col in df.columns]

        # Set FILTER on QUAL if FILTER is missing
        df['FILTER'] = df.apply(qual_to_filter, axis=1)

        df_pass_filter = (df['FILTER'] == 'PASS') | (df['FILTER'] == '.')

        # Filter on FILTER column

        df_filt = df.loc[~ df_pass_filter].copy()

        if df_filt.shape[0] > 0 and filt_file is not None:
            df_filt['FILTER_REASON'] = 'Filter column'
            df_filt.to_csv(filt_file, sep='\t', index=False, header=filt_header)
            filt_file.flush()

            filt_header = False

        df = df.loc[df_pass_filter]

        if df.shape[0] == 0:
            continue

        # Save VCF index
        df['VCF_IDX'] = df.index

        # Separate multiple alleles
        df['VCF_ALT'] = df['VCF_ALT'].apply(lambda val: val.split(','))

        df = df.explode('VCF_ALT').reset_index(drop=True)

        df = df.loc[(df['VCF_ALT'] != '<NON_REF>') & (df['VCF_ALT'] != '.')]

        if df.shape[0] == 0:
            continue

        # Pre-parse callback
        if callback_pre_bed is not None:
            df_reformat = callback_pre_bed(df)

            if df_reformat is None:
                df_filt = df
            else:
                id_pass_set = set(df_reformat.index)
                df_filt = df.loc[[val for val in df.index if val not in id_pass_set]].copy()

            if df_filt.shape[0] > 0 and filt_file is not None:
                df_filt['FILTER_REASON'] = 'Dropped by callback_pre_bed'
                df_filt.to_csv(filt_file, sep='\t', index=False, header=filt_header)
                filt_file.flush()

                filt_header = False

            if df_reformat is None or df_reformat.shape[0] == 0:
                continue

            df = df_reformat

        # Get SV fields
        df_var_fields = svpoplib.pd.apply_parallel(
            df, vcf_fields_to_seq,
            n_part=500, n_core=threads,
            kwds={
                'pos_row': 'VCF_POS',
                'ref_row': 'VCF_REF',
                'alt_row': 'VCF_ALT'
            }
        )

        df = df[[col for col in df.columns if col not in df_var_fields.columns]]

        df = pd.concat([df, df_var_fields], axis=1)

        df = df.loc[df['AC'] > 0].copy()

        if df.shape[0] == 0:
            continue

        # Set index and arrange columns
        df['ALT'] = df['VCF_ALT']
        df['ID'] = get_variant_id(df)

        if out_col_order is None:
            out_col_order = list(order_variant_columns(df, tail_cols=(
            'QUAL', 'FILTER', 'VCF_POS', 'VCF_REF', 'VCF_ALT', 'VCF_IDX', 'SEQ')).columns)

        df = df.loc[:, out_col_order]

        # Rename duplicates
        df['ID'] = version_id(df['ID'], id_set)
        id_set |= set(df['ID'])

        # Write
        if bed_file is not None:
            df.to_csv(bed_file, sep='\t', index=False, header=write_header)
            bed_file.flush()

            write_header = False

        yield df

    # Write empty file header if no records were output
    if write_header and bed_file is not None:

        # Sanity check
        if df is None:
            raise RuntimeError('No dataframes generated by this rule (bug?)')

        if df.shape[0] != 0:
            raise RuntimeError(f'No header written, but the last dataframe is not empty: nrow={df.shape[0]} (bug?)')

        # Add required columns
        for col in ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN']:
            if col not in df.columns:
                df[col] = []

        df = order_variant_columns(df)

        # Write
        if bed_file is not None:
            df.to_csv(bed_file, sep='\t', index=False, header=write_header)
            bed_file.flush()

        yield df
