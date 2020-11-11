"""
Code for merging variant sets.
"""

import collections
import intervaltree
import multiprocessing
import numpy as np
import pandas as pd
import re
import traceback

import analib


def get_merge_def(def_name, config, default_none=False):
    """
    Get a merge definition string from a configured alias, `def_name`. Returns
    `config['merge_def'][def_name]` if it exists and `config` is not `None`, and returns `None` otherwise.

    :param def_name: Definition name to search.
    :param config: Configuration dictionary (or `None`).
    :param default_none: Return `None` if there is no definition alias (default is to return `def_name`).

    :return: Configuration definition or `def_name` if not found.
    """

    if config is None or 'merge_def' not in config:
        if default_none:
            return None

        return def_name

    return config['merge_def'].get(def_name, None if default_none else def_name)


def merge_variants(bed_list, sample_names, strategy, subset_chrom=None, threads=1):
    """
    Merge variants from multiple samples.

    :param bed_list: List of BED files to merge where each BED file is from one samples.
    :param sample_names: List of samples names. Each element matches an element at the same location in
        `bed_list`.
    :param strategy: Describes how to merge variants.
    :param subset_chrom: Merge only records from this chromosome. If `None`, merge all records.
    :param threads: Number of threads to use for intersecting variants.

    :return: A Pandas dataframe of a BED file with the index set to the ID column.
    """

    strategy_tok = strategy.split(':', 1)

    if len(bed_list) != len(sample_names):
        raise RuntimeError('Sample name list length ({}) does not match the input file list length ({})'.format(
            len(sample_names), len(bed_list))
        )

    if len(bed_list) == 0:
        raise RuntimeError('Cannot merge 0 samples')

    if len(strategy_tok) == 1:
        strategy_tok.append(None)

    if strategy_tok[0] == 'nr':
        return merge_variants_nr(bed_list, sample_names, strategy_tok[1], subset_chrom=subset_chrom, threads=threads)

    if strategy_tok[0] == 'fam':
        return merge_variants_fam(bed_list, sample_names, strategy_tok[1], subset_chrom=subset_chrom, threads=threads)

    if strategy_tok[0] == 'nrid':
        return merge_variants_nrid(bed_list, sample_names, strategy_tok[1], subset_chrom=subset_chrom, threads=threads)

    else:
        raise RuntimeError('Unrecognized strategy: {}'.format(strategy))


def merge_variants_nr(bed_list, sample_names, merge_params, subset_chrom=None, threads=1):
    """
    Merge all non-redundant variants from multiple samples.

    Recognized parameters:
    * ro: Match variants by reciprocal overlap.
    * szro: Match variants minimum reciprocal size overlap (max value for ro). This is like ro, but allows the variants
        to be offset. Size overlap is calculated only on variant sizes, and a maximum offset (offset parameter)\
        determines how far apart the variants may be. If szro is specified and ro is not, then ro is implied with the
        same overlap threshold and tried before size/offset overlap.
    * offset: Maximum offset (minimum of start position difference or end position difference)
    * ref, alt, and refalt: If specified, the REF and/or ALT columns must match. This is an additional restriction on
        other merging parameters (is not a merging strategy in itself).

    Merging each BED occurs in these stages:
        1) ID. Same ID is always matched first (no parameters needed to specify). IDs in SV-Pop are assumed to be
            a descriptor of the variant including chrom, position, size. If IDs match, REF and ALT are assumed to also
            match.
        2) RO: If ro and/or szro is set, then a strict recprocal-overlap is run. The minimum overlap is ro if it was
            defined and szro if ro was not defined.
        3) Size-overlap + offset: Try intersecting by reciprocal-size-overlap with a maximum offset.

    :param bed_list: List of BED files to merge where each BED file is from one samples.
    :param sample_names: List of samples names. Each element matches an element at the same location in
        `bed_list`.
    :param merge_params: Overlap percentage or None to merge by the default (50% overlap).
    :param subset_chrom: Merge only records from this chromosome. If `None`, merge all records.
    :param threads: Number of threads to use for intersecting variants.

    :return: A Pandas dataframe of a BED file with the index set to the ID column.
    """

    # Check BED
    if len(bed_list) == 0:
        raise ValueError('BED file list is empty')

    # Check sample names
    if sample_names is None or len(sample_names) == 0:
        raise RuntimeError('Sample names is missing or empty')

    sample_names = [val.strip() for val in sample_names]

    if not all(bool(val) for val in sample_names):
        raise RuntimeError('Found empty sample names')

    if any([re.match('.*\\s.*', val) is not None for val in sample_names]):
        raise RuntimeError('Error: Sample names contain whitespace')

    n_samples = len(sample_names)

    if len(set(sample_names)) != n_samples:
        raise RuntimeError('Sample names may not be duplicated')

    # Parse parameters
    ro_min = None
    szro_min = None
    offset_max = None
    match_ref = False
    match_alt = False
    expand_base = False

    if merge_params is None:
        raise RuntimeError('Cannot merge (strategy=nr) with parameters: None')

    for param_element in merge_params.split(':'):

        # Tokenize
        param_element = param_element.strip()

        if not param_element:
            continue

        param_tok = re.split('\s*=\s*', param_element, 1)

        key = param_tok[0].lower()

        if len(param_tok) > 1:
            val = param_tok[1]
        else:
            val = None

        # Process key
        if key == 'ro':

            if val is None:
                raise ValueError('Missing value for parameter "ro" (e.g. "ro=50"): {}'.format(merge_params))

            if val != 'any':
                ro_min = int(val.strip()) / 100

                if ro_min < 0 or ro_min > 1:
                    raise ValueError(
                        'Overlap length (ro) must be between 0 and 100 (inclusive): {}'.format(param_element)
                    )

            else:
                raise RuntimeError('RO "any" is not yet implemented')
                #ro_min = None

        elif key == 'szro':
            if val is None:
                raise ValueError('Missing value for parameter "szro" (e.g. "szro=50"): {}'.format(merge_params))

            if val != 'any':
                szro_min = int(val.strip()) / 100

                if szro_min < 0 or szro_min > 1:
                    raise ValueError(
                        'Overlap length (szro) must be between 0 and 100 (inclusive): {}'.format(param_element)
                    )

            else:
                szro_min = None

        elif key == 'offset':
            if val is None:
                raise ValueError('Missing value for parameter "offset" (maxiumum offset, e.g. "offset=2000"): {}'.format(merge_params))

            if val != 'any':
                offset_max = int(val.strip())

                if offset_max < 0:
                    raise RuntimeError('Maximum offset (offset parameter) may not be negative: {}'.format(merge_params))

            else:
                offset_max = None

        elif key == 'refalt':
            if val is not None:
                raise RuntimeError('Match-REF/ALT (refalt) should not have an argument')

            match_ref = True
            match_alt = True

        elif key == 'ref':
            if val is not None:
                raise RuntimeError('Match-REF (ref) should not have an argument')

            match_ref = True

        elif key == 'alt':
            if val is not None:
                raise RuntimeError('Match-ALT (alt) should not have an argument')

            match_alt = True

        elif key == 'expand':
            if val is not None:
                raise RuntimeError('Expand-base (expand) should not have an argument')

            expand_base = True

        else:
            raise ValueError('Unknown parameter token in: {}'.format(key))

    # Check parameters
    if szro_min is not None and offset_max is None:
        raise RuntimeError('Parameters "szro" was specified without "offset"')

    # Get merge size threshold
    if ro_min is None and szro_min is not None:
        ro_min = szro_min

    # Note:
    # The first time a variant is found, it is added to the merge set (df), and it is supported by itself. The overlap
    # score is defined as -1. As variants are added that intersect with a variant already in the merge set, SAMPLE and
    # ID are copied from that variant, the supporting variant is noted in SUPPORT_ID and SUPPORT_SAMPLE, and the overlap
    # score is between 0 and 1 (depending on how well they overlap).

    # Initialize a table of variants with the first sample. All variants from this sample are in the merged set.
    sample_name = sample_names[0]

    print('Merging: {}'.format(sample_name))

    df = pd.read_csv(
        bed_list[0], sep='\t', header=0,
        usecols=lambda col: col in {'#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'REF', 'ALT'}
    )

    if 'SVLEN' not in df.columns:
        df['SVLEN'] = df['END'] - df['POS']

    if 'SVTYPE' not in df.columns:
        df['SVTYPE'] = 'RGN'

    df = df.loc[:,
         [col for col in ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'REF', 'ALT') if col in df.columns]]

    if np.any(df['SVLEN'] < 0):
        raise RuntimeError('Negative SVLEN entries in {}'.format(bed_list[0]))

    df.sort_values(['#CHROM', 'POS', 'SVLEN', 'ID'], inplace=True)

    # Check for missing columns
    missing_cols = [col for col in ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN') if col not in df.columns]

    if missing_cols:
        raise RuntimeError('Missing column(s) for source {}: {}: {}'.format(
            sample_name,
            ', '.join(missing_cols),
            bed_list[0]
        ))

    # Subset chromosomes
    if subset_chrom is not None:
        df = df.loc[df['#CHROM'] == subset_chrom]

    # Check for REF and ALT
    if match_ref and 'REF' not in df.columns:
        raise RuntimeError(
            'Parameter Match-REF (ref or refalt) is set and table is missing the REF column: {}'.format(bed_list[0])
        )

    if match_alt and 'ALT' not in df.columns:
        raise RuntimeError(
            'Parameter Match-ALT (alt or refalt) is set and table is missing the ALT column: {}'.format(bed_list[0])
        )

    # Check for unique IDs
    dup_names = [name for name, count in collections.Counter(df['ID']).items() if count > 1]

    if dup_names:
        dup_name_str = ', '.join(dup_names[:3])

        if len(dup_names) > 3:
            dup_name_str += ', ...'

        raise RuntimeError('Found {} IDs with duplicates in sample {}: {}'.format(
            len(dup_names), sample_name, dup_name_str
        ))

    # Add tracking columns
    df['SUPPORT_ID'] = df['ID']

    df['SAMPLE'] = sample_name
    df['SUPPORT_SAMPLE'] = sample_name

    df['SUPPORT_OFFSET'] = -1
    df['SUPPORT_RO'] = -1
    df['SUPPORT_SZRO'] = -1
    df['SUPPORT_OFFSZ'] = -1

    # Set index
    df.set_index(df['ID'], inplace=True, drop=False)
    df.index.name = 'INDEX'

    # Setup a dictionary to translate support sample table columns to the merged table columns
    support_col_rename = {
        'OFFSET': 'SUPPORT_OFFSET',
        'RO': 'SUPPORT_RO',
        'SZRO': 'SUPPORT_SZRO',
        'OFFSZ': 'SUPPORT_OFFSZ',
    }

    # Add each variant to the table
    base_support = list()  # Table of supporting variants if they are not part of subsequent rounds of intersection (defined by "expand").

    for index in range(1, len(bed_list)):

        ## Read ##

        # Read next variant table
        sample_name = sample_names[index]

        print('Merging: {}'.format(sample_name))

        df_next = pd.read_csv(
            bed_list[index], sep='\t', header=0,
            usecols=lambda col: col in {'#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'REF', 'ALT'}
        )

        if 'SVLEN' not in df_next.columns:
            df_next['SVLEN'] = df_next['END'] - df_next['POS']

        if 'SVTYPE' not in df_next.columns:
            df_next['SVTYPE'] = 'RGN'

        if np.any(df_next['SVLEN'] < 0):
            raise RuntimeError('Negative SVLEN entries in {}'.format(bed_list[index]))

        df_next = df_next.loc[:,
             [col for col in ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'REF', 'ALT') if col in df_next.columns]]

        df_next.sort_values(['#CHROM', 'POS', 'SVLEN', 'ID'], inplace=True)

        # Check for missing columns
        missing_cols = [col for col in ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN') if col not in df_next.columns]

        if missing_cols:
            raise RuntimeError('Missing column(s) for source {}: {}: {}'.format(
                sample_name,
                ', '.join(missing_cols),
                bed_list[index]
            ))

        # Subset chromosome
        if subset_chrom is not None:
            df_next = df_next.loc[df_next['#CHROM'] == subset_chrom]

        # Check for REF and ALT
        if match_ref and 'REF' not in df_next.columns:
            raise RuntimeError(
                'Parameter Match-REF (ref or refalt) is set and table is missing the REF column: {}'.format(bed_list[index])
            )

        if match_alt and 'ALT' not in df_next.columns:
            raise RuntimeError(
                'Parameter Match-ALT (alt or refalt) is set and table is missing the ALT column: {}'.format(bed_list[index])
            )

        # Check for unique IDs
        dup_names = [name for name, count in collections.Counter(df_next['ID']).items() if count > 1]

        if dup_names:
            dup_name_str = ', '.join(dup_names[:3])

            if len(dup_names) > 3:
                dup_name_str += ', ...'

            raise RuntimeError('Found {} IDs with duplicates in sample {}: {}'.format(
                len(dup_names), sample_name, dup_name_str
            ))

        # Set IDs
        #df_next['ORG_ID'] = df_next['ID']
        #df_next['ID'] = df_next['ID'].apply(lambda val: sample_name + '_' + val)

        # Set index
        df_next.set_index('ID', inplace=True, drop=False)
        df_next.index.name = 'INDEX'


        ## Build intersect support table ##

        support_table_list = list()

        # Create copies of df and df_next that can be subset (variants removed by each phase of intersection)
        df_sub = df.copy()
        df_next_sub = df_next.copy()

        # INTERSECT: by ID
        support_table_list.append(get_support_table_nrid(df_sub, df_next_sub))

        id_set = set(support_table_list[-1]['ID'])
        id_next_set = set(support_table_list[-1]['TARGET_ID'])

        df_sub = df_sub.loc[df_sub['ID'].apply(lambda val: val not in id_set)]
        df_next_sub = df_next_sub.loc[df_next_sub['ID'].apply(lambda val: val not in id_next_set)]

        # INTERSECT: RO (ro or szro defined)
        if ro_min is not None and df_sub.shape[0] > 0 and df_next_sub.shape[0] > 0:

            df_support_ro = get_support_table(
                df_sub, df_next_sub, threads, None, ro_min, match_ref, match_alt
            )

            id_set = set(df_support_ro['ID'])
            id_next_set = set(df_support_ro['TARGET_ID'])

            df_sub = df_sub.loc[df_sub['ID'].apply(lambda val: val not in id_set)]
            df_next_sub = df_next_sub.loc[df_next_sub['ID'].apply(lambda val: val not in id_next_set)]

            support_table_list.append(df_support_ro)

            del(df_support_ro)

        # INTERSECT: SZRO + OFFSET
        if szro_min is not None and df_sub.shape[0] > 0 and df_next_sub.shape[0] > 0:

            df_support_szro = get_support_table(
                df_sub, df_next_sub, threads, offset_max, szro_min, match_ref, match_alt
            )

            df_support_szro = df_support_szro.loc[(df_support_szro['OFFSET'] <= offset_max) | (df_support_szro['SZRO'] >= szro_min)]

            id_set = set(df_support_szro['ID'])
            id_next_set = set(df_support_szro['TARGET_ID'])

            df_sub = df_sub.loc[df_sub['ID'].apply(lambda val: val not in id_set)]
            df_next_sub = df_next_sub.loc[df_next_sub['ID'].apply(lambda val: val not in id_next_set)]

            support_table_list.append(df_support_szro)

            del(df_support_szro)

        # Concat df_support
        df_support = pd.concat(support_table_list)

        del(df_sub)
        del(df_next_sub)

        # Annotate support
        if df_support.shape[0] > 0:

            # Set Index
            df_support.set_index('TARGET_ID', inplace=True, drop=False)
            df_support.index.name = 'INDEX'

            # Add sample name
            df_support['SAMPLE'] = sample_name
            df_support['SUPPORT_SAMPLE'] = list(df.loc[df_support['ID'], 'SAMPLE'])

            # Fix IDs (ID = new ID, SUPPORT_ID = ID of sample it supports)
            df_support['SUPPORT_ID'] = list(df.loc[df_support['ID'], 'SUPPORT_ID'])
            df_support['ID'] = df_support['TARGET_ID']

            # Rename columns (offset merger column names to support column names)
            df_support.columns = [support_col_rename.get(col, col) for col in df_support.columns]

            # Remove redundant support (if present).
            df_support = df_support.sort_values(
                ['SUPPORT_RO', 'SUPPORT_OFFSET', 'SUPPORT_SZRO'],
                ascending=[False, True, False]
            ).drop_duplicates('ID', keep='first')

            # Arrange columns
            df_support['#CHROM'] = df_next['#CHROM']
            df_support['POS'] = df_next['POS']
            df_support['END'] = df_next['END']
            df_support['SVTYPE'] = df_next['SVTYPE']
            df_support['SVLEN'] = df_next['SVLEN']

            if 'REF' in df_next.columns:
                df_support['REF'] = df_next['REF']

            if 'ALT' in df_next.columns:
                df_support['ALT'] = df_next['ALT']

            df_support = df_support.loc[:, df.columns]

            # Append to existing supporting variants (expand) or save to a list of support (no expand)
            if expand_base:
                df = pd.concat([df, df_support], axis=0)
            else:
                base_support.append(df_support)

        # Read new variants from this sample (variants that do not support an existing call)
        support_id_set = set(df_support['ID'])

        df_new = df_next.loc[df_next['ID'].apply(lambda val: val not in support_id_set)].copy()

        if df_new.shape[0] > 0:
            df_new['SUPPORT_ID'] = df_new['ID']
            df_new['SAMPLE'] = sample_name
            df_new['SUPPORT_SAMPLE'] = sample_name
            df_new['SUPPORT_OFFSET'] = -1
            df_new['SUPPORT_RO'] = -1
            df_new['SUPPORT_SZRO'] = -1
            df_new['SUPPORT_OFFSZ'] = -1

            # Refuse to merge if it creates ID conflicts
            dup_names = sorted(set(df['ID']) & set(df_new['ID']))

            if dup_names:
                dup_name_str = ', '.join(dup_names[:3])

                if len(dup_names) > 3:
                    dup_name_str += ', ...'

                raise RuntimeError('Found {} duplicate IDs when merging new variants from sample {}: {}'.format(
                    len(dup_names), sample_name, dup_name_str
                ))

            # Append new variants
            df = pd.concat([df, df_new.loc[:, df.columns]], axis=0)

        # Sort
        df.sort_values(['#CHROM', 'POS', 'SVLEN', 'ID'], inplace=True)

    # Merge support variants into df (where support tables are kept if df is not expanded)
    if len(base_support) > 0:

        if len(base_support) > 1:
            df_base_support = pd.concat(base_support, axis=0)
        else:
            df_base_support = base_support[0]

        df = pd.concat([df, df_base_support], axis=0)

    # Finalize merged variant set
    if df.shape[0] > 0:

        # Make SAMPLE and SUPPORT_SAMPLE categorical (sort in the same order as they were merged)
        df['SAMPLE'] = pd.Categorical(df['SAMPLE'], sample_names)
        df['SUPPORT_SAMPLE'] = pd.Categorical(df['SUPPORT_SAMPLE'], sample_names)

        # Sort by support (best support first)
        df['SUPPORT_OFFSET'] = df['SUPPORT_OFFSET'].apply(lambda val: np.max([0, val]))
        df['SUPPORT_RO'] = np.abs(df['SUPPORT_RO'])
        df['SUPPORT_SZRO'] = np.abs(df['SUPPORT_SZRO'])
        df['SUPPORT_OFFSZ'] = np.abs(df['SUPPORT_OFFSZ'])

        df['IS_PRIMARY'] = df['SAMPLE'] == df['SUPPORT_SAMPLE']

        df.sort_values(
            ['SAMPLE', 'SUPPORT_RO', 'SUPPORT_OFFSET', 'SUPPORT_SZRO', 'SUPPORT_OFFSZ'],
            ascending=[True, False, True, False, False],
            inplace=True
        )

        # Find best support variant for each mergeset variant (per sample)
        df.drop_duplicates(['SUPPORT_ID', 'SAMPLE', 'SUPPORT_SAMPLE'], keep='first', inplace=True)

        # Re-sort by ID then SAMPLE (for concatenating stats in order)
        df.sort_values(['SUPPORT_ID', 'SAMPLE', 'SUPPORT_SAMPLE'], inplace=True)

        df_support = df.groupby(
            'SUPPORT_ID'
        ).apply(lambda subdf: pd.Series(
            [
                subdf.iloc[0]['SUPPORT_SAMPLE'],
                #subdf.loc[subdf['SUPPORT_ID']].iloc[0].squeeze()['ORG_ID'],
                subdf.iloc[0].squeeze()['ID'],

                subdf.shape[0],
                '{:.4f}'.format(subdf.shape[0] / n_samples),

                ','.join(subdf['SAMPLE']),
                ','.join(subdf['ID']),

                ','.join(['{:.2f}'.format(val) for val in subdf['SUPPORT_RO']]),
                ','.join(['{:.0f}'.format(val) for val in subdf['SUPPORT_OFFSET']]),
                ','.join(['{:.2f}'.format(val) for val in subdf['SUPPORT_SZRO']]),
                ','.join(['{:.2f}'.format(val) for val in subdf['SUPPORT_OFFSZ']]),
            ],
            index=[
                'MERGE_SRC', 'MERGE_SRC_ID',
                'MERGE_AC', 'MERGE_AF',
                'MERGE_SAMPLES', 'MERGE_VARIANTS',
                'MERGE_RO', 'MERGE_OFFSET', 'MERGE_SZRO', 'MERGE_OFFSZ'
            ]
        ))

    else:
        df_support = pd.DataFrame(
            [],
            columns=[
                'MERGE_SRC', 'MERGE_SRC_ID', 'MERGE_AC', 'MERGE_AF',
                'MERGE_SAMPLES', 'MERGE_VARIANTS',
                'MERGE_RO', 'MERGE_OFFSET', 'MERGE_SZRO', 'MERGE_OFFSZ'
            ]
        )

    # Merge original BED files by by df_support
    return merge_sample_by_support(df_support, bed_list, sample_names)


def merge_variants_fam(bed_list, sample_names, merge_params, subset_chrom=None, threads=1):
    """
    Merge all non-redundant variants from multiple samples and perform family-wise annotations.

    :param bed_list: List of BED files to merge where each BED file is from one samples.
    :param sample_names: List of samples names. Each element matches an element at the same location in
        `bed_list`.
    :param merge_params: A colon-delimited list of parameters. Must contain "ro=" for the reciprocal-overlap parameter
        (for non-redundant merge). Must also contain "order=" with a value containing "c" (child), "m" (mother), and
        "f" (father) describing the family members in the bed and sample lists (e.g. "ccmf" means that the samples are
        ordered as child-child-mother-father). A typical parameter list might be "ro=50:ccmf".
    :param subset_chrom: Merge only records from this chromosome. If `None`, merge all records.
    :param threads: Number of threads to use for intersecting variants.

    :return: A Pandas dataframe of a BED file with the index set to the ID column.
    """

    # Merge function aliases
    merge_func_dict = {
        'nr': merge_variants_nr,
        'nrid': merge_variants_nrid
    }

    # Set defaults
    fam_order = None   # Order of samples in bed_list and sample_names. String with "c", "f", and "m" (e.g. "ccmf" for child-child-mother-father).
    merge_func_name = 'nr'  # Algorithm to use for merging.

    # Get parameters
    if merge_params is None:
        raise RuntimeError('Cannot merge family with parameters: None')

    nr_param_list = list()  # List of parameters to be given to the nr merge function

    for param_element in merge_params.split(':'):

        # Tokenize
        param_element = param_element.strip()

        if not param_element:
            continue

        param_tok = re.split('\s*=\s*', param_element, 1)

        if len(param_tok) > 1:
            key, val = param_element.split('=', 1)
        else:
            key = param_element
            val = None

        key = key.lower()

        # Order
        if key == 'order':
            fam_order = val.lower()

        elif key == 'famalg':
            merge_func_name = val.lower()

        else:
            nr_param_list.append(param_element)

    # Check parameters
    #if len(nr_param_list) == 0:
    #    raise RuntimeError('merge_variants_fam(): No parameters for nr merge function')

    if fam_order is None:
        raise RuntimeError('merge_variants_fam(): Missing "order" in parameter list')

    # Get merging algorithm
    if merge_func_name not in merge_func_dict:
        raise RuntimeError('Unrecognized merge function for family-wise merging (famalg=): ' + merge_func_name)

    merge_func = merge_func_dict[merge_func_name]

    # Get parameter string
    nr_param = ':'.join(nr_param_list)

    # Check family
    fam_list = list(fam_order)

    if len(fam_order) != len(sample_names):
        raise RuntimeError(
            'merge_variants_fam(): Number of "order" elements ({}) does not match the number of sample names ({})'.format(
                len(fam_order), len(sample_names)
            ))

    if any([val not in {'c', 'm', 'f'} for val in fam_list]):
        raise RuntimeError(
            'merge_variants_fam(): Found characters that are not in {{c, m, f}} in family order: {}'.format(fam_order)
        )

    if sum([val == 'f' for val in fam_list]) != 1:
        raise RuntimeError(
            'merge_variants_fam(): Expected 1 father (f) in family order: Found {}'.format(sum([val == 'f' for val in fam_list]))
        )

    if sum([val == 'm' for val in fam_list]) != 1:
        raise RuntimeError(
            'merge_variants_fam(): Expected 1 mother (m) in family order: Found {}'.format(sum([val == 'm' for val in fam_list]))
        )

    # Get sample names
    mo_index = [index for val, index in zip(fam_order, range(len(fam_order))) if val == 'm'][0]
    fa_index = [index for val, index in zip(fam_order, range(len(fam_order))) if val == 'f'][0]
    ch_index_list = [index for val, index in zip(fam_order, range(len(fam_order))) if val not in {'m', 'f'}]

    mo_name = sample_names[mo_index]
    fa_name = sample_names[fa_index]

    # Do non-redundant merge
    df = merge_func(
        bed_list=bed_list, sample_names=sample_names, merge_params=nr_param, subset_chrom=subset_chrom, threads=threads
    )

    # Do family-wise annotations
    df_samples = df['MERGE_SAMPLES'].apply(lambda val: val.split(','))

    anno_dict = {
        (False, False): 'DENOVO',
        (False, True): 'FA',
        (True, False): 'MO',
        (True, True): 'MOFA'
    }

    for ch_index in ch_index_list:
        child_name = sample_names[ch_index]

        df['INHERIT_{}'.format(child_name)] = df_samples.apply(
            lambda vals: np.nan if child_name not in vals else anno_dict[(mo_name in vals, fa_name in vals)]
        )

    # Return dataframe
    return df


def merge_variants_nrid(bed_list, sample_names, merge_params, subset_chrom=None, threads=1):
    """
    merge variants by ID. This requires exact matches among the variants.

    :param bed_list: List of BED files to merge.
    :param sample_names: List of sample names.
    :param merge_params: Merge parameters.
    :param subset_chrom: Subset this chromosome.
    :param threads: Number of threads (ignored).

    :return: Merged dataframe.
    """

    # Parse parameters
    if merge_params is None:
        merge_params = ''

    merge_params = merge_params.strip()

    if merge_params:
        raise RuntimeError('Unrecognized options to merged by ID: ' + merge_params)

    # Get a list of supporting sample for each ID
    n_samples = len(sample_names)

    support_sample_list = collections.defaultdict(list)

    for index in range(n_samples):
        sample_name = sample_names[index]

        df_sample = pd.read_csv(
            bed_list[index],
            sep='\t',
            usecols=('#CHROM', 'ID')
        )

        if subset_chrom is not None:
            df_sample = df_sample.loc[df_sample['#CHROM'] == subset_chrom]

        for variant_id in set(df_sample['ID']):
            support_sample_list[variant_id].append(sample_name)

    # Make support table
    df_support_list = list()

    for variant_id, sample_list in support_sample_list.items():
        df_support_list.append(pd.Series(
            [
                sample_list[0],
                variant_id,
                len(sample_list),
                len(sample_list) / n_samples,
                ','.join(sample_list),
                ','.join([variant_id] * len(sample_list))
            ],
            index=['MERGE_SRC', 'MERGE_SRC_ID', 'MERGE_AC', 'MERGE_AF', 'MERGE_SAMPLES', 'MERGE_VARIANTS']
        ))

    del support_sample_list

    if df_support_list:
        df_support = pd.concat(df_support_list, axis=1).T
    else:
        df_support = pd.DataFrame(
            [],
            columns=['MERGE_SRC', 'MERGE_SRC_ID', 'MERGE_AC', 'MERGE_AF', 'MERGE_SAMPLES', 'MERGE_VARIANTS']
        )

    del df_support_list

    return merge_sample_by_support(df_support, bed_list, sample_names)


def merge_annotations(df_merge, anno_tab_list, sample_names, sort_columns=None):
    """
    Merge a table of annotations from several samples into one table.

    :param df_merge: Table of merged structural variants (BED file). Must contain columns 'ID',
        'MERGE_SRC', and 'MERGE_SRC_ID'.
    :param anno_tab_list: List of annotations tables to be merged. Must have an 'ID' column.
    :param sample_names: Names of the samples to be merged. Each element is the name for the data in
        `anno_tab_list` at the same index.
    :param sort_columns: List of column names to sort the merged annotations by or `None` to leave it unsorted.

    :return: Dataframe of merged annotations.
    """

    # Check arguments
    if len(anno_tab_list) != len(sample_names):
        raise RuntimeError('Sample name list length ({}) does not match the input file list length ({})'.format(
            len(sample_names), len(anno_tab_list))
        )

    if len(anno_tab_list) == 0:
        raise RuntimeError('Cannot merge 0 samples')

    df_list = list()

    # Subset from each samples
    for index in range(len(sample_names)):

        # Get samples name and input file name
        sample_name = sample_names[index]
        anno_tab = anno_tab_list[index]

        # Get table of annotations
        df_anno = pd.read_csv(anno_tab, sep='\t', header=0)
        df_anno.index = df_anno['ID']

        # Subset
        df_merge_subset = df_merge.loc[df_merge['MERGE_SRC'] == sample_name]
        id_dict = {row[1]['MERGE_SRC_ID']: row[1]['ID'] for row in df_merge_subset.iterrows()}

        id_subset = set(df_merge_subset['MERGE_SRC_ID'])

        df_anno = df_anno.loc[df_anno['ID'].apply(lambda svid: svid in id_subset)]

        df_anno['ID'] = df_anno['ID'].apply(lambda svid: id_dict[svid])
        df_anno.index = df_anno['ID']

        df_list.append(df_anno)

    # Merge subsets
    df_anno = pd.concat(df_list, axis=0)

    # Resort
    if sort_columns is not None:
        df_anno.sort_values(list(sort_columns), inplace=True)

    # Return
    return df_anno


def get_samples_for_mergeset(mergeset, config):
    """
    Get a list of samples with source variants for a set of merged samples.

    :param mergeset: Merged variant set name.
    :param config: Config (loaded by Snakemake).

    :return: List of samples that must be loaded for `mergeset`.
    """

    # Get list of samples
    sample_list = config['svmerge'][mergeset].get('svsource', None)

    if sample_list is None:
        sample_set_name = config['svmerge'][mergeset]['sampleset']
        sample_list = config['sampleset'][sample_set_name]

        if len(sample_list) == 0:
            raise RuntimeError('Empty input for merge set {} after matching full samples set'.format(mergeset))

    else:
        if len(sample_list) == 0:
            raise RuntimeError('Empty input for merge set {}'.format(mergeset))

    # Return list
    return sample_list


def get_disc_class_by_row(row):
    """
    Get discovery class.

    :param row: A series with fields "MERGE_AF" and "MERGE_AC".

    :return: A class (as string).
    """

    if row['MERGE_AF'] == 1:
        return 'SHARED'

    if row['MERGE_AF'] >= 0.5:
        return 'MAJOR'

    if row['MERGE_AC'] > 1:
        return 'POLY'

    return 'SINGLE'


def get_disc_class(df):
    """
    Get discovery class.

    :param df: A dataframe with columns "MERGE_AF" and "MERGE_AC" or a series with those fields.

    :return: A class (as string) for each row (if dataframe) or a class (as string) if series.
    """

    # If df is series, get value for one row
    if df.__class__ is pd.core.series.Series:
        return get_disc_class_by_row(df)

    # Apply to all rows
    return df.apply(get_disc_class_by_row, axis=1)


def merge_sample_by_support(df_support, bed_list, sample_names):
    """
    Take a support table and generate the merged variant BED file.
    
    The support table must have at least columns:
    1) MERGE_SRC: Sample variant should be extracted from
    2) MERGE_SRC_ID: ID of of the variant to extract from the source. This becomes the record that represents all
        variants merged with it.
    
    The support table may also have:
    1) MERGE_AC: Number of samples variant was found in.
    2) MERGE_AF: MERGE_AC divided by the number of samples merged.
    3) MERGE_SAMPLES: A comma-separated list of samples. If present, number of samples must match MERGE_AC, and the
        first sample must be MERGE_SRC.
    4) MERGE_VARIANTS: A comma-separated list of variants from each sample. If present, number must match MERGE_AC and
        be in the same order as MERGE_SAMPLES.
    5) MERGE_RO: A comma-separated list of reciprocal-overlap values with the variant from each sample against the
        representitive variant. If present, number must match MERGE_AC and be in the same order as MERGE_SAMPLES.
    6) MERGE_OFFSET: A comma-separated list of offset distances with the variant from each sample against the
        representitive variant. If present, number must match MERGE_AC and be in the same order as MERGE_SAMPLES.
    7) MERGE_SZRO: A comma-separated list of size-reciprocal-overlap values with the variant from each sample against
        the representitive variant. If present, number must match MERGE_AC and be in the same order as MERGE_SAMPLES.
    8) MERGE_OFFSZ: A comma-separated list of offset/size values with the variant from each sample against the
        representitive variant. If present, number must match MERGE_AC and be in the same order as MERGE_SAMPLES.

    Other columns are silently ignored.
    
    :param df_support: Variant support table (see above for format). If `None`, treats this as as empty table and
        generates a merged dataframe with no variants.
    :param bed_list: List of input BED files.
    :param sample_names: List of sample names. Must be the same length as `bed_list` where values correspond (name at
        index X is the name of the sample in `bed_list` at index X).
    
    :return: A merged dataframe of variants.
    """

    REQUIRED_COLUMNS = ['MERGE_SRC', 'MERGE_SRC_ID']

    OPT_COL = ['MERGE_AC', 'MERGE_AF', 'MERGE_SAMPLES', 'MERGE_VARIANTS', 'MERGE_RO', 'MERGE_OFFSET', 'MERGE_SZRO', 'MERGE_OFFSZ']

    OPT_COL_DTYPE = {
        'MERGE_AC': np.int32,
        'MERGE_AF': np.float32
    }

    # Check values
    if df_support is None:
        df_support = pd.DataFrame([], columns=REQUIRED_COLUMNS)

    if len(bed_list) != len(sample_names):
        raise RuntimeError(
            'BED list and sample name list lengths differ ({} vs {})'.format(len(bed_list), len(sample_names))
        )

    if df_support.shape[0] > 0:
        missing_cols = [col for col in REQUIRED_COLUMNS if col not in df_support.columns]
    else:
        missing_cols = []

    if missing_cols:
        raise RuntimeError('Missing column(s) in df_support: ' + ', '.join(missing_cols))

    # Get optional columns that are present
    opt_columns = [opt_col for opt_col in OPT_COL if opt_col in df_support.columns]

    # Merged variant IDs should be unique
    dup_id_list = [val for val, count in collections.Counter(df_support['MERGE_SRC_ID']).items() if count > 1]

    if len(dup_id_list) > 0:
        dup_id_str = ', '.join(dup_id_list[:3]) + (', ...' if len(dup_id_list) > 3 else '')
        raise RuntimeError('Duplicate IDs after merging ({}): {}'.format(len(dup_id_list), dup_id_str))

    # Merge original SVs
    col_list = list()  # List of columns in the final output
    merge_df_list = list()  # List of dataframes to be merged into final output

    for index in range(len(bed_list)):
        bed_file_name = bed_list[index]
        sample_name = sample_names[index]

        # Get merged variants for this sample
        df_support_sample = df_support.loc[df_support['MERGE_SRC'] == sample_name].copy()

        df_support_sample.set_index('MERGE_SRC_ID', inplace=True, drop=False)
        df_support_sample.index.name = 'INDEX'

        # Read variants from sample
        df_sample = pd.read_csv(bed_file_name, sep='\t', header=0)

        # Update columns
        col_list.extend([col for col in df_sample.columns if col not in col_list])

        # Get variant names
        if df_support_sample.shape[0] == 0:
            continue

        df_sample.set_index('ID', inplace=True, drop=False)
        df_sample.index.name = 'INDEX'

        df_sample = df_sample.loc[df_support_sample['MERGE_SRC_ID']]

        # Transfer required columns
        df_sample['MERGE_SRC'] = df_support_sample['MERGE_SRC']
        df_sample['MERGE_SRC_ID'] = df_support_sample['MERGE_SRC_ID']

        # Transfer optional columns annotations
        for col in opt_columns:
            df_sample[col] = df_support_sample[col]

        # Append to list
        merge_df_list.append(df_sample)

    # Append reqired and optional columns to the end
    col_list.extend([col for col in REQUIRED_COLUMNS if col not in col_list])
    col_list.extend([col for col in opt_columns if col not in col_list])

    # Create merged dataframe
    if merge_df_list:
        df_merge = pd.concat(merge_df_list, axis=0, sort=False)
    else:
        df_merge = pd.DataFrame([], columns=col_list)

    # Set column data types
    for col in opt_columns:
        if col in OPT_COL_DTYPE:
            df_merge[col] = df_merge[col].astype(OPT_COL_DTYPE[col])

    # Add discovery class
    if 'DISC_CLASS' not in col_list and 'MERGE_AF' in opt_columns and 'MERGE_AC' in opt_columns:
        if df_merge.shape[0] > 0:
            df_merge['DISC_CLASS'] = get_disc_class(df_merge)
            col_list.append('DISC_CLASS')
        else:
            df_merge = pd.DataFrame([], columns=col_list + ['DISC_CLASS'])

    # Sort
    df_merge.sort_values(['#CHROM', 'POS'], inplace=True)

    # Order columns
    head_cols = ['#CHROM', 'POS', 'END', 'ID']

    if 'SVTYPE' in df_merge.columns:
        head_cols += ['SVTYPE']

    if 'SVLEN' in df_merge.columns:
        head_cols += ['SVLEN']

    tail_cols = [col for col in col_list if col not in head_cols]

    df_merge = df_merge.loc[:, head_cols + tail_cols]

    df_merge.reset_index(drop=True, inplace=True)

    return df_merge


def get_support_table(df, df_next, threads, offset_max, ro_szro_min, match_ref, match_alt):

    interval_flank = (offset_max if offset_max is not None else 0) + 1

    if df_next.shape[0] > 0:

        df_support_list_chrom = list()

        chrom_list = sorted(set(df['#CHROM']) | set(df_next['#CHROM']))

        for chrom in chrom_list:

            df_chrom = df.loc[df['#CHROM'] == chrom]
            df_next_chrom = df_next.loc[df_next['#CHROM'] == chrom]

            # Split merge, isolate to overlapping intervals before merging.
            # This strategy limits combinatorial explosion merging large sets.

            # Create an interval tree of records to intersect.
            #
            # For each interval, the data element will be a tuple of two sets:
            #   [0]: Set of source variant IDs in the interval.
            #   [1]: Set of target variant IDs in the interval.
            tree = intervaltree.IntervalTree()

            # Add max intervals for each source variant
            for row_index, row in df_chrom.iterrows():
                tree.addi(
                    row['POS'] - interval_flank,
                    row['END'] + interval_flank,
                    ({row['ID']}, set())
                )

            # For each target variant in turn, merge all source intervals it intersects. Source and target ID sets
            # in the interval data are merged with the intervals.
            for row_index, row in df_next_chrom.iterrows():

                pos = row['POS']
                end = row['END']

                source_rows = set()
                target_rows = {row['ID']}

                pos_set = set()
                end_set = set()

                # Collapse intersecting intervals
                for interval in tree[pos:end]:
                    pos_set.add(interval.begin)
                    end_set.add(interval.end)

                    source_rows |= interval.data[0]
                    target_rows |= interval.data[1]

                    tree.discard(interval)

                # Add new interval if any
                if source_rows:
                    tree.addi(min(pos_set), max(end_set), (source_rows, target_rows))

            # Create a list of tuples where each element is the data from the interval. Include only intervals
            # with at least one target row.
            record_pair_list = [interval.data for interval in tree if len(interval.data[1]) > 0]

            del tree

            # Report
            print('\t* Split ref {} into {} parts'.format(chrom, len(record_pair_list)))

            # Shortcut if no records
            if len(record_pair_list) > 0:

                # Init merged table list (one for each record pair)
                df_support_list = [None] * len(record_pair_list)

                # Setup jobs
                pool = multiprocessing.Pool(threads)

                kwd_args = {
                    'szro_min': ro_szro_min,
                    'offset_max': offset_max,
                    'priority': ['RO', 'OFFSET', 'SZRO'],
                    'threads': 1,
                    'match_ref': match_ref,
                    'match_alt': match_alt
                }

                # Setup callback handler
                def _apply_parallel_cb_result(record_pair_index, df_support_list):
                    """ Get a function to save results. """

                    def callback_handler(subdf):
                        df_support_list[record_pair_index] = subdf

                    return callback_handler

                def _apply_parallel_cb_error(record_pair_index, df_support_list):
                    """Get an error callback function"""

                    def callback_handler(ex):
                        df_support_list[record_pair_index] = ex

                        pool.terminate()

                        print(traceback.format_exc())

                    return callback_handler

                # Submit jobs
                for record_pair_index in range(len(record_pair_list)):
                    pool.apply_async(
                        analib.variant.nearest_by_svlen_overlap,
                        (
                            df_chrom.loc[
                                df_chrom['ID'].apply(lambda var_id: var_id in record_pair_list[record_pair_index][0])
                            ],
                            df_next_chrom.loc[
                                df_next_chrom['ID'].apply(lambda var_id: var_id in record_pair_list[record_pair_index][1])
                            ]
                        ),
                        kwd_args,
                        _apply_parallel_cb_result(record_pair_index, df_support_list),
                        _apply_parallel_cb_error(record_pair_index, df_support_list)
                    )

                # Wait for jobs
                pool.close()
                pool.join()

                # Check for exceptions
                for df_support in df_support_list:
                    if issubclass(df_support.__class__, Exception):
                        raise df_support

                # Check for null output
                n_fail = np.sum([val is None for val in df_support_list])

                if n_fail > 0:
                    raise RuntimeError('Failed merging {} of {} record groups'.format(n_fail, len(df_support_list)))

                # Merge supporting dataframes
                df_support = pd.concat(df_support_list, axis=0, sort=False).reset_index(drop=True)

                # Clean up
                del df_support_list

            else:
                df_support = pd.DataFrame(columns=['ID', 'TARGET_ID', 'OFFSET', 'RO', 'SZRO', 'OFFSZ'])

            # Add to list (one list item per chromosome)
            df_support_list_chrom.append(df_support)

            # Clean up
            del record_pair_list

        # Merge chromosomes and return
        df_support = pd.concat(df_support_list_chrom, axis=0)
        del df_support_list_chrom

    else:
        df_support = analib.variant.nearest_by_svlen_overlap(
            df, df_next,
            ro_szro_min,
            offset_max,
            priority=['RO', 'OFFSET', 'SZRO'],
            threads=threads,
            match_ref=match_ref,
            match_alt=match_alt
        )

    return df_support


def get_support_table_nrid(df, df_next):
    """
    Get an intersect table of exact-match variants by ID.

    :param df: Dataframe.
    :param df_next: Next dataframe.

    :return: Return a support table with exact matches.
    """

    intersect_set = set(df['ID']) & (set(df_next['ID']))

    intersect_list = [id for id in df['ID'] if id in intersect_set]

    intersect_n = len(intersect_list)

    return pd.DataFrame(
        zip(intersect_list, intersect_list, [0] * intersect_n, [1] * intersect_n, [1] * intersect_n, [0] * intersect_n),
        columns=['ID', 'TARGET_ID', 'OFFSET', 'RO', 'SZRO', 'OFFSZ']
    )
