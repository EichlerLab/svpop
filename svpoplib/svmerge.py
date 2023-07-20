"""
Code for merging variant sets.
"""

import collections
import intervaltree
import multiprocessing
import numpy as np
import pandas as pd
import re
import sys
import traceback

import svpoplib
import kanapy


def get_merge_def(def_name, config, default_none=False):
    """
    Get a merge definition string from a configured alias, `def_name`. Returns
    `config['merge_def'][def_name]` if it exists and `config` is not `None`, and returns `None` otherwise.

    :param def_name: Definition name to search. Name must be alpha-numeric and may contain dashes or underscores.
    :param config: Configuration dictionary (or `None`).
    :param default_none: Return `None` if there is no definition alias (default is to return `def_name`).

    :return: Configuration definition or `def_name` if not found.
    """

    # Do not check merge_def if config has no merge_def section or if def_name contains illegal characters
    if config is None or 'merge_def' not in config:
        if default_none:
            return None

        return def_name

    # Sanity checks
    if def_name in config['merge_def']:
        if def_name in {'nr', 'nrsnp', 'nrsnv'}:
            raise RuntimeError(f'Cannot redefine built-in merge strategy: {def_name}')

        if not re.match('^[a-zA-Z0-9\-_]+$', def_name):
            raise RuntimeError(f'Pre-defined merge strategy key contains illegal characters (only allows alhpa, numeric, underscore, and dash): {def_name}')

    # Get pre-defined merge config
    return config['merge_def'].get(def_name, None if default_none else def_name)


def merge_variants(bed_list, sample_names, strategy, fa_list=None, subset_chrom=None, threads=1):
    """
    Merge variants from multiple samples.

    :param bed_list: List of BED files to merge where each BED file is from one samples.
    :param sample_names: List of samples names. Each element matches an element at the same location in
        `bed_list`.
    :param strategy: Describes how to merge variants.
    :param fa_list: List of FASTA files (variant ID is FASTA record ID). Needed if sequences are used during merging.
    :param subset_chrom: Merge only records from this chromosome. If `None`, merge all records. May be a string (single
        chromosome name) or a list, tuple, or set of chromosome names.
    :param threads: Number of threads to use for intersecting variants.

    :return: A Pandas dataframe of a BED file with the index set to the ID column.
    """

    # Check input
    if len(bed_list) != len(sample_names):
        raise RuntimeError('Sample name list length ({}) does not match the input file list length ({})'.format(
            len(sample_names), len(bed_list))
        )

    if len(bed_list) == 0:
        raise RuntimeError('Cannot merge 0 samples')

    # Parse parameters
    merge_config = svpoplib.svmergeconfig.params.get_merge_config(strategy)

    # Run
    if merge_config.strategy in {'nr', 'nrsnv', 'nrsnp'}:
        return merge_variants_nr(
            bed_list, sample_names, merge_config,
            fa_list=fa_list, subset_chrom=subset_chrom, threads=threads
        )
    else:
        raise RuntimeError('Unrecognized merge strategy "{}": {}'.format(merge_config.strategy, strategy))


def merge_variants_nr(bed_list, sample_names, merge_config, fa_list=None, subset_chrom=None, threads=1, verbose=False):
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
        2) RO: If ro and/or szro is set, then a strict reciprocal-overlap is run. The minimum overlap is ro if it was
            defined and szro if ro was not defined.
        3) Size-overlap + offset: Try intersecting by reciprocal-size-overlap with a maximum offset.

    :param bed_list: List of BED files to merge where each BED file is from one samples.
    :param sample_names: List of samples names. Each element matches an element at the same location in
        `bed_list`.
    :param merge_config: Parameter object controlling this merge.
    :param fa_list: List of FASTA files matching `bed_list` and `sample_names`. FASTA files contain sequences for
        variants where each FASTA record ID is the variant ID and the sequence is the variant sequence.
    :param subset_chrom: Merge only records from this chromosome. If `None`, merge all records. May be a string (single
        chromosome name) or a list, tuple, or set of chromosome names.
    :param threads: Number of threads to use for intersecting variants.
    :param verbose: Print progress if True.

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

    if any([',' in val or ';' in val for val in sample_names]):
        raise RuntimeError('Error: Sample names contain commas or semicolons')

    n_samples = len(sample_names)

    if len(set(sample_names)) != n_samples:
        raise RuntimeError('Sample names may not be duplicated')

    # Report parameters
    if verbose:
        print(merge_config.__repr__(pretty=True))

    # Check fa_list if sequences are required
    seq_in_col = False

    if merge_config.read_seq:
        if fa_list is None:
            seq_in_col = True
            fa_list = [None] * n_samples

        elif len(fa_list) != n_samples:
            n_fa = len(fa_list)

            raise RuntimeError(f'Non-redundant merge requires variant sequences, but fa_list is not the same length as the number of samples: expected {n_samples}, fa_list length {n_fa}')

    else:
        fa_list = [None] * n_samples

    # Set required columns for variant DataFrames
    col_list = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN']

    if merge_config.refalt:
        col_list += ['REF']
        col_list += ['ALT']

    if merge_config.read_seq:
        col_list += ['SEQ']


    # Note:
    # The first time a variant is found, it is added to the merge set (df), and it is supported by itself. The overlap
    # score is defined as -1. As variants are added that intersect with a variant already in the merge set, SAMPLE and
    # ID are copied from that variant, the supporting variant is noted in SUPPORT_ID and SUPPORT_SAMPLE, and the overlap
    # score is between 0 and 1 (depending on how well they overlap).

    # Initialize a table of variants with the first sample. All variants from this sample are in the merged set.
    sample_name = sample_names[0]

    if verbose:
        print('Merging: {}'.format(sample_name))

    df = read_variant_table(bed_list[0], sample_name, subset_chrom, fa_list[0], col_list)

    if seq_in_col:
        if 'SEQ' not in df.columns:
            raise RuntimeError(f'Merge requires variant sequences, but none were provided through FASTA files or as SEQ columns: {bed_list[0]}')

    # Add tracking columns
    df['SAMPLE'] = sample_name
    df['SUPPORT_SAMPLE'] = sample_name
    df['SUPPORT_ID'] = df['ID']

    # df['SUPPORT_ID'] = df['ID']
    # df['SUPPORT_SAMPLE'] = sample_name

    # df['SUPPORT_OFFSET'] = -1
    # df['SUPPORT_RO'] = -1
    # df['SUPPORT_SZRO'] = -1
    # df['SUPPORT_OFFSZ'] = -1
    # df['SUPPORT_MATCH'] = -1

    # Setup a dictionary to translate support sample table columns to the merged table columns
    support_col_rename = {
        'TARGET_ID': 'SUPPORT_ID',
        'OFFSET': 'SUPPORT_OFFSET',
        'RO': 'SUPPORT_RO',
        'SZRO': 'SUPPORT_SZRO',
        'OFFSZ': 'SUPPORT_OFFSZ',
        'MATCH': 'SUPPORT_MATCH'
    }

    # Add each variant to the table
    df_support_list = list()

    for index in range(1, len(bed_list)):

        ## Read ##

        # Read next variant table
        sample_name = sample_names[index]

        if verbose:
            print('Merging: {}'.format(sample_name))

        df_next = read_variant_table(bed_list[index], sample_name, subset_chrom, fa_list[index], col_list)
        df_next['SAMPLE'] = sample_name

        if seq_in_col:
            if 'SEQ' not in df_next.columns:
                raise RuntimeError(f'Merge requires variant sequences, but none were provided through FASTA files or as SEQ columns: {bed_list[index]}')

        ## Build intersect support table for each intersect test (exact, RO, etc) ##
        support_table_list = list()

        # Create copies of df and df_next that can be subset (variants removed by each phase of intersection)
        df_sub = df.copy()
        df_next_sub = df_next.copy()

        # Process match types
        for merge_spec in merge_config.spec_list:
            if verbose:
                print(f'* {merge_spec}')

            spec_type = merge_spec.spec_type.lower()

            # Get support table
            if spec_type == 'exact':

                # INTERSECT: Exact match
                df_support = get_support_table_exact(
                    df_sub, df_next_sub,
                    merge_spec.align_match_prop, merge_spec.aligner,
                    merge_spec.refalt, merge_spec.refalt
                )

            elif merge_spec.spec_type == 'ro':

                # INTERSECT: RO (Reciprocal overlap)
                df_support = get_support_table(
                    df=df_sub,
                    df_next=df_next_sub,
                    threads=threads,
                    ro_min=merge_spec.ro,
                    offset_max=merge_spec.dist,
                    match_ref=merge_spec.refalt,
                    match_alt=merge_spec.refalt,
                    align_match_prop=merge_spec.align_match_prop,
                    aligner=merge_spec.aligner,
                    verbose=verbose
                )

            elif merge_spec.spec_type == 'szro':

                # INTERSECT: SZRO (Size reciprocal overlap)
                df_support = get_support_table(
                    df=df_sub,
                    df_next=df_next_sub,
                    threads=threads,
                    szro_min=merge_spec.szro,
                    offset_max=merge_spec.dist,
                    offsz_max=merge_spec.szdist,
                    match_ref=merge_spec.refalt,
                    match_alt=merge_spec.refalt,
                    align_match_prop=merge_spec.align_match_prop,
                    aligner=merge_spec.aligner,
                    verbose=verbose
                )

                pass

            elif merge_spec.spec_type == 'distance':

                # INTERSECT: Distance
                df_support = get_support_table(
                    df=df_sub,
                    df_next=df_next_sub,
                    threads=threads,
                    szro_min=merge_spec.szro,
                    offset_max=merge_spec.dist,
                    offsz_max=merge_spec.szdist,
                    match_ref=merge_spec.refalt,
                    match_alt=merge_spec.refalt,
                    align_match_prop=merge_spec.align_match_prop,
                    aligner=merge_spec.aligner
                )


            else:
                raise RuntimeError(f'Unknown merge specification type: {spec_type}')

            # Append support variants and subset (do not match already matched variants)
            if df_support is not None:
                support_table_list.append(df_support)

                id_set = set(df_support['ID'])
                id_next_set = set(df_support['TARGET_ID'])

                df_sub = df_sub.loc[df_sub['ID'].apply(lambda val: val not in id_set)]
                df_next_sub = df_next_sub.loc[df_next_sub['ID'].apply(lambda val: val not in id_next_set)]

        # Clean
        del(df_sub)
        del(df_next_sub)
        del(df_support)

        # Construct support table
        #
        # Table describes each variant
        #
        # Support table columns:
        # * #CHROM, POS, END, ID, SVTYPE, SVLEN: Columns of the final merged table
        # * SAMPLE:

        if len(support_table_list) > 0 and any([df_support_check.shape[0] > 0 for df_support_check in support_table_list]):
            df_support = pd.concat(support_table_list)
        else:
            df_support = pd.DataFrame([], columns=['ID', 'SUPPORT_ID'])  # Only the ID column is read if the DataFrame is empty

        if df_support.shape[0] > 0:

            # Set Index
            df_support.set_index('ID', inplace=True, drop=False)
            df_support.index.name = 'INDEX'

            # Add sample name
            df_support['SUPPORT_SAMPLE'] = sample_name
            df_support['SAMPLE'] = df['SAMPLE']

            # Rename columns (offset merger column names to support column names)
            df_support.columns = [support_col_rename.get(col, col) for col in df_support.columns]

            # Remove redundant support - multiple sources of support for one variant
            # Should not occur, but drop if it does
            df_support = df_support.sort_values(
                ['SUPPORT_RO', 'SUPPORT_OFFSET', 'SUPPORT_SZRO', 'SUPPORT_MATCH'],
                ascending=[False, True, False, False]
            ).drop_duplicates('ID', keep='first')

            # Arrange columns
            # df_support['#CHROM'] = df_next['#CHROM']
            # df_support['POS'] = df_next['POS']
            # df_support['END'] = df_next['END']
            # df_support['SVTYPE'] = df_next['SVTYPE']
            # df_support['SVLEN'] = df_next['SVLEN']
            #
            # if 'REF' in df_next.columns:
            #     df_support['REF'] = df_next['REF']
            #
            # if 'ALT' in df_next.columns:
            #     df_support['ALT'] = df_next['ALT']

            #df_support = df_support.loc[:, [col for col in df.columns if col != 'SEQ']]

            df_support_list.append(df_support)

        # Read new variants from this sample (variants that do not support an existing call)
        support_id_set = set(df_support['SUPPORT_ID'])

        df_new = df_next.loc[df_next['ID'].apply(lambda val: val not in support_id_set)].copy()

        if df_new.shape[0] > 0:
            df_new['SAMPLE'] = sample_name
            df_new['SUPPORT_SAMPLE'] = sample_name
            df_new['SUPPORT_ID'] = df_new['ID']

            # De-duplicate IDs
            df_new['ID'] = svpoplib.variant.version_id(df_new['ID'], set(df['ID']))
            df_new.set_index('ID', inplace=True, drop=False)

            # Ensure consistent columns
            if set(df.columns) != set(df_new.columns):
                raise RuntimeError('Error merging variants for sample "{}": Columns mismatch in new variant set: Existing = "{}", New = "{}"'.format(
                    sample_name,
                    ', '.join([col for col in df.columns if col not in set(df_new.columns)]),
                    ', '.join([col for col in df_new.columns if col not in set(df.columns)])
                ))

            # Append new variants
            df = pd.concat([df, df_new.loc[:, df.columns]], axis=0)
            df.sort_values(['#CHROM', 'POS', 'END', 'ID'], inplace=True)

    # Remove SEQ
    if 'SEQ' in df.columns:
        del(df['SEQ'])

    # Finalize merged variant set
    if df.shape[0] > 0:

        # Reformat df as a support table
        df['SUPPORT_OFFSET'] = 0.0
        df['SUPPORT_RO'] = 1.0
        df['SUPPORT_SZRO'] = 1.0
        df['SUPPORT_OFFSZ'] = 0.0
        df['IS_PRIMARY'] = True

        # Choose a MATCH default for self variants (1.0 if any variants were matched by sequence, np.nan otherwise)
        if merge_config.any_match():
            df['SUPPORT_MATCH'] = 1.0
        else:
            df['SUPPORT_MATCH'] = np.nan

        df = df[[
            'ID', 'SAMPLE',
            'SUPPORT_ID', 'SUPPORT_SAMPLE',
            'SUPPORT_OFFSET', 'SUPPORT_RO', 'SUPPORT_SZRO', 'SUPPORT_OFFSZ', 'SUPPORT_MATCH',
            'IS_PRIMARY'
        ]]

        # Add support variants
        if len(df_support_list) > 0:
            df_support = pd.concat(df_support_list)

            df_support['IS_PRIMARY'] = False

            df_support = df_support[[
                'ID', 'SAMPLE',
                'SUPPORT_ID', 'SUPPORT_SAMPLE',
                'SUPPORT_OFFSET', 'SUPPORT_RO', 'SUPPORT_SZRO', 'SUPPORT_OFFSZ', 'SUPPORT_MATCH',
                'IS_PRIMARY'
            ]]

            df = pd.concat([df, df_support])

        # Make SAMPLE and SUPPORT_SAMPLE categorical (sort in the same order as they were merged)
        df['SAMPLE'] = pd.Categorical(df['SAMPLE'], sample_names)
        df['SUPPORT_SAMPLE'] = pd.Categorical(df['SUPPORT_SAMPLE'], sample_names)

        # Sort by support (best support first)
        df['SUPPORT_OFFSET'] = df['SUPPORT_OFFSET'].apply(lambda val: np.max([0, val]))
        df['SUPPORT_RO'] = np.abs(df['SUPPORT_RO'])
        df['SUPPORT_SZRO'] = np.abs(df['SUPPORT_SZRO'])
        df['SUPPORT_OFFSZ'] = np.abs(df['SUPPORT_OFFSZ'])
        df['SUPPORT_MATCH'] = np.abs(df['SUPPORT_MATCH'])

        # Sort all support
        df.sort_values(
            ['IS_PRIMARY', 'SAMPLE', 'SUPPORT_RO', 'SUPPORT_OFFSET', 'SUPPORT_SZRO', 'SUPPORT_OFFSZ', 'SUPPORT_MATCH'],
            ascending=[False, True, False, True, False, False, False],
            inplace=True
        )

        # Find best support variant for each mergeset variant (per sample)
        # Duplicates should not occur, but drop them if they do
        df.drop_duplicates(['ID', 'SUPPORT_SAMPLE'], keep='first', inplace=True)

        # Re-sort by ID then SAMPLE (for concatenating stats in order)
        df.sort_values(['ID', 'SUPPORT_SAMPLE'], inplace=True)

        df_support = df.groupby(
            'ID'
        ).apply(lambda subdf: pd.Series(
            [
                ','.join(subdf['SUPPORT_SAMPLE']),
                ','.join(subdf['SUPPORT_ID']),

                ','.join(['{:.2f}'.format(val) for val in subdf['SUPPORT_RO']]),
                ','.join(['{:.0f}'.format(val) for val in subdf['SUPPORT_OFFSET']]),
                ','.join(['{:.2f}'.format(val) for val in subdf['SUPPORT_SZRO']]),
                ','.join(['{:.2f}'.format(val) for val in subdf['SUPPORT_OFFSZ']]),
                ','.join(['{:.2f}'.format(val) for val in subdf['SUPPORT_MATCH']]),
                subdf.shape[0]
            ],
            index=[
                'MERGE_SAMPLES', 'MERGE_VARIANTS',
                'MERGE_RO', 'MERGE_OFFSET', 'MERGE_SZRO', 'MERGE_OFFSZ', 'MERGE_MATCH',
                'MERGE_N'
            ]
        ))

        if not merge_config.any_match():
            df_support['SUPPORT_MATCH'] = np.nan

        df_support.reset_index(inplace=True, drop=False)

    else:
        df_support = pd.DataFrame(
            [],
            columns=[
                'ID',
                'MERGE_SAMPLES', 'MERGE_VARIANTS',
                'MERGE_RO', 'MERGE_OFFSET', 'MERGE_SZRO', 'MERGE_OFFSZ', 'MERGE_MATCH',
                'MERGE_N'
            ]
        )

    # Merge original BED files by by df_support
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


def merge_sample_by_support(df_support, bed_list, sample_names):
    """
    Take a support table and generate the merged variant BED file.
    
    The support table must have at least columns:
    1) ID: ID in the final merged callset. Might have been altered to avoid name clashes (e.g. "XXX.1").
    2) MERGE_SRC: Sample variant should be extracted from
    3) MERGE_SRC_ID: ID of of the variant to extract from the source. This becomes the record that represents all
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

    REQUIRED_COLUMNS = ['ID', 'MERGE_SAMPLES', 'MERGE_VARIANTS']

    OPT_COL = ['MERGE_RO', 'MERGE_OFFSET', 'MERGE_SZRO', 'MERGE_OFFSZ', 'MERGE_MATCH', 'MERGE_SRC', 'MERGE_SRC_ID']

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

    # Merged variant IDs should be unique
    svpoplib.variant.check_unique_ids(df_support, 'Post-merge error')

    # Set Merge source
    if 'MERGE_SRC' not in df_support.columns:
        df_support['MERGE_SRC'] = df_support['MERGE_SAMPLES'].apply(lambda val_list: val_list.split(',')[0])

    if 'MERGE_SRC_ID' not in df_support.columns:
        df_support['MERGE_SRC_ID'] = df_support['MERGE_VARIANTS'].apply(lambda val_list: val_list.split(',')[0])

    # Get optional columns that are present
    # Include MERGE_SRC and MERGE_SRC_ID if they were generated
    opt_columns = [opt_col for opt_col in OPT_COL if opt_col in df_support.columns]

    # Merge original SVs
    col_list = list()  # List of columns in the final output
    merge_df_list = list()  # List of dataframes to be merged into final output

    for index in range(len(bed_list)):
        bed_file_name = bed_list[index]
        sample_name = sample_names[index]

        # Get merged variants for this sample
        df_support_sample = df_support.loc[df_support['MERGE_SRC'] == sample_name].copy()

        # Set index to original ID in this sample
        df_support_sample.set_index('MERGE_SRC_ID', inplace=True, drop=False)
        df_support_sample.index.name = 'INDEX'

        # Read variants from sample
        if type(bed_file_name) != pd.DataFrame:
            df_sample = pd.read_csv(bed_file_name, sep='\t', header=0)  # TODO: low_memory=False
        else:
            df_sample = bed_file_name.copy()

        # Update columns
        col_list.extend([col for col in df_sample.columns if col not in col_list])

        # Get variant names
        if df_support_sample.shape[0] == 0:
            continue

        df_sample.set_index('ID', inplace=True, drop=True)
        df_sample.index.name = 'INDEX'

        df_sample = df_sample.loc[list(df_support_sample['MERGE_SRC_ID'])]

        # Transfer required columns
        df_sample['ID'] = df_support_sample['ID']  # ID in the final merged callset
        df_sample['MERGE_SAMPLES'] = df_support_sample['MERGE_SAMPLES']
        df_sample['MERGE_VARIANTS'] = df_support_sample['MERGE_VARIANTS']
        df_sample['MERGE_SRC'] = df_support_sample['MERGE_SRC']
        df_sample['MERGE_SRC_ID'] = df_support_sample['MERGE_SRC_ID']

        # Transfer optional columns annotations
        for col in opt_columns:
            df_sample[col] = df_support_sample[col]

        # Append to list
        merge_df_list.append(df_sample)

    # Append required and optional columns to the end
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

    # Sort
    df_merge.sort_values(['#CHROM', 'POS', 'END', 'ID'], inplace=True)

    # Get column order
    head_cols = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE']

    if 'SVLEN' in df_merge.columns:
        head_cols += ['SVLEN']

    if 'REF' in df_merge.columns:
        head_cols += ['REF']

    if 'ALT' in df_merge.columns:
        head_cols += ['ALT']

    # Add missing columns (unique fields from a caller with no records, keeps it consistent with the whole callset)
    for col in col_list:
        if col not in df_merge.columns:
            if col in head_cols:
                raise RuntimeError(f'Missing head column while merging sample support (post-merge step to prepare final table): {col}')

            df_merge[col] = np.nan

    # Order columns
    df_merge = svpoplib.variant.order_variant_columns(df_merge, head_cols, ['MERGE_SAMPLES', 'MERGE_VARIANTS'] + opt_columns)

    # Reset index
    df_merge.reset_index(drop=True, inplace=True)

    # Return merged variants
    return df_merge


def read_variant_table(
        bed_file_name,
        sample_name,
        subset_chrom=None,
        fa_file_name=None,
        col_list=('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN')
    ):
    """
    Read a DataFrame of variants and prepare for merging.

    :param bed_file_name: BED file name.
    :param sample_name: Sample name.
    :param subset_chrom: Subset to this chromosome (or None to read all). May be a string (single chromosome name) or
        a list, tuple, or set of chromosome names.
    :param fa_file_name: FASTA file name to read variant sequences (into SEQ column), or `None` if variant sequences
        do not need to be read. FASTA record IDs must match the variant ID column.
    :param col_list: List of columns to be read. Used to order and filter columns.

    :return: Prepared variant DataFrame.
    """


    # Ensure column list contains required columns
    head_cols = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN']
    col_list = head_cols + [col for col in col_list if col not in head_cols]

    # Read variants
    col_set = set(col_list)

    if type(bed_file_name) != pd.DataFrame:
        df = svpoplib.pd.read_csv_chrom(
            bed_file_name, chrom=subset_chrom,
            sep='\t', header=0,
            usecols=lambda col: col in col_set,
            dtype={'#CHROM': str}
        )

    else:
        df = bed_file_name.copy()

        missing_cols = sorted({col for col in col_set if col not in df.columns})

        if missing_cols:
            raise RuntimeError(f'Error getting variant table from an existing DataFrame: Missing columns: {", ".join(missing_cols)}')

        df['#CHROM'] = df['#CHROM'].astype(str)

    # Read SEQ column
    if fa_file_name is not None:
        if 'SEQ' in df.columns:
            raise RuntimeError(f'Duplicate SEQ sources for BED file "{bed_file_name}": BED contains a SEQ column, and read_variant_table() SEQ from a FASTA file')

        df = df.join(svpoplib.seq.fa_to_series(fa_file_name), on='ID', how='left')

        if np.any(pd.isnull(df['SEQ'])):
            id_missing = list(df.loc[pd.isnull(df['SEQ']), 'ID'])

            raise RuntimeError('Missing {} records in FASTA for sample {}: {}{}'.format(
                len(id_missing), sample_name, ', '.join(id_missing[:3]), '...' if len(id_missing) > 3 else ''
            ))

        if 'SEQ' not in col_list:
            col_list = tuple(list(col_list) + ['SEQ'])

    elif 'SEQ' in df.columns:
        if 'SEQ' not in col_list:
            col_list = tuple(list(col_list) + ['SEQ'])

    # Set defaults for missing columns
    if 'SVLEN' not in df.columns:

        if 'END' not in df.columns or 'SVTYPE' not in df.columns:
            raise RuntimeError(f'Missing SVLEN in {bed_file_name}: Need both SVTYPE and END to set automatically')

        if (np.any(df['SVTYPE'].apply(lambda val: val.upper()) == 'INS')):
            raise RuntimeError(f'Missing SVLEN in {bed_file_name}: Cannot compute for insertions (SVTYPE must not be INS for any record)')

        df['SVLEN'] = df['END'] - df['POS']

    if 'SVTYPE' not in df.columns:
        df['SVTYPE'] = 'RGN'

    # Order and subset columns
    try:
        df = svpoplib.variant.order_variant_columns(df, col_list, subset=True)
    except RuntimeError as ex:
        raise RuntimeError(f'Error checking columns in {bed_file_name}: {ex}')

    # Check SVLEN
    if np.any(df['SVLEN'] < 0):
        raise RuntimeError(f'Negative SVLEN entries in {bed_file_name}')

    # Sort
    df.sort_values(['#CHROM', 'POS', 'END', 'ID'], inplace=True)

    # Check for unique IDs
    svpoplib.variant.check_unique_ids(df, 'Error reading merge input')

    # Set index
    df.set_index(df['ID'], inplace=True, drop=False)
    df.index.name = 'INDEX'

    # Return variant DataFrame
    return df


def get_support_table(
        df, df_next,
        threads=1,
        ro_min=None,
        szro_min=None,
        offset_max=None,
        offsz_max=None,
        match_ref=False,
        match_alt=False,
        align_match_prop=None,
        aligner=None,
        verbose=False
    ):
    """
    Get a table describing matched variants between `df` and `df_next` and columns of evidence for support.

    :param df: Set of variants in the accepted set.
    :param df_next: Set of variants to intersect with `df`.
    :param threads: Number of threads to run.
    :param ro_min: Minimum size reciprocal overlap (computed on size only independent of position).
    :param szro_min: Minimum reciprocal-overlap.
    :param offset_max: Max breakpoint offset (offset is the minimum of start position and end position offsets).
    :param offsz_max: Maximum offset/svlen proportion. Computed with the minimum length of the two variants.
    :param match_ref: REF column must match if True.
    :param match_alt: ALT column must match if True.
    :param align_match_prop: Minimum matched base proportion in alignment.
    :param aligner: Configured aligner for matching sequences.
    :param verbose: Print progress if True.

    See `svpoplib.svlenoverlap.nearest_by_svlen_overlap` for a description of the returned columns.

    :return: A table describing variant support between `df` (ID) and `df_next` (TARGET_ID) and supporting evidence
        (OFFSET, RO, SZRO, OFFSZ, MATCH).
    """

    # Determine if the merge can be split into clusters. If True (most merges), then variants are grouped into
    # clusters of events that could possibly merge given location restrictions based on ro_min, offset_max, and
    # offsz_max. szro_min alone is not affected by placement and cannot be used for clustering. Un-clustered merges are
    # degenerate, but could be useful for prioritizing based on attributes (i.e. find all possible matches and choose
    # the best one based on size or other parameters even if they are not used to restrict merges).

    cluster_merge = ro_min is not None or offset_max is not None or offsz_max is not None

    # Base interval flank
    base_interval_flank = offset_max + 1 if offset_max is not None else 1

    if df_next.shape[0] > 0:

        df_support_list_chrom = list()

        chrom_list = sorted(set(df['#CHROM']) | set(df_next['#CHROM']))

        for chrom in chrom_list:

            df_chrom = df.loc[df['#CHROM'] == chrom]
            df_next_chrom = df_next.loc[df_next['#CHROM'] == chrom]

            if cluster_merge:
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
                    if offsz_max is not None:
                        interval_flank = max(base_interval_flank, row['SVLEN'] * offsz_max)
                    else:
                        interval_flank = base_interval_flank

                    tree.addi(
                        row['POS'] - interval_flank,
                        (row['END'] if row['SVTYPE'] != 'INS' else row['POS'] + 1) + interval_flank,
                        ({row['ID']}, set())
                    )

                # For each target variant in turn, merge all source intervals it intersects. Source and target ID sets
                # in the interval data are merged with the intervals.
                for row_index, row in df_next_chrom.iterrows():

                    pos = row['POS']
                    end = (
                        row['END'] if (
                            row['SVTYPE'] != 'INS' or ro_min is None
                        ) else row['POS'] + row['SVLEN']
                    )

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

            else:
                record_pair_list = [
                    (set(df_chrom['ID']), set(df_next_chrom['ID']))
                ]

            # Report
            if verbose:
                print('\t* Split ref {} into {} parts'.format(chrom, len(record_pair_list)))

            # Shortcut if no records
            if len(record_pair_list) > 0:

                # Init merged table list (one for each record pair)
                df_support_list = [None] * len(record_pair_list)

                # Setup jobs
                pool = multiprocessing.Pool(threads)

                kwd_args = {
                    'ro_min': ro_min,
                    'szro_min': szro_min,
                    'offset_max': offset_max,
                    'offsz_max': offsz_max,
                    'priority': ['RO', 'SZRO', 'OFFSET', 'OFFSZ', 'MATCH'],
                    'threads': 1,
                    'match_ref': match_ref,
                    'match_alt': match_alt,
                    'aligner': aligner,
                    'align_match_prop': align_match_prop
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

                        print(f'Failed {record_pair_index}: {ex}', file=sys.stderr)
                        traceback.print_tb(ex.__traceback__)
                        sys.stderr.flush()

                        try:
                            print(f'Terminating: {record_pair_index}', file=sys.stderr)
                            sys.stderr.flush()

                            pool.terminate()

                        except Exception as ex:
                            print(f'Caught error while terminating: {record_pair_index}: {ex}', file=sys.stderr)
                            sys.stderr.flush()

                        print(f'Exiting error handler: {record_pair_index}')

                    return callback_handler

                # Submit jobs
                for record_pair_index in range(len(record_pair_list)):

                    try:
                        pool.apply_async(
                            svpoplib.svlenoverlap.nearest_by_svlen_overlap,
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

                    except:
                        pass

                # Wait for jobs
                if verbose:
                    print('Waiting...')

                sys.stderr.flush()
                sys.stdout.flush()

                pool.close()
                pool.join()
                sys.stderr.flush()
                sys.stdout.flush()

                if verbose:
                    print('Done Waiting.')

                sys.stderr.flush()
                sys.stdout.flush()

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
                df_support = pd.DataFrame(columns=['ID', 'TARGET_ID', 'OFFSET', 'RO', 'SZRO', 'OFFSZ', 'MATCH'])

            # Add to list (one list item per chromosome)
            df_support_list_chrom.append(df_support)

            # Clean up
            del record_pair_list

        # Merge chromosomes and return
        df_support = pd.concat(df_support_list_chrom, axis=0)
        del df_support_list_chrom

    else:
        df_support = svpoplib.svlenoverlap.nearest_by_svlen_overlap(
            df_source=df,
            df_target=df_next,
            ro_min=ro_min,
            szro_min=szro_min,
            offset_max=offset_max,
            offsz_max=offsz_max,
            priority=['RO', 'SZRO', 'OFFSET', 'OFFSZ', 'MATCH'],
            threads=threads,
            match_ref=match_ref,
            match_alt=match_alt
        )

    return df_support if df_support.shape[0] > 0 else None


def get_support_table_exact(df, df_next, align_match_prop=None, aligner=None, match_ref=None, match_alt=None):
    """
    Get an intersect table of exact breakpoint matches.

    :param df: Dataframe.
    :param df_next: Next dataframe.
    :param align_match_prop: Proportion of the max alignment score for matching sequences or None if sequences
        (SEQ column in df and df_next) should not be matched.
    :param aligner: Alignment tool for comparing sequences (svpoplib.aligner.ScoreAligner).
    :param match_ref: Match REF if True (SNVs).
    :param match_alt: Match ALT if True (SNVs).

    :return: Return a support table with exact matches.
    """

    # Set columns to sort by
    sort_cols = ['#CHROM', 'POS', 'SVLEN']

    # Set match_ref
    if match_ref is None:
        match_ref = 'REF' in df.columns or 'REF' in df_next.columns

    if match_ref:
        if 'REF' not in df.columns:
            raise RuntimeError('Cannot match REF for exact match intersect: No REF column in dataframe (df)')

        if 'REF' not in df_next.columns:
            raise RuntimeError('Cannot match REF for exact match intersect: No REF column in dataframe (df_next)')

        sort_cols += ['REF']

    # Set match_alt
    if match_alt is None:
        match_alt = 'ALT' in df.columns or 'ALT' in df_next.columns

    if match_alt:
        if 'ALT' not in df.columns:
            raise RuntimeError('Cannot match ALT for exact match intersect: No ALT column in dataframe (df)')

        if 'ALT' not in df_next.columns:
            raise RuntimeError('Cannot match ALT for exact match intersect: No ALT column in dataframe (df_next)')

        sort_cols += ['ALT']

    # Set match_seq
    match_seq = align_match_prop is not None

    if match_seq:
        if 'SEQ' not in df.columns:
            raise RuntimeError('Cannot match sequences for exact match intersect: No SEQ column in dataframe (df)')

        if 'SEQ' not in df_next.columns:
            raise RuntimeError('Cannot match sequences for exact match intersect: No SEQ column in dataframe (df_next)')

        if aligner is None:
            aligner = svpoplib.aligner.ScoreAligner()

        sort_cols += ['SEQ']

    max_match = np.nan  # Initialize - set during matches if match_seq, left as np.nan for the support table otherwise

    # Check for missing columns
    missing_1 = [col for col in sort_cols if col not in df.columns]
    missing_2 = [col for col in sort_cols if col not in df_next.columns]

    if missing_1 or missing_2:
        raise RuntimeError('Missing columns for exact merging: df="{}", df_next="{}"'.format(
            ', '.join(missing_1), ', '.join(missing_2)
        ))

    # Sort
    df = df.sort_values(sort_cols).copy()
    df_next = df_next.sort_values(sort_cols).copy()

    # Find exact matches
    index_1 = 0
    index_2 = 0

    max_index_1 = df.shape[0]
    max_index_2 = df_next.shape[0]

    df_match_list = list()

    while index_1 < max_index_1 and index_2 < max_index_2:

        # Check by size and position
        cmp_val = is_exact_match_no_seq(df.iloc[index_1], df_next.iloc[index_2], match_ref, match_alt)

        if cmp_val < 0:
            index_2 += 1
            continue

        if cmp_val > 0:
            index_1 += 1
            continue

        # SEQ
        if match_seq:

            # Find all matching rows - search for best align match
            last_index = index_2 + 1

            max_match = aligner.match_prop(df.iloc[index_1]['SEQ'], df_next.iloc[index_2]['SEQ'])
            max_match_index = index_2

            while last_index < max_index_2 and is_exact_match_no_seq(df.iloc[index_1], df_next.iloc[last_index], match_ref, match_alt) == 0:

                match_val = aligner.match_prop(df.iloc[index_1]['SEQ'], df_next.iloc[last_index]['SEQ'])

                if match_val > max_match:
                    max_match = match_val
                    max_match_index = last_index

                last_index += 1

            # No SEQ match for this row in df
            if max_match < align_match_prop:
                index_1 += 1
                continue

            if max_match_index > index_2:
                # Found multiple matches. Swap index_2 and max_match_index
                swap_dict = {
                    index_2: max_match_index,
                    max_match_index: index_2
                }

                df_next = df_next.iloc[[swap_dict.get(i, i) for i in range(max_index_2)]].copy()

        # Found match
        df_match_list.append(pd.Series(
            [
                df.iloc[index_1]['ID'],
                df_next.iloc[index_2]['ID'],
                0, 1, 1, 0, max_match
            ],
            index=['ID', 'TARGET_ID', 'OFFSET', 'RO', 'SZRO', 'OFFSZ', 'MATCH']
        ))

        index_1 += 1
        index_2 += 1

    # Merge match dataframe
    if df_match_list:
        return pd.concat(df_match_list, axis=1).T
    else:
        return None


def is_exact_match_no_seq(row_a, row_b, match_ref, match_alt):
    """
    Determine if two rows are exact matches. This function does not compare sequences, just size and position.

    :param row_a: Row to compare.
    :param row_b: Row to compare.
    :param match_ref: Match REF column. Assumes rows have been checked to verify REF exists.
    :param match_alt: Match ALT column. Assumes rows have been checked to verify ALT exists.

    :return: True if rows are exact matches by size and position (does not evaluate sequence).
    """

    # CHROM
    if row_b['#CHROM'] < row_a['#CHROM']:
        return -1
        # index_2 += 1
        # continue

    elif row_b['#CHROM'] > row_a['#CHROM']:
        return 1
        # index_1 += 1
        # continue

    # POS
    if row_b['POS'] < row_a['POS']:
        return -1

    elif row_b['POS'] > row_a['POS']:
        return 1

    # SVLEN (END)
    if row_b['SVLEN'] < row_a['SVLEN']:
        return -1

    elif row_b['SVLEN'] > row_a['SVLEN']:
        return 1

    # REF
    if match_ref:
        if row_b['REF'] < row_a['REF']:
            return -1

        elif row_b['REF'] > row_a['REF']:
            return 1

    # ALT
    if match_alt:
        if row_b['ALT'] < row_a['ALT']:
            return -1

        elif row_b['ALT'] > row_a['ALT']:
            return 1

    # Match
    return 0
