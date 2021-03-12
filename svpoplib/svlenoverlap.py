"""
Get nearest variant by SVLEN overlap. Used for merging and comparing callsets.
"""

def nearest_by_svlen_overlap(
        df_source, df_target,
        szro_min=None,
        offset_max=None,
        priority=('RO', 'SZRO', 'OFFSET'),
        restrict_samples=False,
        threads=1,
        match_ref=False,
        match_alt=False
):
    """
    For each variant in df_source, get the nearest variant in df_target. Both dataframes must contain fields
    "#CHROM", "POS", "END", and "SVLEN". Return a dataframe with each source variant ("ID") with the best match.

    Match priorities are defined as:
    * RO: Reciprocal overlap
    * OFFSET: Offset (minimum the start position distance or end position distance)
    * SZRO: Reciprocal size overlap. Calculated only on size, offset is not a factor (can be many bp apart).
    * OFFSZ: Offset / size.

    :param df_source: Source dataframe.
    :param df_target: Target dataframe.
    :param szro_min: Reciprocal length proportion of allowed matches.
    :param offset_max: Maximum offset allowed (minimum of start or end postion distance).
    :param priority: A list of match priorities (one or more of OFFSET, RO, SZRO, or OFFSZ, see above).
    :param restrict_samples: If `True` and both dataframes contain a `MERGE_SAMPLES` column, then restrict matches to
        only those that share samples.
    :param threads: Run this many overlap threads in parallel.
    :param match_ref: "REF" column must match between two variants.
    :param match_alt: "ALT" column must match between two variants.

    :return: A dataframe with "ID", "TARGET_ID", "OFFSET", "RO", "SZRO", and "OFFSZ".
    """

    # Return an empty DataFrame if either source or target or empty
    if df_source.shape[0] == 0 or df_target.shape[0] == 0:
        return pd.DataFrame(columns=['ID', 'TARGET_ID', 'OFFSET', 'RO', 'SZRO', 'OFFSZ'])

    # Check for expected columns
    if any(col not in df_source.columns for col in ('#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'ID')):
        raise RuntimeError('Source Dataframe missing at least one column in ("#CHROM", "POS", "END", "SVTYPE", "SVLEN"): {}')

    if any(col not in df_target.columns for col in ('#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'ID')):
        raise RuntimeError('Target Dataframe missing at least one column in ("#CHROM", "POS", "END", "SVTYPE", "SVLEN"): {}')

    # Copy (do not alter original DataFrame)
    df_source = df_source.copy()
    df_target = df_target.copy()

    df_source['SVLEN'] = np.abs(df_source['SVLEN'])
    df_target['SVLEN'] = np.abs(df_target['SVLEN'])

    # IDs must be unique
    if len(set(df_source['ID'])) != df_source.shape[0]:
        raise RuntimeError('Source Dataframe IDs are not unique')

    if len(set(df_target['ID'])) != df_target.shape[0]:
        raise RuntimeError('target Dataframe IDs are not unique')

    # REF and ALT columns must be present if they are used
    if match_ref:
        if 'REF' not in df_source.columns:
            raise RuntimeError('Source table is missing REF column (required when matching reference base)')

        if 'REF' not in df_target.columns:
            raise RuntimeError('Target table is missing REF column (required when matching reference base)')

        df_source['REF'] = df_source['REF'].fillna('').apply(lambda val: val.upper())
        df_target['REF'] = df_target['REF'].fillna('').apply(lambda val: val.upper())

    if match_alt:
        if 'ALT' not in df_source.columns:
            raise RuntimeError('Source table is missing ALT column (required when matching reference base)')

        if 'ALT' not in df_target.columns:
            raise RuntimeError('Target table is missing ALT column (required when matching reference base)')

        df_source['ALT'] = df_source['ALT'].fillna('').apply(lambda val: val.upper())
        df_target['ALT'] = df_target['ALT'].fillna('').apply(lambda val: val.upper())


    # Check priority
    if priority is None:
        raise RuntimeError('Argument "priority" is missing')

    if issubclass(priority.__class__, str):
        priority = [priority]

    elif issubclass(priority.__class__, tuple):
        priority = list(priority)

    elif not issubclass(priority.__class__, list):
        raise RuntimeError('Argument "priority" must be a string, list, or an element name')

    if len(priority) == 0:
        raise RuntimeError('Argument "priority" is empty')

    priority_ascending_cols = {
        'RO': False,      # Highest RO first
        'OFFSET': True,   # Lowest offset first
        'SZRO': False,    # Highest size RO first
        'OFFSZ': True     # Lowest offset / size proportion first
    }

    if any([val not in priority_ascending_cols for val in priority]):
        raise RuntimeError('Argument priority must be a list of set element "{}": {}'.format(
            ', '.join(priority_ascending_cols.keys()), ', '.join(priority)
        ))

    priority_ascending = [priority_ascending_cols[val] for val in priority]

    # Determine if variants are matched on MERGE_SAMPLES
    if restrict_samples and 'MERGE_SAMPLES' in df_source.columns and 'MERGE_SAMPLES' in df_target.columns:
        restrict_samples = True
    else:
        restrict_samples = False

    # Subset and cast to int16 (default int64 uses more memory and is not needed)
    if restrict_samples:
        subset_cols = ['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'ID', 'MERGE_SAMPLES']

    else:
        subset_cols = ['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'ID']

    if match_ref:
        subset_cols += ['REF']

    if match_alt:
        subset_cols += ['ALT']

    df_source = df_source.loc[:, subset_cols]
    df_target = df_target.loc[:, subset_cols]

    # Set column types
    col_types = {
        '#CHROM': np.object,
        'POS': np.int32,
        'END': np.int32,
        'SVLEN': np.int32,
        'ID': np.object
    }

    if restrict_samples:
        col_types['MERGE_SAMPLES'] = np.object

    if df_source.shape[0] > 0:
        df_source = df_source.astype({key: val for key, val in col_types.items() if key in df_source.columns})

    if df_target.shape[0] > 0:
        df_target = df_target.astype({key: val for key, val in col_types.items() if key in df_target.columns})

    # Set index
    df_source.set_index('ID', inplace=True)

    df_target.set_index('ID', inplace=True)
    df_target.index.name = 'TARGET_ID'

    # Set min_len_prop
    if szro_min is not None:
        szro_min = np.float16(szro_min)

        if szro_min <= np.float16(0.0) or szro_min >= np.float16(1.0):
            raise RuntimeError('Length proportion must be between 0 and 1 (inclusive): {}'.format(szro_min))

    # Make set of sample names
    if restrict_samples:
        df_source['MERGE_SAMPLES'] = df_source['MERGE_SAMPLES'].apply(lambda vals: set(vals.split(',')))
        df_target['MERGE_SAMPLES'] = df_target['MERGE_SAMPLES'].apply(lambda vals: set(vals.split(',')))

    # Create an array to save results for parallelized output (one per chrom)
    chrom_list = sorted(set(df_source['#CHROM']))

    df_split_results = [None] * len(chrom_list)
    results_received = {index: False for index in range(len(chrom_list))}

    # Setup arguments
    kwd_args = {
        'df_source': df_source,
        'df_target': df_target,
        'szro_min': szro_min,
        'offset_max': offset_max,
        'restrict_samples': restrict_samples,
        'priority': priority,
        'priority_ascending': priority_ascending,
        'match_ref': match_ref,
        'match_alt': match_alt
    }

    if threads > 1 and len(chrom_list) > 1:

        # Open pool of workers
        pool = multiprocessing.Pool(threads)

        # Submit each split part to the worker pool
        for index in range(len(chrom_list)):
            pool.apply_async(
                _overlap_worker, (chrom_list[index], ), kwd_args,
                _apply_parallel_cb_result(index, df_split_results, results_received),
                _apply_parallel_cb_error(chrom_list[index])
            )

        # Close pool and wait for jobs
        pool.close()
        pool.join()

    else:

        # Fill df_split_results with one thread
        for index in range(len(chrom_list)):
            df_split_results[index] = _overlap_worker(chrom_list[index], **kwd_args)

    # Merge dataframes
    df_match = pd.concat(df_split_results, axis=0)

    df_match = df_match.loc[:, ['ID', 'TARGET_ID', 'OFFSET', 'RO', 'SZRO', 'OFFSZ']]

    # Return dataframe
    return df_match


def _apply_parallel_cb_result(index, df_split_results, results_received):
    """
    Get a function to save results.

    :param index:
    :param df_split_results:
    :return:
    """

    def callback_handler(subdf):
        df_split_results[index] = subdf
        results_received[index] = True

    return callback_handler


def _apply_parallel_cb_error(chrom):

    def callback_handler(subdf):
        raise RuntimeError(
            'Runtime in chromosome {}'.format(chrom)
        ) from subdf

    return callback_handler


def _overlap_worker(
        chrom,
        df_source, df_target,
        szro_min, offset_max,
        restrict_samples,
        priority, priority_ascending,
        match_ref, match_alt
):

    # Get dataframes
    df_source_chr = df_source.loc[df_source['#CHROM'] == chrom].copy()
    df_target_chr = df_target.loc[df_target['#CHROM'] == chrom].copy()

    # Skip if target has no entries for this chromosome
    if df_target_chr.shape[0] == 0 or df_source_chr.shape[0] == 0:
        return pd.DataFrame([], columns=('ID', 'TARGET_ID', 'OFFSET', 'RO', 'SZRO', 'OFFSZ'))

    # Setup list of maximum matches
    overlap_list = list()

    # Get maximum matches
    for index, source_row in df_source_chr.iterrows():

        pos = source_row['POS']
        end = source_row['END']
        svlen = source_row['SVLEN']

        df_target_row_pool = df_target_chr.copy()

        # Filter by REF and ALT
        if match_ref:
            df_target_row_pool = df_target_row_pool.loc[df_target_row_pool['REF'] == source_row['REF']].copy()

        if match_alt:
            df_target_row_pool = df_target_row_pool.loc[df_target_row_pool['ALT'] == source_row['ALT']].copy()

        if df_target_row_pool.shape[0] == 0:
            continue

        # Calculate SZRO (size overlap)
        df_target_row_pool['SZRO'] = df_target_row_pool.apply(
            lambda row: np.min([svlen / row['SVLEN'], row['SVLEN'] / svlen]),
            axis=1
        )

        # Filter by svlen overlap
        if szro_min is not None:
            df_target_row_pool = df_target_row_pool.loc[df_target_row_pool['SZRO'] >= szro_min].copy()

        # Stop if no matches
        if df_target_row_pool.shape[0] == 0:
            continue

        # Calculate offset
        df_target_row_pool['OFFSET'] = df_target_row_pool.apply(
            lambda row: np.min(
                [
                    np.abs(pos - row['POS']),
                    np.abs(end - row['END'])
                ]
            ), axis=1
        )

        # Filter by offset
        if offset_max is not None:
            df_target_row_pool = df_target_row_pool.loc[df_target_row_pool['OFFSET'] <= offset_max].copy()

        # Filter by samples
        if restrict_samples:
            df_target_row_pool = df_target_row_pool.loc[df_target_row_pool['MERGE_SAMPLES'].apply(
                lambda sample_set: bool(sample_set & source_row['MERGE_SAMPLES'])
            )]

        # Stop if no matches
        if df_target_row_pool.shape[0] == 0:
            continue

        # Redefine end points for insertions (make them by-length for reciprocal overlap calculation)
        if source_row['SVTYPE'] == 'INS':
            df_target_row_pool['END'] = df_target_row_pool['POS'] + df_target_row_pool['SVLEN']
            end = pos + svlen

        # Get distance calculations
        df_distance = df_target_row_pool.apply(lambda row: pd.Series(
            [
                row['OFFSET'],
                reciprocal_overlap(pos, end, row['POS'], row['END']),
                row['SZRO'],
                row['OFFSET'] / svlen
            ],
            index=['OFFSET', 'RO', 'SZRO', 'OFFSZ']
        ), axis=1)

        max_row = df_distance.reset_index().sort_values(priority, ascending=priority_ascending).iloc[0]

        max_row['ID'] = index

        # Save match record
        overlap_list.append(max_row)

        # Remove from target pool (cannot support more than one variant)
        df_target_chr.drop(max_row['TARGET_ID'], inplace=True)

        if df_target_chr.shape[0] == 0:
            break  # No more target matches to process

    # Merge and return
    if len(overlap_list) == 0:
        return pd.DataFrame([], columns=('ID', 'TARGET_ID', 'OFFSET', 'RO', 'SZRO', 'OFFSZ'))

    return pd.concat(
        overlap_list, axis=1
    ).T.loc[
       :, ['ID', 'TARGET_ID', 'OFFSET', 'RO', 'SZRO', 'OFFSZ']
    ]
