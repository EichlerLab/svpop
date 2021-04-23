"""
Functions for managing and coordinating caller sets (sets of variants from different callers for the same sample).
"""

import pandas as pd

import snakemake.io

import svpoplib.sm


def get_config_entry(callerset, config):
    """
    Get config entry for this caller set and check it.

    :params callerset: Caller set name.
    :params config: Snakemake config.
    """

    # Get attributes
    if 'callerset' not in config:
        raise RuntimeError('Config has no callerset section')

    callerset_entry = config['callerset'].get(callerset, None)

    if callerset_entry is None:
        raise RuntimeError('Missing definition in config for callerset: {}'.format(callerset))

    callerset_entry = callerset_entry.copy()

    # Check for path_list and name_list
    missing_list = [key for key in {'callsets', 'name_list', 'merge'} if key not in callerset_entry]

    if missing_list:
        raise RuntimeError(
            'Missing attributes in caller set definition: {}: missing = {}'.format(callerset,  ', '.join(missing_list))
        )

    # Check length
    if len(callerset_entry['callsets']) != len(callerset_entry['name_list']):
        raise RuntimeError(
            'Caller set definition has different lengths for "callsets" ({}) and "name_list" ({}): {}'.format(
                callerset,
                len(callerset_entry['callsets']),
                len(callerset_entry['name_list'])
            )
        )

    callerset_entry['n'] = len(callerset_entry['callsets'])

    # Set name and description
    if 'name' not in callerset_entry:
        callerset_entry['name'] = 'Callerset {}'.format(callerset)

    callerset_entry['name'] = callerset_entry['name'].strip()

    callerset_entry['callerset_name'] = callerset

    if 'description' not in callerset_entry or callerset_entry['description'] is None:
        callerset_entry['description'] = None
    else:
        callerset_entry['description'] = callerset_entry['description'].strip()

    # Get params
    callerset_entry['param_string'] = callerset_entry.get('params', None)
    callerset_entry['params'] = svpoplib.util.parse_param_string(callerset_entry['param_string'])

    # Return entry
    return callerset_entry


def get_caller_set_input(callerset, file_pattern, config, wildcards=None, as_tuple=False):
    """
    Get input for a caller set.

    :param file_pattern: File pattern with placeholders for any elements in `wildcards`. For each caller set entry, all
        instances of "{sourcetype}" and "{sourcename}" are replaced in this string.
    :param config: Snakemake config for this project.
    :param wildcards: Wildcards or `None` if file pattern contains no wildcards other than `samples` and `mergeset`
        is set.
    :param mergeset: Name of the mergeset. If `None`, it is set to `wildcards.mergeset`.
    :param as_tuple: Return as a list of tuples linking each sample name to a file "(sample, file)".

    :return: A list of input files with one file per samples.
    """

    # Clone or init wildcards
    wildcards = snakemake.io.Namedlist(toclone=wildcards)

    # Get configuration entry
    callerset_entry = get_config_entry(callerset, config)

    # Make list of input files
    file_list = list()

    for sourcetype, sourcename in callerset_entry['callsets']:
        svpoplib.sm.nlset(wildcards, 'sourcetype', sourcetype)
        svpoplib.sm.nlset(wildcards, 'sourcename', sourcename)

        file_list.append(
            (
                (sourcetype, sourcename),
                file_pattern.format(**wildcards)
            )
        )

    # Return list
    if as_tuple:
        return file_list

    return [file for source_tuple, file in file_list]


def merge_annotations(df_merge, callerset_input, callerset_entry, sort_columns=['#CHROM', 'POS', 'END', 'ID']):
    """
    Merge a table of annotations from several callers into one table.

    :param df_merge: Table of merged structural variants (BED file). Must contain columns 'ID',
        'CALLERSET_SOURCE', and 'CALLERSET_ORG_ID'.
    :param callerset_input: List of input tab or BED files.
    :param callerset_entry: Caller set configuration entry.
    :param sort_columns: List of column names to sort the merged annotations by or `None` to leave it unsorted.

    :return: Dataframe of merged annotations.
    """

    # Check arguments
    if len(callerset_input) != callerset_entry['n']:
        raise RuntimeError(
            'Input entry length ({}) does not match caller set definition length ({})'.format(
                len(callerset_input), callerset_entry['n']
            )
        )

    df_list = list()

    # Subset from each samples
    for index in range(callerset_entry['n']):
        anno_file = callerset_input[index]
        merged_name = callerset_entry['name_list'][index]

        # Get table of annotations
        df_anno = pd.read_csv(anno_file, sep='\t', header=0, low_memory=False)
        df_anno.set_index('ID', inplace=True, drop=False)

        # Subset
        df_merge_subset = df_merge.loc[df_merge['CALLERSET_SRC'] == merged_name]
        id_dict = {row['CALLERSET_SRC_ID']: row['ID'] for index, row in df_merge_subset.iterrows()}

        id_subset = set(df_merge_subset['CALLERSET_SRC_ID'])

        df_anno = df_anno.loc[df_anno['ID'].apply(lambda id: id in id_subset)]

        df_anno['ID'] = df_anno['ID'].apply(lambda svid: id_dict[svid])

        df_list.append(df_anno)

    # Merge subsets
    df_anno = pd.concat(df_list, axis=0)
    df_anno.reset_index(inplace=True, drop=True)

    # Sort by columns
    sort_columns = [col for col in sort_columns if col in df_anno.columns] if sort_columns is not None else []

    if sort_columns:
        df_anno.sort_values(list(sort_columns), inplace=True)

    # Return
    return df_anno


def fa_iter(df, callerset_entry, callerset_input):
    """
    Return FASTA records for callerset merges. Selects the lead variant FASTA record for each merged variant.

    :param df: Dataframe with ID, CALLERSET_SRC_ID and CALLERSET_SRC columns.
    :param callerset_entry: Callerset entry (`svpoplib.callerset.get_config_entry(wildcards.sourcename, config)`).
    :param callerset_input: List of input FASTA file in the same order as `callerset_entry['name_list']`

    :return: An iterator for selected FASTA records.
    """

    # Get parameters
    if 'fa_incomplete' in callerset_entry['params']:
        fa_incomplete = svpoplib.util.as_bool(callerset_entry['params'].get('fa_incomplete'), none_val=True)
    else:
        fa_incomplete = False

    for index in range(len(callerset_entry['name_list'])):
        caller = callerset_entry['name_list'][index]

        df_caller = df.loc[df['CALLERSET_SRC'] == caller]

        id_dict = dict(
            zip(
                df_caller['CALLERSET_SRC_ID'], df_caller['ID']
            )
        )

        for record in svpoplib.seq.fa_to_record_iter(
            callerset_input[index], id_dict, require_all=not fa_incomplete
        ):
            yield record


def cluster_param_cpu(wildcards, config):
    """
    Get number of cores to be allocated for variant merge jobs.
    """

    return int(
        svpoplib.sampleset.get_merge_strategy(
            get_config_entry(wildcards.sourcename, config),
            wildcards.vartype,
            wildcards.svtype
        ).get('cpu', svpoplib.sampleset.DEFAULT_RESOURCES['callerset']['cpu'])
    )


def cluster_param_mem(wildcards, config):
    """
    Get amount of memory to be allocated for variant merge jobs.
    """

    return \
        svpoplib.sampleset.get_merge_strategy(
            get_config_entry(wildcards.sourcename, config),
            wildcards.vartype,
            wildcards.svtype
        ).get('mem', svpoplib.sampleset.DEFAULT_RESOURCES['callerset']['mem'])


def cluster_param_rt(wildcards, config):
    """
    Get cluster runtime to be allocated for variant merge jobs.
    """

    return \
        svpoplib.sampleset.get_merge_strategy(
            get_config_entry(wildcards.sourcename, config),
            wildcards.vartype,
            wildcards.svtype
        ).get('rt', svpoplib.sampleset.DEFAULT_RESOURCES['callerset']['rt'])


def cluster_param_anno_mem(wildcards, config):
    """
    Get amount of memory to be allocated for callerset/sampleset annotation merge jobs.
    """

    return \
        svpoplib.sampleset.get_merge_strategy(
            get_config_entry(wildcards.sourcename, config),
            wildcards.vartype,
            wildcards.svtype
        ).get('anno_mem', svpoplib.sampleset.DEFAULT_RESOURCES['callerset']['anno_mem'])


# def cluster_param_cpu(wildcards, config):
#     """
#     Get number of cores to be allocated for variant merge jobs.
#     """
#
#     return int(get_config_entry(
#         wildcards.sourcename, config
#     ).get('cpu', '4'))
#
#
# def cluster_param_mem(wildcards, config):
#     """
#     Get amount of memory to be allocated for variant merge jobs.
#     """
#
#     return get_config_entry(
#         wildcards.sourcename, config
#     ).get('mem', '2G')
#
#
# def cluster_param_anno_mem(wildcards, config):
#     """
#     Get amount of memory to be allocated for variant merge jobs.
#     """
#
#     return get_config_entry(
#         wildcards.sourcename, config
#     ).get('anno_mem', '4G')
#
#
# def cluster_param_rt(wildcards, config):
#     """
#     Get cluster runtime to be allocated for variant merge jobs.
#     """
#
#     return get_config_entry(
#         wildcards.sourcename, config
#     ).get('rt', '48:00:00')
