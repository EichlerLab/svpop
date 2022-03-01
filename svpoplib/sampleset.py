"""
Functions for managing and coordinating sample sets (sets of variants from different samples).
"""

import gzip
import numpy as np
import pandas as pd
import snakemake.io
import svpoplib.sm

from Bio import SeqIO


# Default values (should be configurable)

# Default resource values
DEFAULT_RESOURCES = {
    'sampleset': {
        'cpu': '4',
        'mem': '12000',
        'rt': '48:00:00',
        'anno_mem': '12000'
    },
    'callerset': {
        'cpu': '4',
        'mem': '12000',
        'rt': '48:00:00',
        'anno_mem': '8000'
    }
}


def get_config_entry(sampleset_name, sampleset_list, config):
    """
    Get config entry for this sampleset and check it.

    :param sampleset_name: Name of the sample set.
    :param sampleset_list: Name of the list of samples. May contain a dash, and everything after the dash is parsed into
        "{id}" in the sample list and sample set configuration dictionary.
    :param config: Pipeline configuration.
    """

    # Split samplelist_name into list name and id (if it contains a dash)
    sampleset_list_tuple = sampleset_list.split('-')

    if len(sampleset_list_tuple) == 1:
        sampleset_list_name = sampleset_list_tuple[0]
        sampleset_id = None

    else:
        sampleset_list_name = sampleset_list_tuple[0]
        sampleset_id = sampleset_list_tuple[1].strip()

    if not sampleset_list_name:
        raise RuntimeError(
            'Sampleset list name is empty (format "name-id" where "-id" is optional): {}'.format(sampleset_list)
        )

    if not sampleset_id:
        sampleset_id = None
    else:
        sampleset_id_dict = {'id': sampleset_id}

    # Get attributes
    if 'sampleset' not in config:
        raise RuntimeError('Config has no sampleset section')

    sampleset_config = config['sampleset'].get(sampleset_name, None)

    if sampleset_config is None:
        raise RuntimeError('Missing definition in config for sampleset: {}'.format(sampleset_name))

    sampleset_config = sampleset_config.copy()

    # Parse entries
    sampleset_entry = dict()

    for key in sampleset_config.keys():
        val = sampleset_config[key]

        if sampleset_id is not None:
            sampleset_entry[key] = val.format(**sampleset_id_dict)
        else:
            sampleset_entry[key] = val

    # Set list of samples
    if 'samplelist' not in config:
        raise RuntimeError('Missing "samplelist" configuration section')

    if sampleset_list_name not in config['samplelist']:
        raise RuntimeError(
            'Missing samplelist definition: {}: Reference by sampleset {}'.format(sampleset_list_name, sampleset_name)
        )

    if sampleset_id is not None:
        sampleset_entry['samples'] = [
            element.format(**sampleset_id_dict) for element in config['samplelist'][sampleset_list_name]
        ]
    else:
        sampleset_entry['samples'] = config['samplelist'][sampleset_list_name]

    sampleset_entry['n'] = len(sampleset_entry['samples'])

    # Check entries
    missing_list = [
        attr for attr in {'sourcetype', 'sourcename', 'merge', 'name', 'description'} if attr not in sampleset_entry.keys()
    ]

    if missing_list:
        raise RuntimeError(
            'Sampleset definition "{}" is missing element(s): {}'.format(sampleset_name, ', '.join(missing_list))
        )

    if "sampleset_name" in sampleset_entry.keys():
        raise RuntimeError('Restricted key "sampleset_name" was found in sampleset {}'.format(sampleset_name))

    sampleset_entry['sampleset_name'] = sampleset_name

    # Get params
    sampleset_entry['param_string'] = sampleset_entry.get('params', None)
    sampleset_entry['params'] = svpoplib.util.parse_param_string(sampleset_entry['param_string'])

    # Return entry
    return sampleset_entry


def get_sample_set_input(sampleset_name, sampleset_list, file_pattern, config, wildcards=None, as_tuple=False):
    """
    Get input for a caller set.

    :param sampleset_name: Name of the sample set.
    :param sampleset_list: Name of the list of samples. May contain a dash, and everything after the dash is parsed into
        "{id}" in the sample list and sample set configuration dictionary.
    :param file_pattern: File pattern with placeholders for any elements in `wildcards`. For each caller set entry, all
        instances of "{sourcetype}", "{sourcename}", and "sample" and any additional wildcards are parsed into the
        string.
    :param config: Snakemake configuration.
    :param wildcards: Wildcards or `None` if file pattern contains no wildcards other than `sourcetype`, `sourcename`,
        and `sample`.
    :param mergeset: Name of the mergeset. If `None`, it is set to `wildcards.mergeset`.
    :param as_tuple: Return as a list of tuples linking each sample name to a file "(sample, file)".

    :return: A list of input files with one file per samples.
    """

    # Clone or init wildcards
    wildcards = snakemake.io.Namedlist(toclone=wildcards)

    # Get configuration entry
    sampleset_entry = get_config_entry(sampleset_name, sampleset_list, config)

    # Make list of input files
    file_list = list()

    svpoplib.sm.nlset(wildcards, 'sourcetype', sampleset_entry['sourcetype'])
    svpoplib.sm.nlset(wildcards, 'sourcename', sampleset_entry['sourcename'])

    for sample_name in sampleset_entry['samples']:
        svpoplib.sm.nlset(wildcards, 'sample', sample_name)

        file_list.append(
            (
                sample_name,
                file_pattern.format(**wildcards)
            )
        )

    # Return list
    if as_tuple:
        return file_list

    return [file for source_tuple, file in file_list]


def merge_annotations(df_merge, sampleset_input, sampleset_entry, sort_columns=['#CHROM', 'POS', 'END', 'ID']):
    """
    Merge a table of annotations from several samples into one table.

    :param df_merge: Table of merged structural variants (BED file). Must contain columns 'ID',
        'MERGE_SOURCE', and 'MERGE_ORG_ID'.
    :param sampleset_input: List of input tab or BED files.
    :param sampleset_entry: Caller set configuration entry.
    :param sort_columns: List of column names to sort the merged annotations by or `None` to leave it unsorted.

    :return: Dataframe of merged annotations.
    """

    # Check arguments
    if len(sampleset_input) != sampleset_entry['n']:
        raise RuntimeError(
            'Input entry length ({}) does not match sample set definition length ({})'.format(
                len(sampleset_input), sampleset_entry['n']
            )
        )

    df_list = list()

    # Subset from each samples
    for index in range(sampleset_entry['n']):
        anno_file = sampleset_input[index]
        sample = sampleset_entry['samples'][index]

        # Get table of annotations
        df_anno = pd.read_csv(anno_file, header=0, sep='\t', low_memory=False)
        df_anno.set_index('ID', inplace=True, drop=False)

        # Subset
        df_merge_subset = df_merge.loc[df_merge['MERGE_SRC'] == sample]
        id_dict = {row['MERGE_SRC_ID']: row['ID'] for index, row in df_merge_subset.iterrows()}

        id_subset = set(df_merge_subset['MERGE_SRC_ID'])

        df_anno = df_anno.loc[df_anno['ID'].apply(lambda id: id in id_subset)]

        df_anno['ID'] = df_anno['ID'].apply(lambda svid: id_dict[svid])

        df_list.append(df_anno)

    # Merge subsets
    df_anno = pd.concat(df_list, axis=0)
    df_anno.reset_index(inplace=True, drop=True)

    # Sort by columns
    if sort_columns is not None:
        sort_columns = [col for col in sort_columns if col in df_anno.columns]

        if sort_columns:
            df_anno.sort_values(list(sort_columns), inplace=True)

    # Return
    return df_anno


def cluster_param_cpu(wildcards, config):
    """
    Get number of cores to be allocated for variant merge jobs.
    """

    return int(
        get_merge_strategy(
            get_config_entry(wildcards.sourcename, wildcards.sample, config),
            wildcards.vartype,
            wildcards.svtype
        ).get('cpu', DEFAULT_RESOURCES['sampleset']['cpu'])
    )


def cluster_param_mem(wildcards, config):
    """
    Get amount of memory to be allocated for variant merge jobs.
    """

    return \
        get_merge_strategy(
            get_config_entry(wildcards.sourcename, wildcards.sample, config),
            wildcards.vartype,
            wildcards.svtype
        ).get('mem', DEFAULT_RESOURCES['sampleset']['mem'])


def cluster_param_rt(wildcards, config):
    """
    Get cluster runtime to be allocated for variant merge jobs.
    """

    return \
        get_merge_strategy(
            get_config_entry(wildcards.sourcename, wildcards.sample, config),
            wildcards.vartype,
            wildcards.svtype
        ).get('rt', DEFAULT_RESOURCES['sampleset']['rt'])


def cluster_param_anno_mem(wildcards, config):
    """
    Get amount of memory to be allocated for callerset/sampleset annotation merge jobs.
    """

    return \
        get_merge_strategy(
            get_config_entry(wildcards.sourcename, wildcards.sample, config),
            wildcards.vartype,
            wildcards.svtype
        ).get('anno_mem', DEFAULT_RESOURCES['sampleset']['anno_mem'])


def get_merge_strategy(sampleset_entry, vartype, svtype):
    """
    Get the merge strategy string for this set entry.

    Function works for both samplesets and callersets.

    :param sampleset_entry: Sampleset or callerset entry.
    :param vartype: Variant type.
    :param svtype: SV type.

    :return: Merge strategy string.
    """

    # Key formats
    # "merge": "strategy-def"
    # "merge": {"vartype1,vartype2:svtype1,svtype2": "strategy-def"}
    # "merge": {"vartype1,vartype2:svtype1,svtype2": {"strategy": "strategy-def", "cpu": "...", "mem": "..."}}
    # "merge": {"vartype1,vartype2": {"strategy": "strategy-def", "cpu": "...", "mem": "..."}}

    if 'sampleset_name' in sampleset_entry:
        set_name = sampleset_entry['sampleset_name']
        set_type = 'sampleset'
    elif 'callerset_name' in sampleset_entry:
        set_name = sampleset_entry['callerset_name']
        set_type = 'callerset'
    else:
        raise RuntimeError('Set entry does not appear to be a sampleset or callerset')

    # Get strategy entry
    merge_entry_dict = sampleset_entry['merge']

    if issubclass(merge_entry_dict.__class__, dict):

        best_match = None
        best_match_level = 0

        # Find best match
        for key in merge_entry_dict.keys():

            if key == 'DEFAULT' and best_match is None:
                best_match = key
                best_match_level = 1
                continue

            key_tok = key.split(':')

            if len(key_tok) > 2:
                raise RuntimeError(
                    'Key token in {} definition "{}" has more than one colon-separated field (expected vartype:svtype or vartype): {}'.format(
                        set_type, set_name, key_tok
                    )
                )

            if vartype in set(key_tok[0].split(',')):
                if len(key_tok) == 2:
                    if svtype in set(key_tok[1].split(',')) and best_match_level < 2:
                        best_match = key
                        best_match_level = 2

                elif best_match_level < 1:
                    best_match = key
                    best_match_level = 1

        # Get best match
        if best_match is None:
            raise RuntimeError(
                'Found no matching {} merge entries in {} (vartype={}, svtype={})'.format(
                    set_type, set_name, vartype, svtype
                )
            )

        merge_entry = merge_entry_dict[best_match]

        if issubclass(merge_entry.__class__, dict):

            if "strategy" not in merge_entry.keys():
                raise RuntimeError('Missing "strategy" in {} "{}" with key "{}"'.format(
                    set_type, set_name, best_match
                ))

            merge_strategy = merge_entry['strategy']

            merge_cpu = merge_entry.get('cpu', DEFAULT_RESOURCES[set_type]['cpu'])
            merge_mem = merge_entry.get('mem', DEFAULT_RESOURCES[set_type]['mem'])
            merge_rt = merge_entry.get('rt', DEFAULT_RESOURCES[set_type]['rt'])
            anno_mem = merge_entry.get('anno_mem', DEFAULT_RESOURCES[set_type]['anno_mem'])

        else:
            merge_strategy = merge_entry
            merge_cpu = DEFAULT_RESOURCES[set_type]['cpu']
            merge_mem = DEFAULT_RESOURCES[set_type]['mem']
            merge_rt = DEFAULT_RESOURCES[set_type]['rt']
            anno_mem = DEFAULT_RESOURCES[set_type]['anno_mem']

    else:
        # Same merge strategy for all vartype/svtype's and default resources
        merge_strategy = merge_entry_dict
        merge_cpu = DEFAULT_RESOURCES[set_type]['cpu']
        merge_mem = DEFAULT_RESOURCES[set_type]['mem']
        merge_rt = DEFAULT_RESOURCES[set_type]['rt']
        anno_mem = DEFAULT_RESOURCES[set_type]['anno_mem']

    merge_dict = {
        'strategy': merge_strategy,
        'cpu': merge_cpu,
        'mem': merge_mem,
        'rt': merge_rt,
        'anno_mem': anno_mem,
        'set_type': set_type,
        'set_name': set_name
    }

    return merge_dict


def fa_write_func(df, wildcards, sampleset_entry, fa_input_pattern, config):
    """
    Function to yield a sequence record iterator for rule variant_sampleset_fa_merge.

    :param df: DataFrame with SAMPLE, ID_SAMPLE (original pre-merged ID), and ID.
    :param wildcards: Rule wildcards.
    :param sampleset_entry: Sampleset entry.
    :param fa_input_pattern: FASTA file input pattern for input into the merge.
    :param config: SV-Pop config.

    :return: SeqRecord iterator.
    """

    # Get format processing items
    fa_dict = {sample: fa_file for sample, fa_file in svpoplib.sampleset.get_sample_set_input(
        wildcards.sourcename,
        wildcards.sample,
        fa_input_pattern,
        config,
        wildcards,
        as_tuple=True
    )}

    # Get parameters
    if 'fa_incomplete' in sampleset_entry['params']:
        fa_incomplete = svpoplib.util.as_bool(sampleset_entry['params'].get('fa_incomplete'), none_val=True)
    else:
        fa_incomplete = False

    # Process samples
    for sample in sampleset_entry['samples']:

        if np.sum(df['SAMPLE'] == sample) == 0:
            continue

        id_dict = dict(df.loc[df['SAMPLE'] == sample, ['ID_SAMPLE', 'ID']].set_index('ID_SAMPLE')['ID'])
        id_set = set(id_dict.keys())

        # Open file
        in_file_name = fa_dict[sample]

        found_ids = set()

        with gzip.open(in_file_name, 'rt') as in_file:
            for record in SeqIO.parse(in_file, 'fasta'):

                found_ids.add(record.id)

                if record.id in id_set:
                    record.id = id_dict[record.id]
                    yield record

        # Check for missing IDs
        if not fa_incomplete:
            missing_ids = id_set - found_ids

            if missing_ids:
                raise RuntimeError('Missing {} variant ID(s) for sample {} when merging FASTAs: {}{}: {}'.format(
                    len(missing_ids), sample,
                    ', '.join(sorted(missing_ids)[:3]), '...' if len(missing_ids) > 3 else '',
                    in_file_name
                ))


def is_read_seq(wildcards, config):
    """
    Determine if merge requires input sequence.

    :param wildcards: Rule wildcards.
    :param config: Configuration.

    :return: `True` if sequences should be read.
    """

    sampleset_entry = svpoplib.sampleset.get_config_entry(wildcards.sourcename, wildcards.sample, config)
    merge_strategy_tok = svpoplib.sampleset.get_merge_strategy(sampleset_entry, wildcards.vartype, wildcards.svtype)['strategy'].split(':', 1)

    if len(merge_strategy_tok) == 1:
        return False

    return svpoplib.svmerge.get_param_set(merge_strategy_tok[1], merge_strategy_tok[0]).read_seq
