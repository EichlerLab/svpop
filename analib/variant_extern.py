"""
Functions for variant calls pulled in from external sources.
"""


import os

import analib.util


def get_config(extern_source, config):
    """
    Get configuration dictionary for an external source and check the entry.

    :params wildcards: Rule wildcards.
    :params config: Pipeline config dict.
    """

    # Check for extern config section
    if 'variant_extern' not in config:
        raise RuntimeError(
            'Missing configuration section "variant_extern" while getting parameters for source: {}'.format(
                extern_source
            )
        )

    if extern_source not in config['variant_extern']:
        raise RuntimeError('Missing configuration in "variant_extern" for source: {}'.format(extern_source))

    # Get config section
    extern_entry = config['variant_extern'][extern_source]

    # Get svtype-combined
    if 'svtype_combined' not in extern_entry:
        extern_entry['svtype_combined'] = False
    else:
        try:
            extern_entry['svtype_combined'] = analib.util.as_bool(extern_entry['svtype_combined'])
        except ValueError as ex:
            raise RuntimeError('Field "svtype_combined" in extern config "{}" is not boolean: {}'.format(extern_source, str(ex)))

    # Check
    if not issubclass(extern_entry.__class__, dict):
        raise RuntimeError('"variant_extern" config entry for source must be a dictionary with a "bed" entry: {}'.format(extern_source))

    if 'bed' not in extern_entry:
        raise RuntimeError('Missing "bed" in "variant_extern" for source: {}'.format(extern_source))

    if '{vartype}' not in extern_entry['bed']:
        raise RuntimeError('Missing "{{vartype}}" wildcard in "bed" entry for external source: {}'.format(extern_source))

    if '{svtype}' not in extern_entry['bed']:
        if not extern_entry['svtype_combined']:
            raise RuntimeError('Missing "{{svtype}}" wildcard in "bed" entry for external source: {}'.format(extern_source))
    else:
        if extern_entry['svtype_combined']:
            raise RuntimeError('Found "{{svtype}}" wildcard in "bed", but "svtype_combined" is also set for external source: {}'.format(extern_source))

    if '{sample}' not in extern_entry['bed']:
        raise RuntimeError('Missing "{{sample}}" wildcard in "bed" entry for external source: {}'.format(extern_source))

    # Set name and description
    if 'name' not in extern_entry:
        extern_entry['name'] = 'Extern-{}'.format(extern_source)

        extern_entry['name'] = extern_entry['name'].strip()

    if 'description' not in extern_entry or extern_entry['description'] is None:
        extern_entry['description'] = None
    else:
        extern_entry['description'] = extern_entry['description'].strip()

    # Return
    return extern_entry


def get_bed(wildcards, config, missing_value=None):
    """
    Get BED file name for an external source.

    :params wildcards: Rule wildcards.
    :params config: Pipeline config dict.
    """

    # Get entry
    extern_entry = get_config(wildcards.extern_source, config)

    # Get BED
    bed_file_name = extern_entry['bed'].format(**wildcards)

    # If BED is missing and an empty flag is set, then allow writing an empty file
    if not os.path.isfile(bed_file_name):

        allow_empty = set([tuple(element) for element in extern_entry.get('empty', [])])

        if (wildcards.vartype, wildcards.svtype) in allow_empty:
            return missing_value

    # Return BED file name
    return bed_file_name


def get_fa(wildcards, config):
    """
    Get FASTA file name for an external source if it is defined ("fasta" config item) or `None` if it is not defined.

    :params wildcards: Rule wildcards.
    :params config: Pipeline config dict.
    """

    # Get entry
    extern_entry = get_config(wildcards.extern_source, config)

    # Get FASTA
    if 'fasta' not in extern_entry:
        return None

    fa_file_name = extern_entry['fasta'].format(**wildcards)

    if not os.path.isfile(fa_file_name):
        raise RuntimeError('Missing extern FASTA file: {}'.format(fa_file_name))

    return fa_file_name
