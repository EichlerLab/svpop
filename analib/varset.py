"""
Rules for managing variant sets.

Variant sets are typically published calls. They may be a "set" (variants from one study) or a "merged" set, which is
several "sets" merged together into a non-redundant set.
"""


import copy

varset_file_pattern_single = 'temp/variant/varset/{varset}/bed/{sample}/all/all/byref/{vartype}_{svtype}.bed'
varset_file_pattern_merged = 'results/variant/varset/{varset}/bed/{sample}/all/all/byref/{vartype}_{svtype}.bed'


def get_config_entry(varset_name, wildcards, config):
    """
    Get variant set configuration entry.

    :param varset_name: Name of the variant set. If the name contains a dash, then everything before the dash is assumed
        to be the varset name with the remainder being some subset defined by the varset.
    :param wildcards: Rule wildcards. If 'None', then the "input_bed" section of the config is None.
    :param config: Snakemake config.

    :return: A dictionary object with configuration details.
    """

    varset_name_full = varset_name.strip()
    varset_name = varset_name.split('-', 1)[0]  # Remove subset

    # Check config
    if 'varset' not in config:
        raise RuntimeError('No "varset" section in config.')

    if 'set' not in config['varset']:
        raise RuntimeError('No "set" in config["varset"]')

    # Get configuration sections
    varset_config_set = config['varset']['set']

    if 'merge' in config['varset']:
        varset_config_merge = config['varset']['merge']
    else:
        varset_config_merge = {}

    # Check for varset name
    if varset_name in varset_config_set:

        varset_entry = copy.deepcopy(varset_config_set[varset_name])
        varset_entry['type'] = 'set'

        if wildcards is not None:
            fmt_dict = dict(wildcards).copy()

            if 'varset' not in fmt_dict:
                fmt_dict['varset'] = varset_name_full

            varset_entry['input_bed'] = [
                varset_file_pattern_single.format(**fmt_dict)
            ]
        else:
            varset_entry['input_bed'] = None

    elif varset_name in varset_config_merge:

        varset_entry = copy.deepcopy(varset_config_merge[varset_name])
        varset_entry['type'] = 'merge'

        if 'varset' not in varset_entry:
            raise RuntimeError(
                'Merged variant set config definition is missing variant sets to merge (element "varset"): {}'.format(
                    varset_name
                )
            )

        varset_list = varset_entry['varset']

        # Break if a varset is inside a varset definition
        for varset_element in varset_list:
            if varset_element in config['varset']['merge']:
                raise RuntimeError(
                    'Recursive variant set name: Found a merged varset name inside another merged varset name: {} in {}'.format(
                        varset_element, varset_name
                    )
                )

        # Save BED list
        if wildcards is not None:
            parse_dict = dict(wildcards)

            input_bed_list = list()

            for varset_element in varset_list:
                parse_dict['varset'] = varset_element
                input_bed_list.append(varset_file_pattern_merged.format(**parse_dict))

            varset_entry['input_bed'] = input_bed_list

        else:
            varset_entry['input_bed'] = None

        # Set mapping strategy
        if 'merge_strategy' not in varset_entry or not varset_entry['merge_strategy'].strip():
            varset_entry['merge_strategy'] = 'nr:szro=50:offset=200'

    else:
        raise RuntimeError('Unrecognized variant set name: {}'.format(varset_name))

    # Number of input files
    if varset_entry['input_bed'] is not None:
        varset_entry['n_input_bed'] = len(varset_entry['input_bed'])
    else:
        varset_entry['n_input_bed'] = None

    # Set name and description
    if 'name' not in varset_entry:
        varset_entry['name'] = varset_name

    if 'description' not in varset_entry:
        varset_entry['description'] = varset_entry['name']

    # Return
    return varset_entry
