"""
Tools for managing Snakemake resources.
"""

import re
import snakemake


def nlset(named_list, key, value=None, wildcards=None):
    """
    Set a value on a named list in Snakemake. This includes objects for wildcards, input, output, parameters, and log.

    :param named_list: Named list (wildcards, input, output, parameters, etc).
    :param key: Key.
    :param value: Value to set. May be a string or an input function taking a single argument, `wildcards`.
        This parameter must be set if `value` is a function. If `None`, split from key (first "=").
    :param wildcards: Format `value` with `wildcards` if set. If `value` is a function, call it with
        `wildcards` as its only parameter.
    """

    if value is None:
        key, value = key.split('=', 1)

        value = value.strip()

        # Remove quotes around the value (parsing key="value" or key='value')
        if len(value) > 2:
            if (value[0] == value[-1]) and (value[0] in {'"', '\''}):
                value = value[1:-1]

    if callable(value):
        # Input function
        if wildcards is None:
            raise RuntimeError('Cannot execute input function with wildcards = None')

        value = value(wildcards)

    else:
        # Format
        if wildcards is not None:
            if type(value) == str:
                value = re.sub('\{([^,\}]+),[^\}]+\}', '{\\1}', value)  # Remove regex qualifiers from wildcards
                value = value.format(**wildcards)  # Format wildcards into value
            elif type(value) == list:
                value = [
                    re.sub('\{([^,\}]+),[^\}]+\}', '{\\1}', item).format(**wildcards) for item in value
                ]

    # Set value
    snake_version_tok = [int(val) for val in snakemake.__version__.split('.')]

    if snake_version_tok[0] > 5 or (snake_version_tok[0] == 5 and snake_version_tok[1] >= 4):
        # Add key if missing
        if key not in named_list.keys():
            named_list.append(value)
            named_list._add_name(key)

        setattr(named_list, key, value)
        named_list[named_list._names[key][0]] = value

    else:
        # Add key if missing
        if key not in named_list.keys().keys():
            named_list.append(None)
            named_list.add_name(key)

        # Set value
        setattr(named_list, key, value)
        named_list[named_list._names[key][0]] = value
