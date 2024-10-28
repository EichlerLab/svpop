"""
Tools for managing Snakemake resources.
"""

import re


def strip_and_format(value, wildcards=None):
    if value is None:
        return value

    # Strip whitespace and commas
    value = value.strip()

    while value.endswith(','):
        value = value[:-1]

    # Remove regex qualifiers
    value = re.sub(r'\{([^,\}]+),[^\}]+\}', r'{\1}', value)

    # Strip quotes
    value = re.sub(r"""^\s*(['"]+)\s*(.+)\s*\1\s*$""", r'\2', value)

    # Remove temp()
    value = re.sub(r"""^(temp\s*\((['"]+)(.+)\s*\2\s*\))?\s*,*\s*$""", r'\3', value)

    # Format
    if wildcards is not None:
        value = value.format(**wildcards)

    return value


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

    if key is None:
        raise RuntimeError('Key cannot be None')

    if value is None and isinstance(key, str) and '=' in key:
        key, value = key.split('=', 1)

        value = strip_and_format(value, wildcards)

    elif value is not None and callable(value):
        # Input function
        if wildcards is None:
            raise RuntimeError('Cannot execute input function with wildcards = None')

        value = value(wildcards)

    else:
        # Format
        if isinstance(value, str):
            strip_and_format(value, wildcards)
        elif isinstance(value, (list, tuple, set)):
            value = [
                strip_and_format(item, wildcards) for item in value
            ]

    # Add key if missing
    if key not in named_list.keys():
        named_list.append(value)
        named_list._add_name(key)

    setattr(named_list, key, value)
    named_list[named_list._names[key][0]] = value
