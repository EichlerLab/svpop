"""
Miscelaneous utilities.
"""

import inspect
import os
import pandas as pd
import string
import sys

import traceback


def is_int(val):
    """
    Determine if value is castable as an integer.

    :param val: Value.

    :return: `True` if `val` is an integer and `False` otherwise.
    """
    try:
        int(val)
    except:
        return False

    return True


def get_install_dir():
    """
    Get the directory where the analysis pipeline is installed.

    :return: Analysis pipeline directory name (full path).
    """
    # Borrowed code from:
    # https://stackoverflow.com/questions/50499/how-do-i-get-the-path-and-name-of-the-file-that-is-currently-executing

    return os.path.dirname(
        os.path.dirname(
            os.path.abspath(
                inspect.getfile(
                    inspect.currentframe()
                )
            )
        )
    )


def get_traceback_details(e):
    """
    Get exception traceback details.

    Code from:
    https://stackoverflow.com/questions/1278705/when-i-catch-an-exception-how-do-i-get-the-type-file-and-line-number

    :param e: Exception.

    :return: Dictionary of details.
    """

    try:
        raise e

    except:

        exc_type, exc_value, exc_traceback = sys.exc_info()

        traceback_details = {
            'filename': exc_traceback.tb_frame.f_code.co_filename,
            'lineno': exc_traceback.tb_lineno,
            'name': exc_traceback.tb_frame.f_code.co_name,
            'type': exc_type.__name__,
            'message': str(exc_value),
            'trace':  traceback.format_exc()
        }

    return traceback_details


def as_bool(val, none_val=None):
    """
    Translate value as a boolean. If `val` is boolean, return `val`. If val is a string, `True` if lower-case string
    is "true", "1", "yes", "t", or "y". `False` if lower-case string is "false", "0", "no", "f", or "n". All other
    values throw a RuntimeError. If `val` is not bool or string, it is converted to a string and the rules above are
    applied.

    :param val: Value to interpret as a boolean.
    :param none_val: Return this value if val is `None`. If this value is also `None`, throw an error.

    :return: `True` or `False` (see above).
    """

    # Handle None
    if val is None:
        if none_val is None:
            raise ValueError('None value cannot be interpreted as a boolean')

        if not issubclass(none_val.__class__, bool):
            raise ValueError('as_bool() got a None value and a non-bool None-replacement value: {}'.format(str(none_val.__class__)))

        return none_val

    # Return if val is already boolean
    if issubclass(val.__class__, bool):
        return val

    # Interepret as bool
    val = str(val).lower()

    if val in {'true', '1', 'yes', 't', 'y'}:
        return True

    if val in {'false', '0', 'no', 'f', 'n'}:
        return False

    # Cannot interpret as bool
    raise ValueError('Cannot interpret as boolean value: {}'.format(val))


def parse_param_string(param_string):
    """
    Parse parameter string into a dictionary. Parameters are semi-colon delimited with attr=value pairs for each entry.
    If the equal sign is missing, `None` becomes the value of the pair. A dictionary of attribute-value pairs is
    returned (attrib = dict key).

    :param param_string: Parameter string. If `None`, NA (Pandas), or an empty string, return an empty dict.

    :return: Dictionary of parameters and values.
    """

    if param_string is None or pd.isnull(param_string) or (issubclass(param_string.__class__, str) and not param_string.strip()):
        return dict()

    return {attr: val for attr, val in
        [
            (tuple((val.strip() for val in element)) if len(element) == 2 else (element[0], None)) for element in
                [attrib_val.split('=', 1) for attrib_val in param_string.split(';')]
        ]
    }


def cmp_ver(ver_a, ver_b):
    """
    Compare version strings element-by-element. Each version string must be composed of integers separated by
    dots (e.g. "1.0", or "1.0.1").

    :param ver_a: Version string.
    :param ver_b: Version string.

    :return: -1 if `ver_a` is smaller, 1 if `ver_b` is smaller, 0 otherwise.
    """

    try:
        ver_a = tuple((int(val) for val in ver_a.split('.')))
    except ValueError as e:
        raise ValueError(f'Error converting version "{ver_a}" to a list of integers: {e}')

    try:
        ver_b = tuple((int(val) for val in ver_b.split('.')))
    except ValueError as e:
        raise ValueError(f'Error converting version "{ver_b}" to a list of integers: {e}')

    return (ver_a > ver_b) - (ver_a < ver_b)

def format_cards(template, **kwargs):
    """
    Format a string with wildcards using kwargs. Any wildcard values missing will be left in the string for the next
    round.

    Credit:
    https://github.com/snakemake/snakemake/issues/124
    https://stackoverflow.com/questions/11283961/partial-string-formatting

    :param template: String template.
    :param kwargs: Wildcard values.

    :return: Parsed string.
    """

    class FormatDict(dict):
        def __missing__(self, key):
            return "{" + key + "}"

    formatter = string.Formatter()

    mapping = FormatDict(**kwargs)

    return formatter.vformat(template, (), mapping)
