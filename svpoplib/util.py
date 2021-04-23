"""
Miscelaneous utilities.
"""

import inspect
import os
import pandas as pd
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
            raise RuntimeError('None value cannot be interpreted as a boolean')

        if not issubclass(none_val.__class__, bool):
            raise RuntimeError('as_bool() got a None value and a non-bool None-replacement value: {}'.format(str(none_val.__class__)))

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
