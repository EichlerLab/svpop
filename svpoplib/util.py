"""
Miscelaneous utilities.
"""

import inspect
import os
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

def as_bool(val):
    """
    Translate value as a boolean. If `val` is boolean, return `val`. If val is a string, `True` if lower-case string
    is "true", "1", "yes", "t", or "y". `False` if lower-case string is "false", "0", "no", "f", or "n". All other
    values throw a RuntimeError. If `val` is not bool or string, it is converted to a string and the rules above are
    applied.

    :param val: Value to interpret as a boolean.

    :return: `True` or `False` (see above).
    """

    if issubclass(val.__class__, bool):
        return val

    val = str(val).lower()

    if val in {'true', '1', 'yes', 't', 'y'}:
        return True

    if val in {'false', '0', 'no', 'f', 'n'}:
        return False

    raise ValueError('Cannot interpret as boolean value: {}'.format(val))
