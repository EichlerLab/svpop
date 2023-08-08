"""
Utilities for pandas dataframes.
"""

import pandas as pd
import numpy as np
import threading
import traceback
import sys

from multiprocessing import Pool

import svpoplib.util

def read_csv_chrom(
        in_file_name, chrom=None, chunksize=5000, chrom_field='#CHROM', **kwargs
):
    """
    Read a DataFrame in chunks and save only records for a chromosome. Prevents reading a whole DataFrame into memory
    when only a chromosome is needed.

    :param in_file_name: Input file name.
    :param chrom: Subset to this chromosome (or None to read all). May be a string (single chromosome name) or
        a list, tuple, or set of chromosome names.
    :param chunksize: Size of chunks to read.
    :param chrom_field: Chromosome field name. Elements of `chrom` are matched against this field.
    :param kwargs: Additional arguments to `pd.read_csv`.

    :return: A Pandas DataFrame.
    """

    # Check for base case (no chrom set, just return the dataframe)
    if chrom is None:
        return pd.read_csv(in_file_name, **kwargs)

    # Check chrom (make into a set)
    if isinstance(chrom, str):
        chrom = {chrom}
    elif isinstance(chrom, list) or isinstance(chrom, tuple):
        chrom = set(chrom)
    elif not isinstance(chrom, set):
        raise RuntimeError(
            f'Unrecognized type for parameter type "chrom": {type(chrom)}: Must be string, list, tuple, or set'
        )

    # Check for dtypes
    if kwargs is None:
        kwargs = dict()

    if 'dtype' not in kwargs:
        kwargs['dtype'] = {'#CHROM': 'str'}
    elif '#CHROM' not in kwargs['dtype']:
        kwargs['dtype']['#CHROM'] = 'str'

    # Read
    df_list = list()
    col_list = None  # Save list of columns and return appropriate columns if no rows are found

    df_iter = pd.read_csv(in_file_name, iterator=True, chunksize=chunksize, **kwargs)

    for df in df_iter:

        if col_list is None:
            col_list = df.columns

        if chrom_field not in df.columns:
            raise RuntimeError(f'Cannot subset "{in_file_name}" by chromosome: No "{chrom_field}" field')

        df_list.append(df.loc[df['#CHROM'].isin(chrom)])

    # Return
    if len(df_list) > 0:
        return pd.concat(df_list, axis=0)
    else:
        return pd.DataFrame([], columns=col_list)


def _apply_parallel_cb_result(index, df_split_results, p_pool, thread_done):
    """
    Get a function to save results.

    :param index:
    :param df_split_results:

    :return: Dataframe results if successful, exception object if not.
    """

    def callback_handler(subdf):

        try:
            df_split_results[index] = subdf

            thread_done[index] = True

            if isinstance(subdf, Exception):
                p_pool.e_list.append(subdf)
                p_pool.has_error = True

            with p_pool.lock:
                p_pool.batch_complete.notify_all()

        except Exception as ex:
            pass

        finally:
            thread_done[index] = True

    return callback_handler


# def _apply_parallel_cb_error(index, df_split, p_pool, thread_done):
#     """
#     Error handler for apply_parallel to accurately report the source of errors.
#
#     :param index: Table index.
#     :param df_split: Split dataframe.
#     :param p_pool: Process pool object.
#     :param thread_done: Thread status array.
#
#     :return: Callback handler function.
#     """
#
#     def callback_handler(e_cause):
#         print('Reporting error in thread {}: {}'.format(index, str(e_cause)))
#         sys.stdout.flush()
#
#         try:
#             with p_pool.lock:
#
#                 # Record error (in protected code to preserve exception order)
#                 traceback = svpoplib.util.get_traceback_details(e_cause)
#
#                 e = RuntimeError(
#                     'Runtime in split {} (index {} - {}): {} [{} line {}]'.format(
#                         index, df_split[index].index[0], df_split[index].index[-1],
#                         traceback['message'], traceback['filename'], traceback['lineno']
#                     )
#                 )
#
#                 e.__cause__ = e_cause
#
#                 p_pool.e_list.append(e)
#
#                 print(traceback['trace'], file=sys.stderr)
#                 sys.stderr.flush()
#
#                 p_pool.has_error = True
#
#                 thread_done[index] = True
#
#                 print('Notifying batch_complete from handler for index {}'.format(index))
#                 p_pool.batch_complete.notify_all()
#                 print('Done notifying batch_complete from handler for index {}'.format(index))
#
#         except Exception as ex:
#             print('Error in callback handler: {}'.format(str(ex)))
#             p_pool.e_list.append(ex)
#
#         finally:
#             thread_done[index] = True
#
#     return callback_handler


def _apply_func(df, func, axis=1, kwds={}):
    """
    Apply a function to a DataFrame. Called by multiprocessing.Pool.apply_async().

    :param df: DataFrame.
    :param func: Function.
    :param axis: Axis to apply function.
    :param kwds: Keyword arguments to `func`.

    :return: Return value from `df.apply()`.
    """

    try:
        if kwds is None:
            kwds = {}

        return df.apply(func, axis=axis, **kwds)

    except Exception as ex:
        traceback.print_exc()  # Print exception and trace
        return ex


class ParallelPool:
    """
    Associate a pool object with an active flag and a lock object.
    """

    def __init__(self, pool):
        self.pool = pool  # Thread pool

        self.lock = threading.RLock()  # Synchronization lock for pool
        self.batch_complete = threading.Condition(self.lock)

        self.has_error = False  # Set to True if an error occurs
        self.e_list = list()  # List of exceptions thrown during analysis


def apply_parallel(df, func, n_part, n_core=None, kwds=None, verbose=False):
    """
    Split a dataframe into `n_part` partitions and apply function `func` to each part concurrently in up to `n_core`
    parallel jobs. Merge the resulting dataframes back into one and return the results.

    From:
    http://www.racketracer.com/2016/07/06/pandas-in-parallel/
    * Viewed 2018-10-08

    :param df: Dataframe `func` will be called on to obtain the final results.
    :param func: A function that takes a dataframe (as its first positional argument) and returns a dataframe.
    :param partitions: Number of partitions `df` is split into before applying `func`.
    :param cores: Maximum number of cores. If `None`, then it is set to `n_part`.
    :params kwds: Keyword arguments to `func`.
    :params verbose: Print thread status information if `True`.

    :return: A transformed and merged dataframe. Dataframe may be out of order and should be re-sorted.
    """

    # Check arguments
    n_part = np.int32(n_part)

    if n_part > df.shape[0]:
        n_part = df.shape[0]

    if n_core is None:
        n_core = n_part
    else:
        n_core = np.int32(n_core)

    if not isinstance(df, pd.DataFrame):
        raise RuntimeError('df is not a DataFrame')

    # Do split
    df_split = np.array_split(df, n_part)

    # Create an array to save results for each split
    df_split_results = [None] * len(df_split)
    thread_done = [False] * len(df_split)

    # Open pool of workers
    p_pool = ParallelPool(Pool(n_core))

    # Submit each split part to the worker pool
    try:
        for index in range(len(df_split)):
            p_pool.pool.apply_async(
                _apply_func, (df_split[index], func, 1, kwds), {},
                _apply_parallel_cb_result(index, df_split_results, p_pool, thread_done),
                None #_apply_parallel_cb_error(index, df_split, p_pool, thread_done)
            )

    except Exception as ex:
        print('Error in submitting jobs: {}'.format(str(ex)))

        with p_pool.lock:
            p_pool.e_list.append(ex)
            p_pool.has_error = True

    while True:
        if verbose:
            print('Checking threads...')

        # Stop waiting if all threads are done
        if np.all(thread_done):
            if verbose:
                print('All threads done')

            break

        with p_pool.lock:
            p_pool.batch_complete.wait(30)

        # if p_pool.has_error:
            # print('Terminating pool...')
            # p_pool.pool.terminate()  # Causes a deadlock in pool
            # print('Done terminating pool.')

            break

    if verbose:
        print('Joining pool...')

    p_pool.pool.close()
    p_pool.pool.join()

    if verbose:
        print('Join done')

    # Check for errors
    if p_pool.e_list:
        raise RuntimeError(f'Parallel job failed with exception: {p_pool.e_list[0]}') from p_pool.e_list[0]

    # Return merged dataframe
    return pd.concat(df_split_results, axis=0)


def concat_frames(frame_list):
    """
    Concatenate data frames and preserve the order of columns. The order of the columns is defined by the first
    dataframe. And new columns in subsequent frames are added to the column order preserving their order if there is
    more than one. This prevents merged dataframes from being sorted randomly after concatenating with disparate
    columns.

    :param frame_list: List of dataframes.

    :return: A concatenated dataframe.
    """

    if len(frame_list) == 0:
        raise RuntimeError('Cannot merge an empty set of data frames')

    # Get the list of columns
    col_list = []

    for frame in frame_list:
        col_list = col_list + [col for col in frame.columns if col not in col_list]

    # Concat
    df = pd.concat(frame_list, sort=False).loc[:, col_list]

    # Sort if at least #CHROM or ID are in the dataframe
    sort_list = [col for col in ('#CHROM', 'POS', 'END', 'ID') if col in df.columns]

    if '#CHROM' in sort_list or 'ID' in sort_list:
        df.sort_values(sort_list, inplace=True)

    # Return
    return df
