"""
Sequence manipulation routines.
"""

import gzip
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def bed_to_seqrecord_iter(df, skip_na=True):
    """
    Get a `SeqRecord` iterator for records in a dataframe. Record IDs use the ID column, and record sequences use
    the SEQ column.

    :param df: Dataframe with `ID` and `SEQ`.
    :param skip_na: Skip NA values without throwing an error.

    :return: An iterator for `SeqRecord` objects.
    """

    # Iterate over dataframe and return SeqRecords

    for index, row in df.iterrows():

        if pd.isnull(row['SEQ']) and skip_na:
            continue

        yield SeqRecord(
            Seq(row['SEQ']),
            id=row['ID'],
            description=''
        )


def fa_to_record_iter(fa_file_name):
    """
    Open a FASTA file, iterate through records. Returns an iterator of SeqIO.Seq objects.

    :param fa_file_name: FASTA file name. May be gzipped.

    :return: SeqIO.Seq iterator.
    """

    with PlainOrGzReader(fa_file_name) as in_file:
        for record in SeqIO.parse(in_file, 'fasta'):
            yield record

def fa_to_series(fa_file_name):
    """
    Read records from a FASTA file and generate a Pandas Series with IDs an keys and sequences as values. FASTA may
    be gzipped.

    :param fa_file_name: FASTA file name.

    :return: Pandas Series with IDs as keys and sequences as values.
    """

    df_series = pd.Series(
        {record.id: str(record.seq) for record in fa_to_record_iter(fa_file_name)}
    )

    df_series.name = 'SEQ'

    return df_series

class PlainOrGzReader:
    """
    Read a plain or a gzipped file using context guard.

    Example:
        with PlainOrGzReader('path/to/file.gz'): ...
    """

    def __init__(self, file_name, mode='rt'):

        self.file_name = file_name

        self.is_gz = file_name.lower().endswith('.gz')
        self.mode = mode

        self.file_handle = None

    def __enter__(self):

        if self.is_gz:
            self.file_handle = gzip.open(self.file_name, self.mode)
        else:
            self.file_handle = open(self.file_name, self.mode)

        return self.file_handle

    def __exit__(self, exc_type, exc_value, traceback):

        if self.file_handle is not None:
            self.file_handle.__exit__()
            self.file_handle = None
