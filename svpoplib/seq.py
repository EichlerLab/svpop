"""
Sequence manipulation routines.
"""

import gzip
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def bed_to_seqrecord_iter(df, record_set=None, require_all=True, skip_na=True, seq_col='SEQ'):
    """
    Get a `SeqRecord` iterator for records in a dataframe. Record IDs use the ID column, and record sequences use
    the SEQ column.

    :param df: Dataframe with `ID` and `SEQ`.
    :param record_set: Set of record names to extract. Ignore all records with an ID not in this set. May be a `dict` if
        this function should translate records from dict keys to dict values; all IDs that are not a dict key will be
        ignored.
    :param require_all: If `True`, throw an error if all records from record_set are not found.
    :param skip_na: Skip NA values without throwing an error.
    :param seq_col: Column in `df` that contains the sequence to be written.

    :return: An iterator for `SeqRecord` objects.
    """

    # Check columns
    missing_col = [val for val in ('ID', seq_col) if val not in df.columns]

    if missing_col:
        raise RuntimeError('Missing required column(s) for BED seqrecord conversion: {}'.format(','.join(missing_col)))

    # Setup record_set and record_dict for translating record names
    if record_set is not None:
        if issubclass(record_set.__class__, dict):
            record_dict = record_set
            record_set = set(record_dict.keys())

        elif issubclass(record_set.__class__, set):
            record_dict = {val: val for val in record_set}

        else:
            raise RuntimeError(f'Parameter record_set must be a set or a dict: {record_set.__class__}')

    else:
        record_set = None
        record_dict = None

    # Iterate over dataframe and return SeqRecords
    found_set = set()

    for index, row in df.iterrows():

        record_id = row['ID']
        record_seq = row[seq_col]

        record_id_org = record_id

        # Translate record names or skip if not an expected record
        if record_dict:
            record_id = record_dict.get(record_id, None)

            # Skip if
            if record_id is None:
                continue

        # Skip or error if no sequence exists for record
        if pd.isnull(record_seq):
            if skip_na:
                continue
            else:
                raise RuntimeError('Error converting record to sequence: Missing sequence: {}'.format(record_id_org))

        # Check for duplicates
        if record_id in found_set:
            raise RuntimeError('Duplicate BED entries for record: {}{}'.format(
                record_id, f' (original ID {record_id_org})' if record_id != record_id_org else ''
            ))

        # Yield record
        found_set.add(record_id)

        yield SeqRecord(
            Seq(record_seq),
            id=record_id,
            description=''
        )

    # Check for required records
    if require_all and record_set is not None and found_set != record_set:
        missing_set = record_set - found_set

        raise RuntimeError('Missing {} records when parsing {} {}: {}{}'.format(
            len(missing_set),
            input_format.upper(),
            fa_file_name,
            ', '.join(sorted(missing_set)[:3]),
            '...' if len(missing_set) > 3 else ''
        ))


def fa_to_record_iter(fa_file_name, record_set=None, require_all=True, input_format='fasta'):
    """
    Open a FASTA file, iterate through records. Returns an iterator of SeqIO.Seq objects.

    :param fa_file_name: FASTA file name. May be gzipped.
    :param record_set: Set of record names to extract. Ignore all records with an ID not in this set. May be a `dict` if
        this function should translate records from dict keys to dict values; all IDs that are not a dict key will be
        ignored.
    :param require_all: If `True`, throw an error if all records from record_set are not found.
    :param input_format: Input file format. Must be "fasta" or "fastq" (case sensitive).

    :return: SeqIO.Seq iterator.
    """

    # Check file format
    if input_format not in {'fasta', 'fastq'}:
        raise RuntimeError(f'Unrecognized input format: "{input_format}"')

    # Setup record_set and record_dict for translating record names
    if record_set is not None:
        if issubclass(record_set.__class__, dict):
            record_dict = record_set
            record_set = set(record_dict.keys())

        elif issubclass(record_set.__class__, set):
            record_dict = {val: val for val in record_set}

        else:
            raise RuntimeError(f'Parameter record_set must be a set or a dict: {record_set.__class__}')

    else:
        record_set = None
        record_dict = None

    # Open and parse records
    found_set = set()

    with PlainOrGzReader(fa_file_name) as in_file:
        for record in SeqIO.parse(in_file, input_format):

            record_id_org = record.id

            # Translate record names or skip if not an expected record
            if record_dict:
                record_id = record_dict.get(record.id, None)

                # Skip if
                if record_id is None:
                    continue

                record.id = record_id

            # Check for duplicate records
            if record.id in found_set:
                raise RuntimeError('Duplicate {} entries for record: {}{}'.format(
                    input_format.upper(),
                    record.id,
                    f' (original ID {record_id_org})' if record.id != record_id_org else ''
                ))

            # Yield record
            found_set.add(record.id)

            yield record

    # Check for required records
    if require_all and record_set is not None and found_set != record_set:
        missing_set = record_set - found_set

        raise RuntimeError('Missing {} records when parsing {} {}: {}{}'.format(
            len(missing_set),
            input_format.upper(),
            fa_file_name,
            ', '.join(sorted(missing_set)[:3]),
            '...' if len(missing_set) > 3 else ''
        ))


def gfa_to_record_iter(gfa_file_name, record_set=None, require_all=True):
    """
    Open a GFA file and parse "S" lines into sequence records. Returns an iterator of SeqIO.Seq objects.

    :param gfa_file_name: GFA file name. May be gzipped.
    :param record_set: Set of record names to extract. Ignore all records with an ID not in this set. May be a `dict` if
        this function should translate records from dict keys to dict values; all IDs that are not a dict key will be
        ignored.
    :param require_all: If `True`, throw an error if all records from record_set are not found.

    :return: SeqIO.Seq iterator.
    """

    found_set = set()

    # Setup record_set and record_dict for translating record names
    if record_set is not None:
        if issubclass(record_set.__class__, dict):
            record_dict = record_set
            record_set = set(record_dict.keys())

        elif issubclass(record_set.__class__, set):
            record_dict = {val: val for val in record_set}

        else:
            raise RuntimeError(f'Parameter record_set must be a set or a dict: {record_set.__class__}')

    else:
        record_set = None
        record_dict = None

    # Open and parse records
    with PlainOrGzReader(gfa_file_name) as in_file:

        line_count = 0

        for line in in_file:
            line_count += 1

            # Get record ID and sequence
            tok = line.split('\t')

            if tok[0] != 'S':
                continue

            if len(tok) < 3:
                raise RuntimeError(f'Error reading GFA "S" record at line {line_count}: Expected at least 3 tab-separated columns: {gfa_file_name}')

            record_id = tok[1].strip()
            record_seq = tok[2].strip()

            if not record_seq:
                continue

            # Check for duplicate IDs
            if record_id in found_set:
                raise RuntimeError('Duplicate GFA entries for record: {}'.format(record_id))

            found_set.add(record_id)

            # Translate record names
            if record_dict:
                record_id = record_dict.get(record_id, record_id)

            # Yield record
            yield SeqRecord(
                Seq(record_seq),
                id=record_id,
                description=''
            )

    # Check for required records
    if require_all and record_set is not None and found_set != record_set:
        missing_set = record_set - found_set

        raise RuntimeError('Missing {} records when parsing GFA {}: {}{}'.format(
            len(missing_set),
            gfa_file_name,
            ', '.join(sorted(missing_set)[:3]),
            '...' if len(missing_set) > 3 else ''
        ))


def fa_to_series(fa_file_name):
    """
    Read records from a FASTA file and generate a Pandas Series with IDs and keys and sequences as values. FASTA may
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
