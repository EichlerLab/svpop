"""
Fast Smith-waterman implementation.
"""

import collections
import numpy as np

import kanapy

# Op codes
OP_NONE = 0
OP_MATCH = 1
OP_MISMATCH = 2
OP_GAP_SUB = 3  # Gap subject, insertion
OP_GAP_QRY = 4  # Gap query, deletion

OP_CODE = ['*', '=', 'X', 'I', 'D']  # Indexed by OP_ constants

# Score sources
SOURCE_NONE = 0
SOURCE_ALN = 1
SOURCE_GAP_SUB = 2
SOURCE_GAP_QRY = 3


class ScoreTraceNode:
    def __init__(self, op_code=OP_NONE, score=0):
        self.op_code = op_code
        self.score = score

    def __str__(self):
        return f'TraceNode({OP_CODE[self.op_code]}, score={self.score}, id={id(self):x}'


def get_kmer_count(seq, k_size):

    counter = collections.Counter()

    seq = seq.strip().upper()

    for index in range(len(seq) - k_size + 1):
        counter[seq[index:index + k_size]] += 1

    return counter


def jaccard_distance(seq_a, seq_b, k_size):
    """
    Get the Jaccard distance between k-merized sequences. This Jaccard distance is computed on the total number of
    k-mers including multiplicity (the same kmer may appear more than once). For example, if a k-mer is in A twice
    and B once, one instance of the k-mer matches and one does not match.

    :param seq_a: Sequence A (string).
    :param seq_b: Sequence B (string).
    :param k_size: K-mer size.

    :return: Jaccard distance account for multiplicity.
    """

    count1 = get_kmer_count(seq_a, k_size)
    count2 = get_kmer_count(seq_b, k_size)

    key_set = set(count1.keys()) | set(count2.keys())

    return np.sum(
        [np.min([count1[key], count2[key]]) for key in key_set]  # Matching k-mers
    ) / np.sum(
        [np.max([count1[key], count2[key]]) for key in key_set]  # All k-mers
    ) if len(count1) > 0 and len(count2) > 0 else 0


class ScoreAligner:
    """
    An engine for aligning two sequences to obtain an alignment score only, no alignment is generated. This is intended
    for determining sequence similarity through an alignment where the alignment itself is not important.

    :param match: Base match score (> 0.0).
    :param mismatch: Base mismatch score (< 0.0).
    :param gap_open: Gap open cost (<= 0.0).
    :param gap_extend: Gap extend cost (<= 0.0).
    :param map_limit: Maximum sequence size before falling back to Jaccard index in match_prop().
    :param jaccard_kmer: Jaccard k-mer size for comparisons falling back to Jaccard index in match_prop().
    """

    def __init__(self, match=2.0, mismatch=-1.0, gap_open=-1.0, gap_extend=-0.25, map_limit=20000, jaccard_kmer=9):

        # Check and assign values
        try:
            self.__match = float(match)
        except ValueError:
            raise RuntimeError(f'ScoreAligner(): Parameter match must be numeric: {match}')

        try:
            self.__mismatch = float(mismatch)
        except ValueError:
            raise RuntimeError(f'ScoreAligner(): Parameter mismatch must be numeric: {mismatch}')

        try:
            self.__gap_o = float(gap_open)
        except ValueError:
            raise RuntimeError(f'ScoreAligner(): Parameter gap_open must be numeric: {gap_open}')

        try:
            self.__gap_e = float(gap_extend)
        except ValueError:
            raise RuntimeError(f'ScoreAligner(): Parameter gap_extend must be numeric: {gap_extend}')

        if map_limit is not None:
            try:
                self.__map_limit = int(map_limit)
            except ValueError:
                raise RuntimeError(f'ScoreAligner(): Parameter map_limit must be an integer or None (no limit): {map_limit}')
        else:
            self.__map_limit = None

        try:
            self.__jaccard_kmer = int(jaccard_kmer)
        except ValueError:
            raise RuntimeError(f'ScoreAligner(): Parameter jaccard_kmer must be an integer: {jaccard_kmer}')

        # Check argument ranges
        if self.__match <= 0.0:
            raise RuntimeError(f'Match score must be > 0.0: {self.__match}')

        if self.__mismatch >= 0.0:
            raise RuntimeError(f'Mismatch score must be < 0.0: {self.__match}')

        if self.__gap_o > 0.0:
            raise RuntimeError(f'Gap-open score must be <= 0.0: {self.__match}')

        if self.__gap_e > 0.0:
            raise RuntimeError(f'Gap-extend score must be <= 0.0: {self.__match}')

        if self.__map_limit is not None and self.__map_limit < 0:
            raise RuntimeError(f'Map must be >= 0.0 or None (no limit): {self.__map_limit}')

        if self.__jaccard_kmer <= 0:
            raise RuntimeError(f'Jaccard k-mer size must be > 0: {self.jaccard_kmer}')

        return

    def score_align(self, seq_a, seq_b):
        """
        Get max score aligning two sequences. Only returns the score, not the alignment.

        :param seq_a: Subject sequence.
        :param seq_b: Query sequence.

        :return: Maximum alignment score.
        """

        # Compute for convenience
        gap_1bp = self.__gap_o + self.__gap_e

        # Scrub sequences
        seq_a = seq_a.upper().strip()
        seq_b = seq_b.upper().strip()

        # Get length
        len_a = len(seq_a) + 1
        len_b = len(seq_b) + 1

        trace_matrix = [ScoreTraceNode() for i in range(len_a)]
        trace_matrix_last = [ScoreTraceNode() for i in range(len_a)]

        # Max values
        global_max_score = 0.0  # Max score

        # Iterate bases in seq_b
        for j in range(1, len_b):

            # Swap trace arrays
            trace_matrix, trace_matrix_last = trace_matrix_last, trace_matrix

            # Iterate bases in seq_a
            for i in range(1, len_a):

                # Aligned bases
                if seq_a[i - 1] == seq_b[j - 1]:
                    score_max = trace_matrix_last[i - 1].score + self.__match
                    op_code = OP_MATCH

                else:
                    score_max = trace_matrix_last[i - 1].score + self.__mismatch
                    op_code = OP_MISMATCH

                # Gap subject (insertion)
                score_gap_sub = trace_matrix_last[i].score + (
                    self.__gap_e if trace_matrix_last[i].op_code == OP_GAP_SUB else gap_1bp
                )

                if score_gap_sub > score_max:
                    score_max = score_gap_sub
                    op_code = OP_GAP_SUB

                # Gap query (deletion)
                score_gap_qry = trace_matrix[i - 1].score + (
                    self.__gap_e if trace_matrix[i - 1].op_code == OP_GAP_QRY else gap_1bp
                )

                if score_gap_qry > score_max:
                    score_max = score_gap_qry
                    op_code = OP_GAP_SUB

                # Update trace matrix
                if score_max > 0:
                    trace_matrix[i].score = score_max
                    trace_matrix[i].op_code = op_code

                    if op_code == OP_MATCH:

                        # Check for new global max
                        if score_max >= global_max_score:
                            global_max_score = score_max

                else:
                    trace_matrix[i].op_code = OP_NONE
                    trace_matrix[i].score = 0

        return global_max_score

    def match_prop(self, seq_a, seq_b):
        """
        Get the alignment score proportion over the max possible score between two sequences. To ollow tandem
        duplications to map correctly, seq_b is duplicated head-to-tail (seq_b + seq_b) and seq_a is aligned to it.

        The max possible score is achieved if seq_a and seq_b are the same size and seq_a aligns to seq_b + seq_b with
        all seq_a bases matching (function returns 1.0).

        The numerator is the alignment score and the denominator is the size of the larger sequence times the match
        score:

        min(score_align(seq_a, seq_b + seq_b), min_len(seq_a, seq_b) * match) / (max_len(seq_a, seq_b) * match)

        :param seq_a: Subject sequence.
        :param seq_b: Query sequence.

        :return: Alignment proportion with seq_b duplicated head to tail.
        """

        max_len = np.max([len(seq_a), len(seq_b)])
        min_len = np.min([len(seq_a), len(seq_b)])

        if self.__map_limit is None or max_len <= self.__map_limit:
            return min([
                    np.min([self.score_align(seq_a, seq_b + seq_b), min_len * self.__match]) / (max_len * self.__match),
                    1.0
            ])

        elif min_len > self.__jaccard_kmer:
            return jaccard_distance(seq_a, seq_b, self.__jaccard_kmer)

        else:
            return 1 if seq_a.upper() == seq_b.upper() else 0
