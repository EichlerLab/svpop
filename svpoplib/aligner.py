"""
Fast Smith-waterman implementation.
"""

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


class ScoreAligner:
    """
    An engine for aligning two sequences to obtain an alignment score only, no alignment is generated. This is intended
    for determining sequence similarity through an alignment where the alignment itself is not important.

    :param match: Base match score (> 0.0).
    :param mismatch: Base mismatch score (< 0.0).
    :param gap_open: Gap open cost (<= 0.0).
    :param gap_extend: Gap extend cost (<= 0.0).
    """

    def __init__(self, match=2.0, mismatch=-1.0, gap_open=-1.0, gap_extend=-1.0):
        self.__match = match
        self.__mismatch = mismatch
        self.__gap_o = gap_open
        self.__gap_e = gap_extend

        # Check arguments
        if self.__match <= 0.0:
            raise RuntimeError(f'Match score must be > 0.0: {self.__match}')

        if self.__mismatch >= 0.0:
            raise RuntimeError(f'Mismatch score must be < 0.0: {self.__match}')

        if self.__gap_o > 0.0:
            raise RuntimeError(f'Gap-open score must be <= 0.0: {self.__match}')

        if self.__gap_e > 0.0:
            raise RuntimeError(f'Gap-extend score must be <= 0.0: {self.__match}')

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
        seq_a = seq_a.strip()
        seq_b = seq_b.strip()

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
                    self.__gap_e if trace_matrix_last[i].op_code == OP_GAP_QRY else gap_1bp
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
