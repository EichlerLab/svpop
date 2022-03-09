"""
Parameters for svmerge.
"""

import re

import svpoplib

ALIGN_PARAM_FIELD_LIST = [
    ('SCORE', 0.8, float),      # Minimum proportion of the maximum possible score (all bases aligned)
    ('MATCH', 2.0, float),      # Match score
    ('MISMATCH', -1.0, float),  # Mismatch score
    ('OPEN', -1.0, float),      # Gap open score
    ('EXTEND', -0.25, float),   # Gap extend score
    ('LIMIT', 20000, int),      # Fall back to Jaccard index after this size.
    ('KSIZE', 9, int)           # K-mer size for Jaccard index (Using Jasmine default)
]

# Map parameter string position to parameter name
ALIGN_PARAM_KEY = {
    i: ALIGN_PARAM_FIELD_LIST[i][0] for i in range(len(ALIGN_PARAM_FIELD_LIST))
}

# Default parameter values keyed by parameter name
ALIGN_PARAM_DEFAULT = {
    ALIGN_PARAM_FIELD_LIST[i][0]: ALIGN_PARAM_FIELD_LIST[i][1] for i in range(len(ALIGN_PARAM_FIELD_LIST))
}

# Data type for each parameter (keyed by parameter name)
ALIGN_PARAM_TYPE = {
    ALIGN_PARAM_FIELD_LIST[i][0]: ALIGN_PARAM_FIELD_LIST[i][2] for i in range(len(ALIGN_PARAM_FIELD_LIST))
}


class SVMergeParams:
    """
    Parse parameters into a parameter object. Throws exceptions if there are any problems with the parameters.

    :param merge_params: Parameter string (strategy prefix removed).
    :param strategy: Full merge strategy for error reporting.

    :return: An object with fields set.
    """
    
    def __init__(self, merge_params, strategy=None):

        if strategy is None:
            strategy = merge_params

        # General parameters
        self.ro_min = None        # Reciprocal overlap threshold
        self.szro_min = None      # Size reciprocal-overlap threshold
        self.offset_max = None    # Max variant offset
        self.match_ref = False    # Match REF (exact match) if True
        self.match_alt = False    # Match ALT (exact match) if True

        # Sequence match parameters
        self.match_seq = False        # If True, match sequences by aligning and applying a match threshold
        self.aligner = None           # Configured alignger with match/mismatch/gap parameters set
        self.align_match_prop = None  # Match proprotion threshold (matching bases / sequence length)
        self.align_param = None       # Alignment parameters used to configure aligner

        # Read sequence
        self.read_seq = False  # If set, a parameter requires sequence-resolution as a SEQ column (flag to load SEQ from FASTA)

        # Strategy and params
        self.strategy = strategy
        self.merge_params = merge_params

        # Check parameters
        if merge_params is None:
            raise RuntimeError(f'Cannot merge (strategy={strategy}) with parameters: None')

        # Split and parse
        for param_element in merge_params.split(':'):

            # Tokenize
            param_element = param_element.strip()

            if not param_element:
                continue

            param_tok = re.split('\s*=\s*', param_element, 1)

            key = param_tok[0].lower()

            if len(param_tok) > 1:
                val = param_tok[1]
            else:
                val = None

            # Process key
            if key == 'ro':

                if val is None:
                    raise ValueError(f'Missing value for parameter "ro" (e.g. "ro=50"): {merge_params}')

                if val != 'any':
                    self.ro_min = int(val.strip()) / 100

                    if self.ro_min < 0 or self.ro_min > 1:
                        raise ValueError(
                            f'Overlap length (ro) must be between 0 and 100 (inclusive): {param_element}'
                        )

                else:
                    raise RuntimeError('RO "any" is not yet implemented')

            elif key == 'szro':
                if val is None:
                    raise ValueError(f'Missing value for parameter "szro" (e.g. "szro=50"): {merge_params}')

                if val != 'any':
                    self.szro_min = int(val.strip()) / 100

                    if self.szro_min < 0 or self.szro_min > 1:
                        raise ValueError(
                            f'Overlap length (szro) must be between 0 and 100 (inclusive): {param_element}'
                        )

                else:
                    self.szro_min = None

            elif key == 'offset':
                if val is None:
                    raise ValueError(f'Missing value for parameter "offset" (maxiumum offset, e.g. "offset=2000"): {merge_params}')

                if val != 'any':
                    self.offset_max = int(val.strip())

                    if self.offset_max < 0:
                        raise RuntimeError(f'Maximum offset (offset parameter) may not be negative: {merge_params}')

                else:
                    self.offset_max = None

            elif key == 'refalt':
                if val is not None:
                    raise RuntimeError(f'Match-REF/ALT (refalt) should not have an argument: {merge_params}')

                self.match_ref = True
                self.match_alt = True

            elif key == 'ref':
                if val is not None:
                    raise RuntimeError(f'Match-REF (ref) should not have an argument: {merge_params}')

                self.match_ref = True

            elif key == 'alt':
                if val is not None:
                    raise RuntimeError(f'Match-ALT (alt) should not have an argument: {merge_params}')

                self.match_alt = True

            elif key == 'match':
                # Align arguments:
                # 0: Min score proportion
                # 1: Match
                # 2: Mismatch
                # 3: Gap open
                # 4: Gap extend
                # 5: Map limit ("NA" or "UNLIMITED" sets no limit). Fall back to Jaccard index after this limit.
                # 6: Jaccard k-mer size

                # Default parameters
                self.align_param = ALIGN_PARAM_DEFAULT.copy()

                # Split val into tokens
                if not val:
                    val = ''

                val = val.strip()

                val_split = val.split(',')

                # Process each value token
                keyword_param = False  # Set to true for the first keyword parameter (positional parameters are not allowed after)

                for i in range(len(val_split)):

                    # Check positional limit
                    if len(val_split) > len(self.align_param) and not keyword_param:
                        raise RuntimeError('Positional parameters in "match" argument count {} exceeds max: {}: {}'.format(len(val_split), len(self.align_param), merge_params))

                    if '=' in val_split[i]:
                        field_name, tok = val_split[i].split('=', 1)

                        field_name = field_name.strip().upper()
                        tok = tok.strip()

                        keyword_param = True

                    else:
                        if keyword_param:
                            raise RuntimeError('Found positional parameter "{}" after keyword parameters (key=...) in parameter string: {}'.format(
                                val_split[i], merge_params
                            ))

                        tok = val_split[i].strip()
                        field_name = ALIGN_PARAM_KEY[i]

                    # Get parameter type
                    param_type = ALIGN_PARAM_TYPE[field_name]

                    if not tok:
                        continue  # Leave default unchanged

                    # Get token value
                    if field_name == 'LIMIT':

                        # Stop at unlimited
                        if tok.lower in {'unlimited'}:
                            tok_val = None

                        else:
                            # Process multipiler
                            tok_val = tok.lower()

                            if tok_val.endswith('k'):
                                multiplier = int(1e3)
                                tok_val = tok_val[:-1]

                            elif tok_val.endswith('m'):
                                multiplier = int(1e6)
                                tok_val = tok_val[:-1]

                            elif tok_val.endswith('g'):
                                multiplier = int(1e9)
                                tok_val = tok_val[:-1]

                            else:
                                multiplier = 1

                            # To int
                            try:
                                tok_val = int(tok_val) * multiplier
                            except ValueError:
                                raise RuntimeError(f'Alignment parameter {i} in "match" type mismatch: Expected int: {tok_val}')

                            # Check
                            if tok_val < 0:
                                raise RuntimeError(f'Alignment parameter {i} ({field_name}) in "match" must not be negative: {tok}')

                    else:
                        try:
                            tok_val = param_type(tok)

                        except ValueError:
                            raise RuntimeError(f'Alignment parameter {i} in "match" type mismatch: Expected {param_type}: {tok}')

                        if field_name == 'SCORE':
                            if tok_val <= 0.0 or tok_val > 1:
                                raise RuntimeError(f'Alignment parameter {i} ({field_name}) in "match" must be between 0 (exclusive) and 1 (inclusive): {tok}')

                        elif field_name in {'MATCH', 'KSIZE'}:
                            if tok_val <= 0:
                                raise RuntimeError(f'Alignment parameter {i} ({field_name}) in "match" must be positive: {tok}')

                        elif field_name in {'MISMATCH', 'OPEN', 'EXTEND'}:
                            if tok_val > 0.0:
                                raise RuntimeError(f'Alignment parameter {i} ({field_name}) in "match" must not be positive: {tok}')

                        else:
                            raise RuntimeError(f'PROGRAM BUG: Unrecognized alignment parameter {i} ({field_name}): {tok}')

                    # Assign
                    self.align_param[field_name] = tok_val

            else:
                raise ValueError(f'Unknown parameter token: {key}')

        # Check parameters
        if self.szro_min is not None and self.offset_max is None:
            raise RuntimeError('Parameters "szro" was specified without "offset"')

        # Get merge size threshold
        if self.ro_min is None and self.szro_min is not None:
            self.ro_min = self.szro_min

        # Set match_seq and aligner
        if self.align_param is not None:
            self.match_seq = True

            self.aligner = svpoplib.aligner.ScoreAligner(
                match=self.align_param['MATCH'],
                mismatch=self.align_param['MISMATCH'],
                gap_open=self.align_param['OPEN'],
                gap_extend=self.align_param['EXTEND']
            )

            self.align_match_prop = self.align_param['SCORE']

        else:
            self.match_seq = False
            self.aligner = None
            self.align_match_prop = None

        # Set read_seq
        self.read_seq = self.match_seq  # Future: read_seq may be set for other reasons, keep as a separate flag

    def format_string(self, multi_line=True):
        """
        Get a formatted printable string describing the parameters.

        :param multi_line: If True, print on multiple lines, otherwise, comma-separated on one line.

        :return: Formatted string.
        """

        fmt_dict = dict(self.__dict__)

        if multi_line:
            fmt_dict['sep'] = '\n'
            fmt_dict['prefix'] = '\t'
        else:
            fmt_dict['sep'] = ', '
            fmt_dict['prefix'] = ''

        if fmt_dict['align_param'] is not None:
            key_list = [tok[0] for tok in ALIGN_PARAM_FIELD_LIST]
            fmt_dict['align_param'] = '{' + ', '.join(['{}: {}'.format(key, fmt_dict['align_param'][key]) for key in key_list]) + '}'

        return (
            """SVMergeParams{{{sep}"""
            """{prefix}ro_min: {ro_min}{sep}"""
            """{prefix}szro_min: {szro_min}{sep}"""
            """{prefix}offset_max: {offset_max}{sep}"""
            """{prefix}match_ref: {match_ref}{sep}"""
            """{prefix}match_alt: {match_alt}{sep}"""
            """{prefix}match_seq: {match_seq}{sep}"""
            """{prefix}aligner: {aligner}{sep}"""
            """{prefix}align_match_prop: {align_match_prop}{sep}"""
            """{prefix}align_param: {align_param}{sep}"""
            """{prefix}read_seq: {read_seq}{sep}"""
            """{prefix}strategy: {strategy}{sep}"""
            """{prefix}merge_params: {merge_params}{sep}"""
            """}}"""
        ).format(**fmt_dict)
