"""
General definitions used in multiple places.
"""

def get_svtype_label(svtype):
    """
    Get a label for the SV type.

    :param svset: SV type key.

    :return: SV type label for plots.
    """

    svtype = svtype.strip()

    if not svtype:
        return 'Unknown Variant Type'

    if svtype == 'ins':
        return 'Insertions'

    elif svtype == 'del':
        return 'Deletions'

    elif svtype == 'inv':
        return 'Inversions'

    elif svtype == 'dup':
        return 'Duplications'

    elif svtype == 'insdelinv':
        return 'INS/DEL/INV'

    elif svtype == 'insdel':
        return 'INS/DEL'

    elif svtype == 'snv':
        return 'SNVs'

    else:
        raise RuntimeError('Unrecognized SV type: {}'.format(svtype))

def get_svset_label(svset):
    """
    Get a label for the SV set.

    :param svset: SV set key.

    :return: SV set label for plots.
    """

    svset = svset.strip()

    if not svset:
        return 'Unknown Variant Set'

    if svset == 'notr':
        return 'No TR'

    elif svset == 'intr':
        return 'In TR'

    elif svset == 'notrsd':
        return 'No TR/SD'

    elif svset == 'nosd':
        return 'No SD'

    elif svset == 'insd':
        return 'In SD'

    elif svset == 'alu':
        return 'Alu'

    elif svset == 'aluy':
        return 'AluY'

    elif svset == 'l1':
        return 'L1'

    elif svset == 'sva':
        return 'SVA'

    elif svset == 'all':
        return 'All Variants'

    else:
        return svset.upper()


def get_sample_name(sourcetype, sourcename, sample, prefix_sourcename=False):
    """
    Get the name of a sample by its source type (e.g. "caller" or "sampleset"), source name (e.g. "smrtsv"), and
    sample name.

    :param sourcetype: Source type.
    :param sourcename: Source name.
    :param sample: Sample name.
    :param prefix_sourcename: Prefix sample with the source name. Set this if the sample is being compared with another
        possibly different source.
    """

    # Get sample name
    if sourcetype == 'sampleset':
        sampleset_entry = svpoplib.sampleset.get_config_entry(sourcename, sample, config)

        sample = sampleset_entry['name']

    # Return if no sample prefix
    if not prefix_sourcename:
        return sample

    # Correct for altdup
    source_suffix = None

    if sourcetype == 'caller' and sourcename in SAMPLE_TABLE.index and SAMPLE_TABLE.loc[sourcename, 'TYPE'].squeeze() == 'altdup':
        tok = SAMPLE_TABLE.loc[sourcename, 'DATA'].squeeze().split(',')
        tok = [v.strip() for v in tok]

        if len(tok) >= 2:
            sourcetype = tok[0]
            sourcename = tok[1]
            source_suffix = 'ALT-DUP'

    # Set prefix
    sample_prefix = None

    if sourcetype == 'sampleset':
        sample_prefix = 'Sampleset'

    elif sourcetype == 'callerset':
        sample_prefix = 'Callerset'

    elif sourcetype == 'caller':
        if sourcename in SAMPLE_TABLE.index:
            sample_type = SAMPLE_TABLE.loc[sourcename, 'TYPE'].squeeze()
        else:
            sample_type = sourcename

        caller_name = {
            'pbsv': 'PBSV',
            'deepvariant': 'DeepVariant',
            'bionano': 'Bionano',
            'gatk': 'GATK',
            'sniffles': 'Sniffles',
            'longshot': 'Longshot',
            'vcf': 'VCF Callset',
            'pav': 'PAV',
            'pavbed': 'PAV',
            'svim': 'SVIM',
            'svim-asm': 'SVIM-ASM',
            'sniffles': 'Sniffles'
        }.get(sample_type, sample_type)

        sample_prefix = caller_name

        if source_suffix is not None:
            sample_prefix = f'{sample_prefix} ({source_suffix})'

    # Set sample name and return
    return '{} - {}'.format(sample_prefix, sample)

# Variant type to string
SVTYPE_TO_STRING = collections.defaultdict(lambda: 'Unknown')
SVTYPE_TO_STRING['ins'] = 'Insertion'
SVTYPE_TO_STRING['del'] = 'Deletion'
SVTYPE_TO_STRING['inv'] = 'Inversion'


# CIGAR operations
CIGAR_OP_STR = 'MIDNSHP=X'
CIGAR_INT = {op: code for op, code in zip(list(CIGAR_OP_STR), range(len(CIGAR_OP_STR)))}
CIGAR_CHAR = {code: op for op, code in zip(list(CIGAR_OP_STR), range(len(CIGAR_OP_STR)))}

# Pysam constants
BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8
BAM_CBACK = 9
BAM_NM = 10
