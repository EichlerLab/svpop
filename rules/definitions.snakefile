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

    elif svtype == 'all':
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

def get_seq_set_name(seq_set):
    """
    Get name of a sequence set (after "-" is sourcename).

    :param seq_set: Sequence set name.

    :return: Formatted sequence set name or `seq_set` if Unknown.
    """

    if seq_set is None:
        return None

    return {
        'hifi': 'HiFi',
        'clr': 'CLR',
        'illumina': 'Illumina',
        'ont': 'ONT',
        'solve': 'Solve',  # Bionano
        'rvp': 'RVP' # Bionano
    }.get(seq_set, seq_set)

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

    if sourcetype == 'varset' and sample == 'all':
        sample = 'All'

    # Return if no sample prefix
    if not prefix_sourcename:
        return sample

    # Set prefix
    sample_prefix = None

    if sourcetype == 'sampleset':
        sample_prefix = 'Sampleset'

    elif sourcetype == 'callerset':
        sample_prefix = 'Callerset'

    elif sourcetype == 'caller':

        tok = sourcename.split('-', 1)

        caller_name = tok[0]
        seq_set_name = get_seq_set_name(tok[1] if len(tok) > 1 else None)

        if sourcename.startswith('extern-'):
            sample_prefix = svpoplib.variant_extern.get_config(sourcename[len('extern-'):], config)['name']

        else:
            caller_name = {
                'smrtsv': 'SMRT-SV',
                'pbsv': 'PBSV',
                'deepvariant': 'DeepVariant',
                'bionano': 'Bionano',
                'gatk': 'GATK',  # from std_vcf
                'sniffles': 'Sniffles',
                'longshot': 'Longshot', # from std_vcf
                'gtgq': 'Callset',  # Unknown callset from std_vcf
            }.get(caller_name, caller_name)

            if seq_set_name is not None:

                seq_set_name = {
                    'hifi': 'HiFi',
                    'clr': 'CLR',
                    'illumina': 'Illumina',
                    'ont': 'ONT',
                    'solve': 'Solve',  # Bionano
                    'rvp': 'RVP' # Bionano
                }.get(seq_set_name, seq_set_name)

                sample_prefix = '{} {}'.format(caller_name, seq_set_name)

            else:
                sample_prefix = caller_name

    elif sourcetype == 'varset':
        varset_entry = svpoplib.varset.get_config_entry(sourcename, None, config)

        sample_prefix = varset_entry['name']

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
