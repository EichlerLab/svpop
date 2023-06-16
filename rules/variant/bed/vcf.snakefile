"""
Extract variant callset from standard VCFs.
"""

CALLER_VCF_STD_FIELDS = {
    'gatk': {
        'info': [],
        'format': ['GT', 'GQ', 'DP', 'AD']
    },
    'longshot': {
        'info': [],
        'format': ['GT', 'GQ', 'DP']
    },
    'deepvariant': {
        'info': [],
        'format': ['GT', 'GQ', 'DP', 'AD', 'VAF', 'PL']
    },
    'sniffles': {
        'info': ['SVTYPE', 'SVLEN', 'END'],
        'format': ['GT', 'GQ', 'DR', 'DV']
    },
    'sniffles2': {
        'info': ['SVTYPE', 'SVLEN', 'END'],
        'format': ['GT', 'GQ', 'DR', 'DV']
    },
    'svimasm0': {
        'info': ['SVTYPE', 'END', 'SVLEN'],
        'format': ['GT']
    },
    'svimasm': {
        'info': ['SVTYPE', 'END', 'SVLEN'],
        'format': ['GT']
    },
    'clair3': {
        'info': ['P', 'F'],
        'format': ['GT', 'GQ', 'DP', 'AD', 'PL', 'AF']
    },
    'cutesv': {
        'info': ['PRECISE', 'IMPRECISE', 'SVTYPE', 'SVLEN', 'END', 'CIPOS', 'CILEN'],
        'format': ['DR', 'DV', 'PL', 'GQ']
    }
}

VARIANT_BED_VCF_TYPE_PATTERN = '(vcf|{})'.format('|'.join(CALLER_VCF_STD_FIELDS.keys()))

DEFAULT_FILTER_PASS_SET = {'PASS', '.'}

#
# Helper functions
#

def variant_bed_vcf_get_param_dict(wildcards):
    """
    Get dictionary of parameters

    :param wildcards: Rule wildcards.
    :return:
    """

    # Get sample entry
    sample_entry = svpoplib.rules.sample_table_entry(
        wildcards.sourcename,
        SAMPLE_TABLE,
        wildcards=wildcards,
        caller_type=wildcards.callertype
    )

    # Initialize INFO and FORMAT field lists
    param_dict = {
        'info_list': list(),    # List of INFO fields to include
        'format_list': list(),  # List of FORMAT fields to gather
        'pass': None      # Set of passing FILTER values. Set to default value if None
    }

    keyword_set = list()  # Keywords already processed

    # Add type keyword
    if wildcards.callertype != 'vcf':

        # Check for recognized caller type
        if wildcards.callertype not in CALLER_VCF_STD_FIELDS.keys():
            raise RuntimeError('Error processing entry "{}" (type "{}", sample "{}"): Unrecognized type "{}"'.format(
                wildcards.sourcename, wildcards.callertype, wildcards.sample, wildcards.callertype
            ))

        # Add to INFO and FORMAT
        param_dict['info_list'] = CALLER_VCF_STD_FIELDS[wildcards.callertype]['info']
        param_dict['format_list'] = CALLER_VCF_STD_FIELDS[wildcards.callertype]['format']

        # Add filter
        if 'pass' in CALLER_VCF_STD_FIELDS[wildcards.callertype]:
            param_dict['pass'] = CALLER_VCF_STD_FIELDS[wildcards.callertype]['pass']

    # Process options
    param_string = sample_entry['PARAM_STRING']

    if not (param_string is None or pd.isnull(param_string)):
        param_string = param_string.strip()
    else:
        param_string = ''

    if param_string:
        for avp in param_string.split(';'):

            avp = avp.strip()

            # Process keywords
            if avp in CALLER_VCF_STD_FIELDS.keys():

                if avp not in keyword_set:
                    std_entry = CALLER_VCF_STD_FIELDS[avp]

                    param_dict['info_list'] += [val for val in std_entry['info'] if val not in param_dict['info_list']]
                    param_dict['format_list'] += [val for val in std_entry['format'] if val not in param_dict['format_list']]
                    keyword_set.append(avp)

                continue

            if '=' in avp:
                attrib, value = avp.split('=', 1)
                value = value.strip()
            else:
                attrib, value = avp, None

            # Split and process attribute/value
            attrib = attrib.strip()

            if attrib == '':
                # Ignore empty entries
                if value is not None and value != '':
                    raise RuntimeError(f'Found value with an empty attribute (blank before "=": {avp}')

            elif attrib == 'info':
                if value is None:
                    raise RuntimeError(f'Missing value for "info" in callset configuration: {avp}')

                param_dict['info_list'] += [val for val in value.split(',') if val not in param_dict['info_list']]

            elif attrib == 'format':
                if value is None:
                    raise RuntimeError(f'Missing value for "format" in callset configuration: {avp}')

                param_dict['format_list'] += [val for val in value.split(',') if val not in param_dict['format_list']]

            elif attrib == 'pass':
                # Append to pass filters if "=" is in the avp
                # Accept all filters if "=" is not in the avp (override existing list)

                if value is not None:
                    if param_dict['pass'] is None:
                        param_dict['pass'] = set()

                    new_filter_set = set([filter_val.strip() for filter_val in value.split(',') if filter_val.strip()])
                    param_dict['pass'] |= new_filter_set

                else:
                    param_dict['pass'] = set()

            else:
                # Add as a list of values, one for each attribute appearance.

                if 'attrib' not in param_dict:
                    param_dict['attrib'] = list()

                param_dict['attrib'].append(value)

    # Set default values
    if param_dict['pass'] is None:
        param_dict['pass'] = DEFAULT_FILTER_PASS_SET

    return param_dict


def variant_bed_vcf_get_bcftools_query(wildcards):
    """
    Get the bcftools query string for this entry.

    :param wildcards: Rule wildcards.

    :return: bcftools format string.
    """

    param_dict = variant_bed_vcf_get_param_dict(wildcards)

    # Construct query string for bcftools
    query_string = r'"%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER'

    for element in param_dict['info_list']:
        query_string += f'\\t%INFO/{element}'

    for element in param_dict['format_list']:
        query_string += f'[\\t%{element}]'

    query_string += '\\n"'

    # Return formatted query string
    return query_string

# def variant_bed_vcf_fix_sniffles2(df):
#     """
#     Correct poor VCF formatting in Sniffles2.
#     """
#
#     df_list = list()
#
#     # Transform records to symbolic ALTs where SEQ was stuck into the ALT column by Sniffles2.
#     for index, row in df.iterrows():
#         if row['VCF_REF'].upper() == 'N':
#             if re.match('^(?![Nn])[ACGTNacgtn]+(?<![Nn])$', row['VCF_ALT']) is not None:
#                 row['SEQ'] = row['VCF_ALT']
#                 row['VCF_REF'] = '.'
#                 row['VCF_ALT'] = '<INS>'
#
#         elif row['VCF_ALT'].upper() == 'N':
#             if re.match('^(?![Nn])[ACGTNacgtn]+(?<![Nn])$', row['VCF_REF']) is not None:
#                 row['SEQ'] = row['VCF_REF']
#                 row['VCF_REF'] = '.'
#                 row['VCF_ALT'] = '<DEL>'
#
#         else:
#             row['SEQ'] = np.nan
#
#         df_list.append(row)
#
#     return pd.concat(df_list, axis=1).T


CALLER_CALLBACK_PRE_BED = {
#     'sniffles2': variant_bed_vcf_fix_sniffles2
}

#
# Rules
#

# variant_vcf_bed_fa
#
# Make final BED and FASTA of ins/del/inv sequences.
rule variant_vcf_bed_fa:
    input:
        bed='temp/variant/caller/{callertype}/{sourcename}/{sample}/tsv/variants.bed.gz'
    output:
        bed=temp('temp/variant/caller/{callertype}/{sourcename}/{sample}/all/all/bed/pre_filter/{vartype}_{svtype}.bed.gz'),
        fa=temp('temp/variant/caller/{callertype}/{sourcename}/{sample}/all/all/bed/pre_filter/fa/{vartype}_{svtype}.fa.gz')
    wildcard_constraints:
        svtype='ins|del|inv|snv|dup|sub|rgn',
        callertype=VARIANT_BED_VCF_TYPE_PATTERN
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        df = df.loc[(df['VARTYPE'] == wildcards.vartype.upper()) & (df['SVTYPE'] == wildcards.svtype.upper())]

        # Sort
        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Write FASTA
        if 'SEQ' in df.columns and df.shape[0] > 0:
            df_seq = df[['ID', 'SEQ']].copy()

            df_seq = df_seq.loc[~ pd.isnull(df_seq['SEQ'])]

            if df_seq.shape[0] == 0:
                df_seq = None
        else:
            df_seq = None

        if df_seq is not None:

            # Write FASTA
            with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
                SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(df_seq), out_file, 'fasta')

            shell("""samtools faidx {output.fa}""")

            # Remove SEQ column
            del(df_seq)

        else:
            # No sequence output. Make empty files.
            shell("""touch {output.fa}""")

        if 'SEQ' in df.columns:
            del(df['SEQ'])

        # Write
        df.to_csv(output.bed, sep='\t', index=False)

# variant_dv_tsv_to_bed
#
# TSV to BED.
rule variant_bed_vcf_tsv_to_bed:
    input:
        tsv='temp/variant/caller/{callertype}/{sourcename}/{sample}/tsv/bcftools_table.tsv.gz'
    output:
        bed=temp('temp/variant/caller/{callertype}/{sourcename}/{sample}/tsv/variants.bed.gz'),
        tsv_filt='results/variant/caller/{sourcename}/{sample}/all/all/bed/filtered/vcf_fail_{callertype}.tsv.gz'
    wildcard_constraints:
        callertype=VARIANT_BED_VCF_TYPE_PATTERN
    params:
        mem='6000',
        chunk_size=20000  # DataFrame chunk size
    threads: 6
    run:

        param_dict = variant_bed_vcf_get_param_dict(wildcards)

        with gzip.open(output.bed, 'wt') as bed_file:
            with gzip.open(output.tsv_filt, 'wt') as filt_file:
                for df in svpoplib.variant.vcf_tsv_to_bed(
                    tsv_in=input.tsv,
                    sample=wildcards.sample,
                    bed_file=bed_file, filt_file=filt_file,
                    chunk_size=params.chunk_size,
                    threads=threads,
                    callback_pre_bed=CALLER_CALLBACK_PRE_BED.get(wildcards.callertype, None),
                    filter_pass_set=param_dict['pass']
                ):
                    pass  # Iterate through all records



# variant_dv_vcf_to_tab
#
# VCF to TAB file.
rule variant_bed_vcf_to_tsv:
    input:
        vcf=lambda wildcards: svpoplib.rules.sample_table_entry(wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type=wildcards.callertype)['DATA']
    output:
        tsv=temp('temp/variant/caller/{callertype}/{sourcename}/{sample}/tsv/bcftools_table.tsv.gz')
    params:
        query_string=variant_bed_vcf_get_bcftools_query
    wildcard_constraints:
        callertype=VARIANT_BED_VCF_TYPE_PATTERN
    shell:
        """bcftools query -H -f{params.query_string} {input.vcf} | """
        """gzip """
        """> {output.tsv}"""
