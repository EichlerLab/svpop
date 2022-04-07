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
    }
}

VARIANT_BED_VCF_TYPE_PATTERN = '(vcf|{})'.format('|'.join(CALLER_VCF_STD_FIELDS.keys()))

def variant_bed_vcf_get_bcftools_query(wildcards):
    """
    Get the bcftools query string for this entry.

    :param wildcards: Rule wildcards.

    :return: bcftools format string.
    """

    # Get sample entry
    sample_entry = svpoplib.rules.sample_table_entry(
        wildcards.sourcename,
        SAMPLE_TABLE,
        wildcards=wildcards,
        type=wildcards.callertype
    )

    # Initialize INFO and FORMAT field lists
    info_list = list()    # List of INFO fields to include
    format_list = list()  # List of FORMAT fields to gather
    keyword_set = list()  # Keywords already processed

    # Add type keyword
    if wildcards.callertype != 'vcf':

        # Check for recognized caller type
        if wildcards.callertype not in CALLER_VCF_STD_FIELDS.keys():
            raise RunimeError('Error processing entry "{}" (type "{}", sample "{}"): Unrecognized type "{}"'.format(
                wildcards.sourcename, wildcards.callertype, wildcards.sample, wildcards.callertype
            ))

        # Add to INFO and FORMAT
        info_list = CALLER_VCF_STD_FIELDS[wildcards.callertype]['info']
        format_list = CALLER_VCF_STD_FIELDS[wildcards.callertype]['format']

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

                    info_list += [val for val in std_entry['info'] if val not in info_list]
                    format_list += [val for val in std_entry['format'] if val not in format_list]
                    keyword_set.add(avp)

                continue

            if '=' not in avp:
                raise RunimeError('Error processing entry "{}" (type "{}", sample "{}"): Unrecognized bare keyword "{}" (no "=") in parameter string "{}"'.format(
                    wildcards.sourcename, wildcards.callertype, wildcards.sample, avp, param_string
                ))

            # Split and process attribute/value
            attrib, value = avp.split('=', 1)

            attrib = attrib.strip()
            value = value.strip()

            if attrib == 'info':
                info_list += [val for val in value.split(',') if val not in info_list]

            elif attrib == 'format':
                format_list += [val for val in value.split(',') if val not in format_list]

            else:
                pass
                # Ignore, let the VCF parser handle unknown attribute values
                # raise RunimeError('Error processing entry "{}" (type "{}", sample "{}"): Unrecognized keyword "{}" in parameter string "{}"'.format(
                #     wildcards.sourcename, wildcards.callertype, wildcards.sample, attrib, param_string
                # ))

    # Construct query string for bcftools
    query_string = r'"%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER'

    for element in info_list:
        query_string += f'\\t%INFO/{element}'

    for element in format_list:
        query_string += f'[\\t%{element}]'

    query_string += '\\n"'

    # Return formatted query string
    return query_string


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
            del(df['SEQ'])

        else:
            # No sequence output. Make empty files.
            shell("""touch {output.fa}""")

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
        cpu=6,
        mem='6000',
        df_chunk_size=5000  # DataFrame chunk size
    run:

        # Definitions
        COL_RENAME = {
            'CHROM': '#CHROM',
            'REF': 'VCF_REF',
            'POS': 'VCF_POS',
            'ALT': 'VCF_ALT'
        }

        REQUIRED_COLS = {'CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'GT'}

        # Read (as iterator)
        df_iter = pd.read_csv(input.tsv, sep='\t', header=0, low_memory=False, iterator=True, chunksize=params.df_chunk_size)

        # Process input in chunks
        col_names = None      # Input column names (after first regex transform)
        out_col_order = None  # Output column order

        id_set = set()  # Set of IDs already seen (used for ensuring IDs are unique)

        write_header = True  # Write header for this chunk (only true for first chunk)
        filt_header = True   # Write header for this chunk's dropped SV records (only true for first chunk)

        with gzip.open(output.bed, 'wt') as out_file:
            with gzip.open(output.tsv_filt, 'wt') as filt_file:
                for df in df_iter:

                    # Remove prefixes added by bcftools (e.g. "# [1]CHROM" to "CHROM", "[2]POS" to "POS"
                    df.columns = [re.sub('^#?\s*\[[^\]]+\]', '', col) for col in df.columns]

                    if col_names is None:

                        # Check samples. If single-sample VCF, accept it regardless of the sample name. If multi-sample, must contain
                        # wildcards.sample as a sample in the VCF.
                        samples = {col.split(':', 1)[0] for col in df.columns if re.match('.+:.+', col)}

                        if len(samples) > 1:

                            # Multi-samples and none match wildcards.sample
                            if wildcards.sample not in samples:
                                raise RuntimeError(
                                    'Detected multiple samples in VCF with no samples matching wildcards.sample: {}'.format(
                                        ', '.join(sorted(samples))
                                    )
                                )

                            # Subset to non-sample columns (CHROM, POS, REF, ALT, etc) and sample columns (GT, etc)
                            col_names = [col for col in df.columns if not re.match('.+:.+', col) or col.startswith(wildcards.sample + ':')]

                        else:
                            col_names = df.columns

                    # Subset to required columns
                    df = df[col_names].copy()

                    # Remove sample name from columns
                    df.columns = [(col if ':' not in col else col.split(':', 1)[1]) for col in df.columns]

                    if len(samples) > 1 and 'GT' not in df.columns:
                        raise RuntimeError('Multi-sample VCF missing GT column')

                    # Check columns and rename
                    missing_cols = REQUIRED_COLS - set(df.columns)

                    if missing_cols:
                        raise RuntimeError('Missing VCF columns: {}'.format(', '.join(sorted(missing_cols))))

                    df.columns = [COL_RENAME[col] if col in COL_RENAME else col for col in df.columns]

                    # Set FILTER on QUAL if FILTER is missing
                    df['FILTER'] = df.apply(svpoplib.variant.qual_to_filter, axis=1)

                    df_pass_filter = (df['FILTER'] == 'PASS') | (df['FILTER'] == '.')

                    # Filter on FILTER column

                    df_filt = df.loc[~ df_pass_filter]

                    if df_filt.shape[0] > 0:
                        df_filt.to_csv(filt_file, sep='\t', index=False, header=filt_header)
                        filt_file.flush()

                        filt_header = False

                    df = df.loc[df_pass_filter]

                    if df.shape[0] == 0:
                        continue

                    # Get SV fields
                    df_var_fields = svpoplib.pd.apply_parallel(
                        df, svpoplib.variant.vcf_fields_to_seq,
                        n_part=500, n_core=params.cpu,
                        kwds={
                            'pos_row': 'VCF_POS',
                            'ref_row': 'VCF_REF',
                            'alt_row': 'VCF_ALT'
                        }
                    )

                    df = df[[col for col in df.columns if col not in df_var_fields.columns]]

                    df = pd.concat([df, df_var_fields], axis=1)

                    if 'GT' in df.columns and np.min(df['AC'] > 0):
                        df = df.loc[df['AC'] > 0]

                    if df.shape[0] == 0:
                        continue

                    # Save VCF index
                    df['VCF_IDX'] = df.index

                    # Set index and arrange columns
                    df['ALT'] = df['VCF_ALT']
                    df['ID'] = svpoplib.variant.get_variant_id(df)

                    if out_col_order is None:
                        out_col_order = list(svpoplib.variant.order_variant_columns(df, tail_cols=('QUAL', 'FILTER', 'VCF_POS', 'VCF_REF', 'VCF_ALT', 'VCF_IDX', 'SEQ')).columns)

                    df = df.loc[:, out_col_order]

                    # Rename duplicates
                    df['ID'] = svpoplib.variant.version_id(df['ID'], id_set)
                    id_set |= set(df['ID'])

                    # Write
                    df.to_csv(out_file, sep='\t', index=False, header=write_header)
                    out_file.flush()

                    write_header = False


# variant_dv_vcf_to_tab
#
# VCF to TAB file.
rule variant_bed_vcf_to_tsv:
    input:
        vcf=lambda wildcards: svpoplib.rules.sample_table_entry(wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, type=wildcards.callertype)['DATA']
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
