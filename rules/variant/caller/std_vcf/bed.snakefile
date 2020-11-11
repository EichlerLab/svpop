"""
Extract variant callset from standard VCFs.
"""

CALLER_STDVCF_VARSOURCE_PATTERN = '(gatk|longshot|gtgq)'

STD_VCF_BCFTOOLS_QUERY_STRING = {
    'gatk': '"%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT][\t%GQ][\t%DP][\t%AD]\n"',
    'longshot': '"%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/DP[\t%GT][\t%GQ][\t%DP]\n"',
    'gtgq': '"%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT][\t%GQ]\n"'  # Generic with GT & GQ
}

# variant_vcf_bed_fa
#
# Make final BED and FASTA of ins/del/inv sequences.
rule variant_vcf_bed_fa:
    input:
        bed='temp/variant/caller/{var_source}-{seq_set}/bed/{sample}/tsv/variants_usort.bed.gz'
    output:
        bed='results/variant/caller/{var_source}-{seq_set}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        fa='results/variant/caller/{var_source}-{seq_set}/fasta/{sample}/all/{filter}/{vartype}_{svtype}.fa.gz'
    wildcard_constraints:
        svtype='ins|del|snv',
        var_source=CALLER_STDVCF_VARSOURCE_PATTERN
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        df = df.loc[(df['VARTYPE'] == wildcards.vartype.upper()) & (df['SVTYPE'] == wildcards.svtype.upper())]

        # Sort
        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Write FASTA
        if 'SEQ' in df.columns and df.shape[0] > 0:

            # Write FASTA
            with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
                SeqIO.write(analib.seq.bed_to_seqrecord_iter(df), out_file, 'fasta')

            shell("""samtools faidx {output.fa}""")

            # Remove SEQ column
            del(df['SEQ'])

        else:
            # No sequence output. Make empty files.
            shell("""touch {output.fa}""")

        # Write
        df.to_csv(output.bed, sep='\t', index=False)

# variant_dv_tsv_to_bed
#
# TSV to BED.
rule variant_vcf_tsv_to_bed:
    input:
        tsv='temp/variant/caller/{var_source}-{seq_set}/bed/{sample}/tsv/variants.tsv.gz'
    output:
        bed=temp('temp/variant/caller/{var_source}-{seq_set}/bed/{sample}/tsv/variants_usort.bed.gz'),
        tsv_filt='results/variant/caller/{var_source}-{seq_set}/bed/{sample}/all/all/filtered/vcf_fail.tsv.gz'
    wildcard_constraints:
        var_source=CALLER_STDVCF_VARSOURCE_PATTERN
    params:
        cpu=6,
        mem='1G',
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

        write_header = True  # Write header for this chunk (only true for first chunk)
        filt_header = True   # Write header for this chunk's dropped SV records (only true for first chunk)

        cumulative_id_set = set()  # Set of variants seen in previous chunks for duplicate removal

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

                    # Check columns and rename
                    missing_cols = REQUIRED_COLS - set(df.columns)

                    if missing_cols:
                        raise RuntimeError('Missing VCF columns: {}'.format(', '.join(sorted(missing_cols))))

                    df.columns = [COL_RENAME[col] if col in COL_RENAME else col for col in df.columns]

                    # Set FILTER on QUAL if FILTER is missing
                    df['FILTER'] = df.apply(analib.variant.qual_to_filter, axis=1)

                    # Filter on FILTER column

                    df_filt = df.loc[df['FILTER'] != 'PASS']

                    if df_filt.shape[0] > 0:
                        df_filt.to_csv(filt_file, sep='\t', index=False, header=filt_header)
                        filt_file.flush()

                        filt_header = False

                    df = df.loc[df['FILTER'] == 'PASS']

                    if df.shape[0] == 0:
                        continue

                    # Get SV fields
                    df = pd.concat(
                        [
                            df,
                            analib.pd.apply_parallel(
                                df, analib.variant.vcf_fields_to_seq,
                                n_part=500, n_core=params.cpu,
                                kwds={
                                    'pos_row': 'VCF_POS',
                                    'ref_row': 'VCF_REF',
                                    'alt_row': 'VCF_ALT'
                                }
                            )
                        ],
                        axis=1
                    )

                    df = df.loc[df['AC'] > 0]

                    if df.shape[0] == 0:
                        continue

                    # Save VCF index
                    df['VCF_IDX'] = df.index

                    # Set index and arrange columns
                    df['ALT'] = df['VCF_ALT']
                    df['ID'] = analib.variant.get_variant_id(df)

                    if out_col_order is None:
                        out_col_order = list(analib.variant.order_variant_columns(df, tail_cols=('REF', 'ALT', 'QUAL', 'FILTER', 'VCF_POS', 'VCF_REF', 'VCF_ALT', 'VCF_IDX', 'SEQ')).columns)

                    df = df.loc[:, out_col_order]

                    # Drop duplicates
                    df.drop_duplicates('ID', inplace=True)

                    df = df.loc[df['ID'].apply(lambda val: val not in cumulative_id_set)]

                    cumulative_id_set |= set(df['ID'])

                    # Write
                    df.to_csv(out_file, sep='\t', index=False, header=write_header)
                    out_file.flush()

                    write_header = False


# variant_dv_vcf_to_tab
#
# VCF to TAB file.
rule variant_vcf_to_tab:
    input:
        vcf=lambda wildcards: variant_global_sample_table_entry(wildcards.var_source, wildcards=wildcards)['DATA']
    output:
        tsv=temp('temp/variant/caller/{var_source}-{seq_set}/bed/{sample}/tsv/variants.tsv.gz')
    params:
        query_string=lambda wildcards: STD_VCF_BCFTOOLS_QUERY_STRING[wildcards.var_source]
    wildcard_constraints:
        var_source=CALLER_STDVCF_VARSOURCE_PATTERN
    shell:
        """bcftools query -H -f{params.query_string} {input.vcf} | """
        """gzip """
        """> {output.tsv}"""
