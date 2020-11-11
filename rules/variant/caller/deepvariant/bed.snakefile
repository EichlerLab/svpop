"""
Extract variant callset from DeepVariant VCFs.
"""

# variant_pbsv_bed_fa
#
# Make final BED and FASTA of ins/del/inv sequences.
rule variant_dv_bed_fa:
    input:
        bed='temp/variant/caller/deepvariant-{seq_set}/bed/{sample}/tsv/variants_sv.bed.gz'
    output:
        bed='results/variant/caller/deepvariant-{seq_set}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        fa='results/variant/caller/deepvariant-{seq_set}/fasta/{sample}/all/{filter}/{vartype}_{svtype}.fa.gz'
    wildcard_constraints:
        svtype='ins|del|snv'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        df = df.loc[(df['VARTYPE'] == wildcards.vartype.upper()) & (df['SVTYPE'] == wildcards.svtype.upper())]

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
rule variant_dv_tsv_to_bed:
    input:
        tsv='temp/variant/caller/deepvariant-{seq_set}/bed/{sample}/tsv/variants_sv.tsv.gz'
    output:
        bed=temp('temp/variant/caller/deepvariant-{seq_set}/bed/{sample}/tsv/variants_sv.bed.gz'),
        tsv_filt='results/variant/caller/deepvariant-{seq_set}/bed/{sample}/all/all/filtered/vcf_fail.tsv.gz'
    params:
        cpu=6,
        mem='2g'
    run:

        # Definitions
        COL_RENAME = {
            'CHROM': '#CHROM',
            'REF': 'VCF_REF',
            'POS': 'VCF_POS',
            'ALT': 'VCF_ALT'
        }

        REQUIRED_COLS = {'CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'GT', 'GQ', 'DP', 'AD'}

        # Read
        df = pd.read_csv(input.tsv, sep='\t', header=0, low_memory=False)

        # Remove prefixes added by bcftools (e.g. "# [1]CHROM" to "CHROM", "[2]POS" to "POS"
        df.columns = [re.sub('^#?\s*\[[^\]]+\]', '', col) for col in df.columns]

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
            col_set = {col for col in df.columns if not re.match('.+:.+', col) or col.startswith(wildcards.sample + ':')}

            df = df.loc[[col for col in df.columns if col in col_set]]

        # Remove sample name from columns
        df.columns = [(col if ':' not in col else col.split(':', 1)[1]) for col in df.columns]

        # Check columns and rename
        missing_cols = REQUIRED_COLS - set(df.columns)

        if missing_cols:
            raise RuntimeError('Missing VCF columns: {}'.format(', '.join(sorted(missing_cols))))

        df.columns = [COL_RENAME[col] if col in COL_RENAME else col for col in df.columns]

        # Filter on FILTER column
        df.loc[df['FILTER'] != 'PASS'].to_csv(output.tsv_filt, sep='\t', index=False, compression='gzip')

        df = df.loc[df['FILTER'] == 'PASS']

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

        # Save VCF index
        df['VCF_IDX'] = df.index

        # Set index and arrange columns
        df['ALT'] = df['VCF_ALT']
        df['ID'] = analib.variant.get_variant_id(df)

        df = analib.variant.order_variant_columns(df, tail_cols=('REF', 'ALT', 'QUAL', 'FILTER', 'VCF_POS', 'VCF_REF', 'VCF_ALT', 'VCF_IDX', 'SEQ'))

        # Drop duplicates
        df.drop_duplicates('ID', inplace=True)

        # Write all variants
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# variant_dv_vcf_to_tab
#
# VCF to TSV.
rule variant_dv_vcf_to_tab:
    input:
        vcf=lambda wildcards: variant_global_sample_table_entry('deepvariant', wildcards=wildcards)['DATA']
    output:
        tsv=temp('temp/variant/caller/deepvariant-{seq_set}/bed/{sample}/tsv/variants_{vartype}.tsv.gz')
    shell:
        """
        bcftools query -H -f"%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT][\t%GQ][\t%DP][\t%AD]\n" {input.vcf} | """
        """gzip """
        """> {output.tsv}"""
