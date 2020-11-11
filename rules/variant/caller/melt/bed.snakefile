"""
Extract variant callset from DeepVariant VCFs.
"""

# variant_melt_bed_fa
#
# Make final BED and FASTA of ins/del/inv sequences.
rule variant_melt_bed_fa:
    input:
        bed='temp/variant/caller/melt-{seq_set}/bed/{sample}/all/{filter}/byref/sv_{svtype}.bed'
    output:
        bed='results/variant/caller/melt-{seq_set}/bed/{sample}/all/{filter}/byref/sv_{svtype}.bed',
        fa='results/variant/caller/melt-{seq_set}/fasta/{sample}/all/{filter}/sv_{svtype}.fa.gz'
    wildcard_constraints:
        svtype='ins|del'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        df = df.loc[df['SVTYPE'] == wildcards.svtype.upper()]

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


# variant_pbsv_bed_filter_region
#
# Apply a BED filter to SVs.
rule variant_melt_bed_filter_region:
    input:
        bed='temp/variant/caller/melt-{seq_set}/bed/{sample}/tsv/variants_sv.bed.gz',
        filter=lambda wildcards: analib.variant.get_filter_bed(wildcards.filter, UCSC_REF_NAME, config, SVPOP_DIR)
    output:
        bed=temp('temp/variant/caller/melt-{seq_set}/bed/{sample}/all/{filter}/byref/sv_{svtype}.bed'),
        bed_filt='variant/caller/melt-{seq_set}/bed/{sample}/all/{filter}/byref/removed/sv_{svtype}_filtered.bed'
    wildcard_constraints:
        filter='((?!all)[^\/]*)|all.+',
        svtype='ins|del|inv|dup'
    run:

        if wildcards.filter != 'all':
            shell(
                """bedtools intersect -wa -v -sorted -a {input.bed} -b {input.filter} -header > {output.bed}; """
                """bedtools intersect -wa -u -sorted -a {input.bed} -b {input.filter} -header > {output.bed_filt}; """
            )
        else:
            shell(
                """cp {input.bed} {output.bed}; """
                """touch {output.bed_filt}; """
            )


# variant_melt_tsv_to_bed
#
# TSV to BED.
rule variant_melt_tsv_to_bed:
    input:
        tsv='temp/variant/caller/melt-{seq_set}/bed/{sample}/tsv/variants_sv.tsv.gz'
    output:
        bed=temp('temp/variant/caller/melt-{seq_set}/bed/{sample}/tsv/variants_sv.bed.gz'),
        tsv_filt='results/variant/caller/melt-{seq_set}/bed/{sample}/all/all/filtered/vcf_fail.tsv.gz'
    params:
        cpu=6,
        mem='2g'
    run:

        # Definitions
        REQUIRED_COLS = {'CHROM', 'POS', 'ALT', 'QUAL', 'FILTER', 'SVTYPE', 'SVLEN', 'GT'}

        COL_RENAME = {
            'CHROM': '#CHROM',
            'SVTYPE': 'MECLASS'
        }

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

            df = df[[col for col in df.columns if col in col_set]]

        # Remove sample name from columns
        df.columns = [(col if ':' not in col else col.split(':', 1)[1]) for col in df.columns]

        # Check columns and rename
        missing_cols = REQUIRED_COLS - set(df.columns)

        if missing_cols:
            raise RuntimeError('Missing VCF columns: {}'.format(', '.join(sorted(missing_cols))))

        df.columns = [COL_RENAME[col] if col in COL_RENAME else col for col in df.columns]

        # Get alleles
        df['GT'].fillna('./.', inplace=True)
        df = df.loc[df['GT'].apply(analib.variant.gt_to_ac) > 0]

        # Filter on FILTER column
        df.loc[df['FILTER'] != 'PASS'].to_csv(output.tsv_filt, sep='\t', index=False, compression='gzip')

        df = df.loc[df['FILTER'] == 'PASS']

        # Alter fields
        df['SVTYPE'] = df['ALT'].apply(lambda val: re.sub(r'<([A-Z]+)(:.*)?>', '\\1', val))

        df['SVLEN'] = np.abs(df['SVLEN'])

        df['END'] = df['POS'] + df['SVLEN']

        # Save VCF index
        df['VCF_IDX'] = df.index

        # Set index and arrange columns
        df['ID'] = analib.variant.get_variant_id(df)

        df = analib.variant.order_variant_columns(df, tail_cols=('ALT', 'QUAL', 'FILTER', 'VCF_IDX'))

        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Drop duplicates
        df.drop_duplicates('ID', inplace=True)

        # Write all variants
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# variant_melt_vcf_to_tab
#
# VCF to TSV.
rule variant_melt_vcf_to_tab:
    input:
        vcf=lambda wildcards: variant_global_sample_table_entry('melt', wildcards=wildcards)['DATA']
    output:
        tsv=temp('temp/variant/caller/melt-{seq_set}/bed/{sample}/tsv/variants_sv.tsv.gz')
    shell:
        """bcftools query -H -f"%CHROM\\t%POS\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/ASSESS\\t%INFO/TSD\\t%INFO/INTERNAL\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t%INFO/MEINFO\\t%INFO/DIFF\\t%INFO/LP\\t%INFO/RP\\t%INFO/PRIOR\\t%INFO/SR[\\t%GT][\\t%GL]\n" {input.vcf} | """
        """gzip """
        """> {output.tsv}"""
