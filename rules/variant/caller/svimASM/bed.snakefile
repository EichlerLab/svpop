"""
Process svimASM variant calls.
"""

# variant_svimASM_bed
#
# Write filtered BED and empty FASTA files.
rule variant_svimASM_bed:
    input:
        bed='temp/variant/caller/svimASM-{seq_set}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed'
    output:
        bed='results/variant/caller/svimASM-{seq_set}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        fa='results/variant/caller/svimASM-{seq_set}/fasta/{sample}/all/{filter}/{vartype}_{svtype}.fa.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup'
    run:

        df = pd.read_csv(input.bed, sep='\t')

        df.to_csv(output.bed, sep='\t', index=False)

        with open(output.fa, 'w'):
            pass

# variant_svimASM_bed_filter_region
#
# Apply a BED filter to SVs.
rule variant_svimASM_bed_filter_region:
    input:
        bed='temp/variant/caller/svimASM-{seq_set}/bed/{sample}/all/all/{vartype}_{svtype}.bed',
        filter=lambda wildcards: analib.variant.get_filter_bed(wildcards.filter, UCSC_REF_NAME, config, SVPOP_DIR)
    output:
        bed=temp('temp/variant/caller/svimASM-{seq_set}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed'),
        bed_filt='results/variant/caller/svimASM-{seq_set}/bed/{sample}/all/{filter}/byref/removed/{vartype}_{svtype}_filtered.bed'
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

# variant_dv_tsv_to_bed
#
# TSV to BED.
rule variant_svimASM_tsv_to_bed:
    input:
        tsv='temp/variant/caller/svimASM-{seq_set}/bed/{sample}/tsv/variants.tsv.gz'
    output:
        sv_ins=temp('temp/variant/caller/svimASM-{seq_set}/bed/{sample}/all/all/sv_ins.bed'),
        sv_del=temp('temp/variant/caller/svimASM-{seq_set}/bed/{sample}/all/all/sv_del.bed'),
        sv_inv=temp('temp/variant/caller/svimASM-{seq_set}/bed/{sample}/all/all/sv_inv.bed'),
        indel_del=temp('temp/variant/caller/svimASM-{seq_set}/bed/{sample}/all/all/indel_del.bed'),
        indel_ins=temp('temp/variant/caller/svimASM-{seq_set}/bed/{sample}/all/all/indel_ins.bed'),
        tsv_filt='results/variant/caller/svimASM-{seq_set}/bed/{sample}/all/all/filtered/vcf_fail.tsv.gz'
    run:

        # Definitions
        COL_RENAME = {
            'CHROM': '#CHROM'
        }

        REQUIRED_COLS = {'CHROM', 'POS', 'QUAL', 'FILTER', 'GT'}

        # Read
        df = pd.read_csv(input.tsv, sep='\t', header=0, low_memory=False)

        # Remove prefixes added by bcftools (e.g. "# [1]CHROM" to "CHROM", "[2]POS" to "POS"
        df.columns = [re.sub('^#?\s*\[[^\]]+\]', '', col) for col in df.columns]

        # Save VCF index
        df['VCF_IDX'] = df.index

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

        # Drop records without at least one ALT allele
        df = df.loc[df['GT'].apply(
            lambda val: np.any([False if val in {'0', '.'} else True for val in re.split('/|\|', val)])
        )]

        # Check columns and rename
        missing_cols = REQUIRED_COLS - set(df.columns)

        if missing_cols:
            raise RuntimeError('Missing VCF columns: {}'.format(', '.join(sorted(missing_cols))))

        df.columns = [COL_RENAME[col] if col in COL_RENAME else col for col in df.columns]

        # Filter on FILTER column
        df.loc[df['FILTER'] != 'PASS'].to_csv(output.tsv_filt, sep='\t', index=False, compression='gzip')

        df = df.loc[df['FILTER'] == 'PASS']
        df = df.loc[(df['SVTYPE'] == 'INS') | (df['SVTYPE'] == 'INV') | (df['SVTYPE'] == 'DEL')]

        # Set VCF fields
        for index in df.index:
            try:
                df.at[index, 'SVLEN'] = np.abs(int(df.at[index,'SVLEN']))
            except:
                df.at[index,'SVLEN'] = np.abs((int(df.at[index, 'END']) - int(df.at[index,'POS'])))

        # df['END'] = df['POS'] + df.apply(lambda row: 1 if row['SVTYPE'] == 'INS' else abs(row['SVLEN']), axis=1)

        df['ID'] = analib.variant.get_variant_id(df)

        del(df['FILTER'])
        del(df['REF'])
        del(df['ALT'])

        df = analib.variant.order_variant_columns(df)

        # Drop duplicates
        df.drop_duplicates('ID', inplace=True)

        # Write all variants
        df.sort_values(['#CHROM', 'POS'], inplace=True)

        df_sub = df.loc[df['SVLEN'] < 50]
        df_sub.loc[df_sub['SVTYPE'] == 'INS'].to_csv(output.indel_ins, sep='\t', na_rep='NA', index=False)
        df_sub.loc[df_sub['SVTYPE'] == 'DEL'].to_csv(output.indel_del, sep='\t', na_rep='NA', index=False)

        # Write SVs
        df_sub = df.loc[df['SVLEN'] >= 50]
        df_sub.loc[df_sub['SVTYPE'] == 'INS'].to_csv(output.sv_ins, sep='\t', na_rep='NA', index=False)
        df_sub.loc[df_sub['SVTYPE'] == 'DEL'].to_csv(output.sv_del, sep='\t', na_rep='NA', index=False)
        df_sub.loc[df_sub['SVTYPE'] == 'INV'].to_csv(output.sv_inv, sep='\t', na_rep='NA', index=False)

# variant_bn_vcf_to_tsv
#
# VCF to TSV.
rule variant_svimASM_vcf_to_tsv:
    input:
        vcf=lambda wildcards: variant_global_sample_table_entry('svimASM', wildcards=wildcards)['DATA']
    output:
        tsv=temp('temp/variant/caller/svimASM-{seq_set}/bed/{sample}/tsv/variants.tsv.gz')
    shell:
        """
        bcftools query -H -f"%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/SVTYPE\\t%INFO/END\\t%INFO/SVLEN[\\t%GT\\t%CN]\n" {input.vcf} | """
        """gzip """
        """> {output.tsv}"""



