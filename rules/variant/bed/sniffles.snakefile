"""
Parse PBSV variants to BED.
"""

###################
### Definitions ###
###################

# Minimum duplication size
MIN_DUP = 1000

# Minimum inversion size
MIN_INV = 1000

# bcftools query string
VARIANT_BED_SNIFFLES_BCFTOOLS_QUERY = '"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/END\t%INFO/RE\t%INFO/SUPTYPE\t%INFO/ZMW[\t%DV]\n"'



#############
### Rules ###
#############

# variant_bed_sniffles_tsv_to_bed
#
# Parse SV/indel variants to a BED file.
rule variant_bed_sniffles_tsv_to_bed:
    input:
        tsv='temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/bed/{sample}/tsv/variants.tsv.gz'
    output:
        sv_ins=temp('temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/{sample}/all/all/bed/pre_filter/sv_ins.bed.gz'),
        sv_del=temp('temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/{sample}/all/all/bed/pre_filter/sv_del.bed.gz'),
        sv_inv=temp('temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/{sample}/all/all/bed/pre_filter/sv_inv.bed.gz'),
        sv_dup=temp('temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/{sample}/all/all/bed/pre_filter/sv_dup.bed.gz'),
        indel_ins=temp('temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/{sample}/all/all/bed/pre_filter/indel_ins.bed.gz'),
        indel_del=temp('temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/{sample}/all/all/bed/pre_filter/indel_del.bed.gz'),
        fa_sv_ins=temp('temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/{sample}/all/all/bed/pre_filter/fa/sv_ins.fa.gz'),
        fa_sv_del=temp('temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/{sample}/all/all/bed/pre_filter/fa/sv_del.fa.gz'),
        fa_sv_inv=temp('temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/{sample}/all/all/bed/pre_filter/fa/sv_inv.fa.gz'),
        fa_sv_dup=temp('temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/{sample}/all/all/bed/pre_filter/fa/sv_dup.fa.gz'),
        fa_indel_ins=temp('temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/{sample}/all/all/bed/pre_filter/fa/indel_ins.fa.gz'),
        fa_indel_del=temp('temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/{sample}/all/all/bed/pre_filter/fa/indel_del.fa.gz'),
        filtered='results/variant/caller/sniffles/{sourcename_base}-{seq_set}/{sample}/all/all/bed/filtered/filtered_sv_insdel.bed.gz'
    wildcard_constraints:
        vartype='sv|indel'
    run:

        # Read and format/subset to this sample
        df = pd.read_csv(input.tsv, sep='\t')

        df = svpoplib.varbed.bcftools_query_to_tsv(df, wildcards.sample, filter_gt=False)

        # Assign complex type
        df['COMPLEX'] = df['SVTYPE'].apply(lambda val:
            val if val.upper() not in {'INS', 'DEL', 'INV', 'DUP', 'BND', 'SNV'} else np.nan
        )

        df['SVTYPE'] = df.apply(lambda row:
            row['SVTYPE'] if pd.isnull(row['COMPLEX']) else (
                'DEL' if 'DEL' in row['SVTYPE'].upper() else (
                    'DUP' if 'DUP' in row['SVTYPE'].upper() else (
                        'INV' if 'INV' in row['SVTYPE'].upper() else (
                            'INS' if 'INS' in row['SVTYPE'].upper() else (
                                'BND' if row['SVTYPE'].upper() == 'BND' else \
                                    np.nan
                            )
                        )
                    )
                )
            ), axis=1
        )

        # Move ID
        df['SNIF_ID'] = df['ID']
        del(df['ID'])

        # Positive SVLEN
        df['SVLEN'] = np.abs(df['SVLEN'])

        # Filter
        filter_pass = df['FILTER'] != 'PASS'
        filter_any = filter_pass

        filter_min_inv = df.apply(lambda row: (row['SVTYPE'] == 'INV') & (row['SVLEN'] < MIN_INV), axis=1) & ~ filter_any

        filter_any |= filter_min_inv

        filter_min_dup = df.apply(lambda row: (row['SVTYPE'] == 'INV') & (row['SVLEN'] < MIN_INV), axis=1) & ~ filter_any

        filter_any |= filter_min_dup

        df['FAIL_REASON'] = np.nan

        df_filter = df.loc[filter_any]

        df_filter.loc[filter_pass, 'FAIL_REASON'] = 'FILTER != PASS'
        df_filter.loc[filter_min_inv, 'FAIL_REASON'] = 'SVLEN < MIN_INV'
        df_filter.loc[filter_min_dup, 'FAIL_REASON'] = 'SVLEN < MIN_DUP'

        df_filter = df.loc[filter_any].to_csv(output.filtered, sep='\t', index=False, compression='gzip')

        df = df.loc[~ filter_any].copy()

        del(df['FAIL_REASON'])

        # Drop BND and fix END (BND END is "." changing the dtype from int to object)
        df = df.loc[df['SVTYPE'] != 'BND'].copy()
        df['END'] = df['END'].astype(int)

        # Set SEQ
        df['SEQ'] = df.apply(lambda row: row['REF'] if row['SVTYPE'] == 'DEL' else (row['ALT'] if row['SVTYPE'] == 'INS' else np.nan), axis=1)

        # Set ID
        df['ID'] = svpoplib.variant.get_variant_id(df)
        df = svpoplib.variant.order_variant_columns(df)

        # Count and drop duplicates (ins/del)
        id_count = df.groupby('ID').apply(lambda subdf: subdf.shape[0])

        df['SNIF_ID_COUNT'] = df['ID'].apply(lambda id: id_count[id])

        df.drop_duplicates(subset='ID', keep='first', inplace=True)

        # Normalize columns
        del(df['REF'], df['ALT'], df['FILTER'])

        no_rename = {'#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'SEQ'}

        df.columns = [('SNIF_' + col if (col not in no_rename and not col.startswith('SNIF_')) else col) for col in df.columns]

        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Write FASTA files
        for vartype in ('sv', 'indel'):

            if vartype == 'sv':
                id_set_vartype = set(df.loc[df['SVLEN'] >= 50, 'ID'])
            else:
                id_set_vartype = set(df.loc[df['SVLEN'] < 50, 'ID'])

            for svtype in ('ins', 'del', 'inv', 'dup'):

                if svtype in {'inv', 'dup'} and vartype == 'indel':
                    continue

                id_set_svtype = set(df.loc[df['SVTYPE'] == svtype.upper(), 'ID'])

                # Write FASTA
                out_file_name = output[f'fa_{vartype}_{svtype}']

                if svtype not in {'inv', 'dup'}:

                    with Bio.bgzf.BgzfWriter(out_file_name, 'wb') as out_file:
                        SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(df, id_set_vartype & id_set_svtype), out_file, 'fasta')

                    shell("""samtools faidx {out_file_name}""")
                else:
                    with open(out_file_name, 'w'):
                        pass

        # Remove SEQ column
        del(df['SEQ'])

        # Write indels
        df_sub = df.loc[df['SVLEN'] < 50]
        df_sub.loc[df_sub['SVTYPE'] == 'INS'].to_csv(output.indel_ins, sep='\t', na_rep='NA', index=False, compression='gzip')
        df_sub.loc[df_sub['SVTYPE'] == 'DEL'].to_csv(output.indel_del, sep='\t', na_rep='NA', index=False, compression='gzip')

        # Write SVs
        df_sub = df.loc[df['SVLEN'] >= 50]
        df_sub.loc[df_sub['SVTYPE'] == 'DEL'].to_csv(output.sv_del, sep='\t', na_rep='NA', index=False, compression='gzip')
        df_sub.loc[df_sub['SVTYPE'] == 'INS'].to_csv(output.sv_ins, sep='\t', na_rep='NA', index=False, compression='gzip')
        df_sub.loc[df_sub['SVTYPE'] == 'INV'].to_csv(output.sv_inv, sep='\t', na_rep='NA', index=False, compression='gzip')
        df_sub.loc[df_sub['SVTYPE'] == 'DUP'].to_csv(output.sv_dup, sep='\t', na_rep='NA', index=False, compression='gzip')


# variant_bed_sniffles_vcf_to_tsv
#
# VCF to TSV file.
rule variant_bed_sniffles_vcf_to_tsv:
    input:
        vcf=lambda wildcards: svpoplib.rules.sample_table_entry(wildcards.sourcename_base, SAMPLE_TABLE, wildcards=wildcards, type='sniffles')['DATA']
    output:
        tsv=temp('temp/variant/caller/sniffles/{sourcename_base}-{seq_set}/bed/{sample}/tsv/variants.tsv.gz')
    params:
        query_string=VARIANT_BED_SNIFFLES_BCFTOOLS_QUERY  # Could adjust by caller version if needed
    shell:
        """bcftools query -H -f{params.query_string} {input.vcf} | """
        """gzip """
        """> {output.tsv}"""

