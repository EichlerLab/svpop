"""
Parse PBSV variants to BED.
"""

###################
### Definitions ###
###################

# Discard variants with a proportion of N's in SEQ of this value or higher.
N_PROP_CUTOFF = 0.1

# Minimum duplication size
MIN_DUP = 1000

# Minimum inversion size
MIN_INV = 1000

# bcftools query string
VARIANT_BED_PBSV_BCFTOOLS_QUERY = '"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO/SVTYPE\t%INFO/END\t%INFO/SVLEN\t%INFO/CIPOS\t%INFO/MATEID\t%INFO/MATEDIST\t%INFO/IMPRECISE[\t%GT][\t%AD][\t%DP][\t%SAC][]\n"'

def variant_bed_pbsv_tsv_inv_dup(wildcards):
    """
    Get input to tsv-to-bed rule. If all variants were output to a single VCF, find all in "variants_sv.tsv.gz",
    otherwise, look in "variants_{vartype}.tsv.gz".

    :param wildcards: Rule wildcards.

    :return: Path to input TSV to convert to BED.
    """

    sample_entry = svpoplib.rules.sample_table_entry('pbsv', SAMPLE_TABLE, wildcards=wildcards)

    if wildcards.varsvtype not in {'sv_inv', 'dup_dup'}:
        raise RuntimeError('PBSV Parser input function received varsvtype = {varsvtype}: Expected "sv_inv" or "dup_dup"'.format(**wildcards))

    if wildcards.varsvtype == 'sv_inv' or 'vartype' not in sample_entry['WILDCARDS']:
        return 'temp/variant/caller/pbsv-{seq_set}/bed/{sample}/tsv/variants_sv.tsv.gz'.format(**wildcards)

    else:
        return 'temp/variant/caller/pbsv-{seq_set}/bed/{sample}/tsv/variants_dup.tsv.gz'.format(**wildcards)

#############
### Rules ###
#############

# variant_pbsv_bed_tab_to_bed_dupinv
#
# Parse duplication and inversion variants to a BED file.
rule variant_pbsv_bed_tab_to_bed_dupinv:
    input:
        tsv=variant_bed_pbsv_tsv_inv_dup
    output:
        bed=temp('temp/variant/caller/pbsv-{seq_set}/{sample}/all/all/bed/pre_filter/{varsvtype}.bed.gz'),
        fa=temp('temp/variant/caller/pbsv-{seq_set}/{sample}/all/all/bed/pre_filter/fa/{varsvtype}.fa.gz'),
        filtered='results/variant/caller/pbsv-{seq_set}/{sample}/all/all/bed/filtered/filtered_{varsvtype}.bed.gz'
    wildcard_constraints:
        varsvtype='sv_inv|dup_dup'
    run:

        vartype, svtype = (val.upper() for val in wildcards.varsvtype.split('_'))

        # Read
        df = pd.read_csv(input.tsv, sep='\t', header=0)

        df = svpoplib.varbed.bcftools_query_to_tsv(df, wildcards.sample)

        df = df.loc[df['SVTYPE'] == svtype]

        df['END'] = df['END'].astype(int)

        df['SVLEN'] = df['END'] - df['POS']

        del(
            df['CIPOS'], df['IMPRECISE'], df['MATEDIST'], df['MATEID'],
            df['REF'], df['ALT']
        )

        # Set ID
        df['PBSV_ID'] = df['ID']
        df['ID'] = svpoplib.variant.get_variant_id(df)

        # Arrange columns
        df = svpoplib.variant.order_variant_columns(df)

        # Get SEQ
        # df['SEQ'] = svpoplib.ref.get_ref_region(df, config['reference'])

        # Fail set: Proportion of N's in SEQ is too high
        # fail_set_n = set(df.loc[df.apply(
        #     lambda row: (
        #             np.sum([val == 'N' for val in row['SEQ'].upper()]) / row['SVLEN']
        #         ) >= N_PROP_CUTOFF if not pd.isnull(row['SEQ']) else False,
        #     axis=1
        # )].index)

        # Fail set: Not PASS in FILTER
        df['FAIL_REASON'] = np.nan
        fail_set = set(df.loc[df['FILTER'] != 'PASS'].index)

        df.loc[fail_set, 'FAIL_REASON'] = 'FILTER != PASS'

        # Remove filtered variants
        df.loc[
            fail_set
        ].sort_values(
            ['#CHROM', 'POS']
        ).to_csv(
            output.filtered, sep='\t', index=False, compression='gzip'
        )

        df = df.loc[[val not in fail_set for val in df.index]].copy()

        del(df['FAIL_REASON'])

        # Count and drop duplicates (ins/del)
        id_count = df.groupby('ID').apply(lambda subdf: subdf.shape[0])

        df['PBSV_ID_COUNT'] = df['ID'].apply(lambda id: id_count[id]).copy()

        df.drop_duplicates(subset='ID', keep='first', inplace=True)

        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Normalize columns
        no_rename = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'SEQ', 'PBSV_ID', 'PBSV_ID_COUNT']

        df.columns = [('PBSV_' + col if col not in no_rename else col) for col in df.columns]

        # Empty FASTA (INV and DUP are not sequence-resolved in PBSV)
        with open(output.fa, 'w') as out_file:
            pass

        # Write
        df.to_csv(output.bed, sep='\t', na_rep='NA', index=False)

# variant_pbsv_bed_tsv_to_bed_sv
#
# Parse SV/indel variants to a BED file.
rule variant_pbsv_bed_tsv_to_bed_sv:
    input:
        tsv='temp/variant/caller/pbsv-{seq_set}/bed/{sample}/tsv/variants_sv.tsv.gz'
    output:
        indel_ins=temp('temp/variant/caller/pbsv-{seq_set}/{sample}/all/all/bed/pre_filter/indel_ins.bed.gz'),
        indel_del=temp('temp/variant/caller/pbsv-{seq_set}/{sample}/all/all/bed/pre_filter/indel_del.bed.gz'),
        sv_ins=temp('temp/variant/caller/pbsv-{seq_set}/{sample}/all/all/bed/pre_filter/sv_ins.bed.gz'),
        sv_del=temp('temp/variant/caller/pbsv-{seq_set}/{sample}/all/all/bed/pre_filter/sv_del.bed.gz'),
        fa_indel_ins=temp('temp/variant/caller/pbsv-{seq_set}/{sample}/all/all/bed/pre_filter/fa/indel_ins.fa.gz'),
        fa_indel_del=temp('temp/variant/caller/pbsv-{seq_set}/{sample}/all/all/bed/pre_filter/fa/indel_del.fa.gz'),
        fa_sv_ins=temp('temp/variant/caller/pbsv-{seq_set}/{sample}/all/all/bed/pre_filter/fa/sv_ins.fa.gz'),
        fa_sv_del=temp('temp/variant/caller/pbsv-{seq_set}/{sample}/all/all/bed/pre_filter/fa/sv_del.fa.gz'),
        filtered='results/variant/caller/pbsv-{seq_set}/{sample}/all/all/bed/filtered/filtered_sv_insdel.bed.gz'
    run:

        # Read and format/subset to this sample
        df = pd.read_csv(input.tsv, sep='\t')

        df = svpoplib.varbed.bcftools_query_to_tsv(df, wildcards.sample)

        # Subset on SVTYPE
        df = df.loc[df['SVTYPE'].apply(lambda val: val in {'INS', 'DEL'})]

        del(
            df['CIPOS'],
            df['MATEID'],
            df['MATEDIST']
        )

        # Get SV POS, END, SVLEN, SVTYPE, and SEQ
        df_coord = df.apply(svpoplib.variant.vcf_fields_to_seq, axis=1)

        df['POS'] = df_coord['POS']
        df['END'] = df_coord['END']
        df['SVLEN'] = df_coord['SVLEN']
        df['SVTYPE'] = df_coord['SVTYPE']
        df['SEQ'] = df_coord['SEQ']
        df['FAIL_REASON'] = np.nan

        # Set ID
        df['PBSV_ID'] = df['ID']
        df['ID'] = svpoplib.variant.get_variant_id(df)

        # Arrange columns
        df = svpoplib.variant.order_variant_columns(df)

        # Fail set: Proportion of N's in SEQ is too high
        fail_set_nseq = set(df.loc[df.apply(
            lambda row: (
                    np.sum([val == 'N' for val in row['SEQ'].upper()]) / row['SVLEN']
                ) >= N_PROP_CUTOFF if not pd.isnull(row['SEQ']) else False,
            axis=1
        )].index)

        df.loc[fail_set_nseq, 'FAIL_REASON'] = 'SEQ N proportion'

        # Fail set: INS in N
        fail_set_ins_n = set(df.loc[(df['SVTYPE'] == 'INS') & (df['REF'].apply(lambda val: val.upper() == 'N'))].index)

        df.loc[fail_set_ins_n, 'FAIL_REASON'] = 'INS in N'

        # Fail set: Not PASS in FILTER
        fail_set_filter = set(df.loc[df['FILTER'] != 'PASS'].index)

        df.loc[fail_set_filter, 'FAIL_REASON'] = 'FILTER != PASS'

        fail_set = fail_set_nseq | fail_set_ins_n | fail_set_filter

        # Remove filtered variants
        df.loc[
            fail_set
        ].sort_values(
            ['#CHROM', 'POS']
        ).to_csv(
            output.filtered, sep='\t', index=False, compression='gzip'
        )

        df = df.loc[[val not in fail_set for val in df.index]].copy()

        del(df['FAIL_REASON'])

        # Count and drop duplicates (ins/del)
        id_count = df.groupby('ID').apply(lambda subdf: subdf.shape[0])

        df['PBSV_ID_COUNT'] = df['ID'].apply(lambda id: id_count[id])

        df.drop_duplicates(subset='ID', keep='first', inplace=True)

        # Normalize columns
        del(df['REF'], df['ALT'], df['FILTER'])

        no_rename = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'SEQ', 'PBSV_ID', 'PBSV_ID_COUNT']

        df.columns = [('PBSV_' + col if col not in no_rename else col) for col in df.columns]
        
        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Write FASTA files
        for vartype in ('sv', 'indel'):

            if vartype == 'sv':
                id_set_vartype = set(df.loc[df['SVLEN'] >= 50, 'ID'])
            else:
                id_set_vartype = set(df.loc[df['SVLEN'] < 50, 'ID'])

            for svtype in ('ins', 'del'):

                id_set_svtype = set(df.loc[df['SVTYPE'] == svtype.upper(), 'ID'])

                # Write FASTA
                out_file_name = output[f'fa_{vartype}_{svtype}']

                with Bio.bgzf.BgzfWriter(out_file_name, 'wb') as out_file:
                    SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(df, id_set_vartype & id_set_svtype), out_file, 'fasta')

                shell("""samtools faidx {out_file_name}""")

        # Remove SEQ column
        del(df['SEQ'])

        # Write indels
        df_sub = df.loc[df['SVLEN'] < 50]
        df_sub.loc[df_sub['SVTYPE'] == 'INS'].to_csv(output.indel_ins, sep='\t', na_rep='NA', index=False, compression='gzip')
        df_sub.loc[df_sub['SVTYPE'] == 'DEL'].to_csv(output.indel_del, sep='\t', na_rep='NA', index=False, compression='gzip')

        # Write SVs
        df_sub = df.loc[df['SVLEN'] >= 50]
        df_sub.loc[df_sub['SVTYPE'] == 'INS'].to_csv(output.sv_ins, sep='\t', na_rep='NA', index=False, compression='gzip')
        df_sub.loc[df_sub['SVTYPE'] == 'DEL'].to_csv(output.sv_del, sep='\t', na_rep='NA', index=False, compression='gzip')


# variant_dv_vcf_to_tsv
#
# VCF to TSV file.
rule variant_pbsv_bed_vcf_to_tsv:
    input:
        vcf=lambda wildcards: svpoplib.rules.sample_table_entry('pbsv', SAMPLE_TABLE, wildcards=wildcards)['DATA']
    output:
        tsv=temp('temp/variant/caller/pbsv-{seq_set}/bed/{sample}/tsv/variants_{vartype}.tsv.gz')
    params:
        query_string=VARIANT_BED_PBSV_BCFTOOLS_QUERY  # Could adjust by PBSV version if needed
    shell:
        """bcftools query -H -f{params.query_string} {input.vcf} | """
        """gzip """
        """> {output.tsv}"""

