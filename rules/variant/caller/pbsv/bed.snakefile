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


#############
### Rules ###
#############

#
# Make final BED and FASTA
#

# variant_pbsv_bed_fa
#
# Make final BED and FASTA of ins/del/inv sequences.
rule variant_pbsv_bed_fa:
    input:
        bed='temp/variant/caller/pbsv-{seq_set}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed'
    output:
        bed='results/variant/caller/pbsv-{seq_set}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        fa='results/variant/caller/pbsv-{seq_set}/fasta/{sample}/all/{filter}/{vartype}_{svtype}.fa.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup'
    run:
        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

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


#
# Variant calls
#

# variant_pbsv_bed_filter_region
#
# Apply a BED filter to SVs.
rule variant_pbsv_bed_filter_region:
    input:
        bed='temp/variant/caller/pbsv-{seq_set}/bed/{sample}/all/all/nofilter/{vartype}_{svtype}.bed',
        filter=lambda wildcards: analib.variant.get_filter_bed(wildcards.filter, UCSC_REF_NAME, config, SVPOP_DIR)
    output:
        bed=temp('temp/variant/caller/pbsv-{seq_set}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed'),
        bed_filt='results/variant/caller/pbsv-{seq_set}/bed/{sample}/all/{filter}/byref/removed/{vartype}_{svtype}_filtered.bed'
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

# variant_pbsv_bed_tab_to_bed_dupinv
#
# Parse duplication and inversion variants to a BED file.
rule variant_pbsv_bed_tab_to_bed_dupinv:
    input:
        tab=lambda wildcards:
            'temp/variant/caller/pbsv-{seq_set}/bed/{sample}/tab/variants_dup.tab.gz'.format(**wildcards) if wildcards.svtype == 'dup' else
            'temp/variant/caller/pbsv-{seq_set}/bed/{sample}/tab/variants_sv.tab.gz'.format(**wildcards)
    output:
        bed=temp('temp/variant/caller/pbsv-{seq_set}/bed/{sample}/all/all/nofilter/sv_{svtype}.bed'),
        filtered='results/variant/caller/pbsv-{seq_set}/bed/{sample}/all/all/byref/filtered/filtered_in_vcf_{svtype}.bed.gz'
    wildcard_constraints:
        svtype='inv|dup'
    run:

        # Read
        df = pd.read_csv(input.tab, sep='\t', header=0)

        df = df.loc[df['SVTYPE'] == wildcards.svtype.upper()]

        df['SVLEN'] = df['END'] - df['POS']

        del(
            df['QUAL'], df['CIPOS'], df['IMPRECISE'], df['MATEDIST'], df['MATEID'], df['SHADOWED'], df['SVANN'],
            df['REF'], df['ALT']
        )

        # Set ID
        df['PBSV_ID'] = df['ID']
        df['ID'] = analib.variant.get_variant_id(df)

        # Arrange columns
        df = analib.variant.order_variant_columns(df)

        # Get SEQ
        df['SEQ'] = analib.ref.get_ref_region(df, config['reference'])

        # Fail set: Proportion of N's in SEQ is too high
        fail_set = set(df.loc[df.apply(
            lambda row: (
                    np.sum([val == 'N' for val in row['SEQ'].upper()]) / row['SVLEN']
                ) >= N_PROP_CUTOFF if not pd.isnull(row['SEQ']) else False,
            axis=1
        )].index)

        # Fail set: Not PASS in FILTER
        fail_set |= set(df.loc[df['FILTER'] != 'PASS'].index)

        # Remove filtered variants
        df.loc[
            fail_set
        ].sort_values(
            ['#CHROM', 'POS']
        ).to_csv(
            output.filtered, sep='\t', index=False, compression='gzip'
        )

        df = df.loc[[val not in fail_set for val in df.index]].copy()

        # Count and drop duplicates (ins/del)
        id_count = df.groupby('ID').apply(lambda subdf: subdf.shape[0])

        df['PBSV_CALL_COUNT'] = df['ID'].apply(lambda id: id_count[id]).copy()

        df.drop_duplicates(subset='ID', keep='first', inplace=True)

        # Normalize columns
        no_rename = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'SEQ', 'PBSV_ID', 'PBSV_CALL_COUNT']

        df.columns = [('PBSV_' + col if col not in no_rename else col) for col in df.columns]

        # Write
        df.to_csv(output.bed, sep='\t', na_rep='NA', index=False)

# variant_pbsv_bed_tab_to_bed_sv
#
# Parse SV/indel variants to a BED file.
rule variant_pbsv_bed_tab_to_bed_sv:
    input:
        tab='temp/variant/caller/pbsv-{seq_set}/bed/{sample}/tab/variants_sv.tab.gz'
    output:
        indel_ins=temp('temp/variant/caller/pbsv-{seq_set}/bed/{sample}/all/all/nofilter/indel_ins.bed'),
        indel_del=temp('temp/variant/caller/pbsv-{seq_set}/bed/{sample}/all/all/nofilter/indel_del.bed'),
        sv_ins=temp('temp/variant/caller/pbsv-{seq_set}/bed/{sample}/all/all/nofilter/sv_ins.bed'),
        sv_del=temp('temp/variant/caller/pbsv-{seq_set}/bed/{sample}/all/all/nofilter/sv_del.bed'),
        filtered='results/variant/caller/pbsv-{seq_set}/bed/{sample}/all/all/byref/filtered/filtered_in_vcf_sv.bed.gz'
    run:

        # Read
        df = pd.read_csv(input.tab, sep='\t', header=0)

        df = df.loc[df['SVTYPE'].apply(lambda val: val in {'INS', 'DEL'})]

        del(
            df['CIPOS'], df['IMPRECISE'], df['MATEDIST'], df['MATEID'], df['SHADOWED'], df['SVANN']
        )

        # Get SV POS, END, SVLEN, SVTYPE, and SEQ
        df_coord = df.apply(analib.variant.vcf_fields_to_seq, axis=1)

        df['POS'] = df_coord['POS']
        df['END'] = df_coord['END']
        df['SVLEN'] = df_coord['SVLEN']
        df['SVTYPE'] = df_coord['SVTYPE']
        df['SEQ'] = df_coord['SEQ']

        # Set ID
        df['PBSV_ID'] = df['ID']
        df['ID'] = analib.variant.get_variant_id(df)

        # Arrange columns
        df = analib.variant.order_variant_columns(df)

        # Fail set: Proportion of N's in SEQ is too high
        fail_set = set(df.loc[df.apply(
            lambda row: (
                    np.sum([val == 'N' for val in row['SEQ'].upper()]) / row['SVLEN']
                ) >= N_PROP_CUTOFF if not pd.isnull(row['SEQ']) else False,
            axis=1
        )].index)

        # Fail set: INS in N
        fail_set |= set(df.loc[(df['SVTYPE'] == 'INS') & (df['REF'].apply(lambda val: val.upper() == 'N'))].index)

        # Fail set: Not PASS in FILTER
        fail_set |= set(df.loc[df['FILTER'] != 'PASS'].index)

        # Remove filtered variants
        df.loc[
            fail_set
        ].sort_values(
            ['#CHROM', 'POS']
        ).to_csv(
            output.filtered, sep='\t', index=False, compression='gzip'
        )

        df = df.loc[[val not in fail_set for val in df.index]].copy()

        # Count and drop duplicates (ins/del)
        id_count = df.groupby('ID').apply(lambda subdf: subdf.shape[0])

        df['PBSV_CALL_COUNT'] = df['ID'].apply(lambda id: id_count[id])

        df.drop_duplicates(subset='ID', keep='first', inplace=True)

        # Normalize columns
        del(df['AC'], df['REF'], df['ALT'], df['FILTER'], df['AN'], df['NS'])

        no_rename = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'SEQ', 'PBSV_ID', 'PBSV_CALL_COUNT']

        df.columns = [('PBSV_' + col if col not in no_rename else col) for col in df.columns]

        # Write indels
        df_sub = df.loc[df['SVLEN'] < 50]
        df_sub.loc[df_sub['SVTYPE'] == 'INS'].to_csv(output.indel_ins, sep='\t', na_rep='NA', index=False)
        df_sub.loc[df_sub['SVTYPE'] == 'DEL'].to_csv(output.indel_del, sep='\t', na_rep='NA', index=False)

        # Write SVs
        df_sub = df.loc[df['SVLEN'] >= 50]
        df_sub.loc[df_sub['SVTYPE'] == 'INS'].to_csv(output.sv_ins, sep='\t', na_rep='NA', index=False)
        df_sub.loc[df_sub['SVTYPE'] == 'DEL'].to_csv(output.sv_del, sep='\t', na_rep='NA', index=False)


# variant_pbsv_bed_vcf_to_tab
#
# VCF to tab.
rule variant_pbsv_bed_vcf_to_tab:
    input:
        vcf=lambda wildcards: variant_global_sample_table_entry('pbsv', wildcards=wildcards)['DATA']
    output:
        tab=temp('temp/variant/caller/pbsv-{seq_set}/bed/{sample}/tab/variants_{vartype}.tab.gz')
    shell:
        """vcfkeepsamples {input.vcf} {wildcards.sample} | """
        """vcffixup - | """
        """vcffilter -f "AC > 0" | """
        """vcf2tsv -n 'NA' | """
        """gzip """
        """> {output.tab}"""
