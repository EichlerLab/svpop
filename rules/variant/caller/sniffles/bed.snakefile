"""
Process Sniffles variant calls into BED files.
"""

#
# Indel
#

# variant_sniffles_bed_indel
#
# Make emtpy indel BED (temporary rule).
rule variant_sniffles_bed_indel:
    output:
        bed_ins='results/variant/caller/sniffles/bed/{sample}/all/{filter}/byref/indel_ins.bed',
        bed_del='results/variant/caller/sniffles/bed/{sample}/all/{filter}/byref/indel_del.bed',
        fa_ins='results/variant/caller/sniffles/fasta/{sample}/all/{filter}/indel_ins.fa.gz',
        fa_del='results/variant/caller/sniffles/fasta/{sample}/all/{filter}/indel_del.fa.gz'
    shell:
        """echo -e "#CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN" > {output.bed_ins}; """
        """echo -e "#CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN" > {output.bed_del}; """
        """touch {output.fa_ins} {output.fa_del}"""

#
# SV
#

# variant_sniffles_bed_make_fa
#
# Write FASTA of SV sequence.
rule variant_sniffles_bed_make_fa:
    input:
        bed='results/variant/caller/sniffles/bed/{sample}/all/{filter}/byref/sv_{svtype}.bed',
        ref=config['reference']
    output:
        fa='results/variant/caller/sniffles/fasta/{sample}/all/{filter}/sv_{svtype}.fa.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup'
    run:

        if wildcards.svtype != 'ins':
            # Read
            df = pd.read_csv(input.bed, sep='\t', header=0)

            # Get sequences
            with pysam.FastaFile(input.ref) as in_file:
                df['SEQ'] = df.apply(lambda row: in_file.fetch(row['#CHROM'], row['POS'], row['END']), axis=1)

            # Write FASTA
            with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
                SeqIO.write(analib.seq.bed_to_seqrecord_iter(df), out_file, 'fasta')

            shell("""samtools faidx {output.fa}""")

        else:
            shell("""touch {output.fa}""")


# variant_sniffles_bed_filter_region_sv
#
# Apply a BED filter to SVs.
rule variant_sniffles_bed_filter_region_sv:
    input:
        bed='temp/variant/caller/sniffles/bed/{sample}/all/all/byref/sv_{svtype}.bed',
        filter=lambda wildcards: analib.variant.get_filter_bed(wildcards.filter, UCSC_REF_NAME, config, SVPOP_DIR)
    output:
        bed='results/variant/caller/sniffles/bed/{sample}/all/{filter}/byref/sv_{svtype}.bed',
        sv_filt='results/variant/caller/sniffles/bed/{sample}/all/{filter}/byref/removed/sv_{svtype}_in_filter.bed'
    wildcard_constraints:
        filter='((?!all)[^\/]*)|all.+'
    run:

        if wildcards.filter != 'all':
            shell(
                """bedtools intersect -wa -v -sorted -a {input.bed} -b {input.filter} -header > {output.bed}; """
                """bedtools intersect -wa -u -sorted -a {input.bed} -b {input.filter} -header > {output.sv_filt}"""
            )
        else:
            shell(
                """cp {input.bed} {output.bed}; """
                """touch {output.sv_filt}; """
            )

# variant_sniffles_bed_sv_bed
#
# Make SV BED.
rule variant_sniffles_bed_sv_bed:
    input:
        tab='temp/variant/caller/sniffles/bed/{sample}/all/all/byref/noqc/sv_all.tab'
    output:
        bed_anno='results/variant/caller/sniffles/bed/{sample}/all/all/byref/info/caller_annotations.bed',
        sv_ins=temp('temp/variant/caller/sniffles/bed/{sample}/all/all/byref/sv_ins.bed'),
        sv_del=temp('temp/variant/caller/sniffles/bed/{sample}/all/all/byref/sv_del.bed'),
        sv_inv=temp('temp/variant/caller/sniffles/bed/{sample}/all/all/byref/sv_inv.bed'),
        sv_dup=temp('temp/variant/caller/sniffles/bed/{sample}/all/all/byref/sv_dup.bed'),
        bpt_bnd='results/variant/caller/sniffles/bed/{sample}/all/all/byref/bpt_bnd.bed'
    run:

        # Read
        df = pd.read_csv(input.tab, sep='\t', header=0)

        # Remove uninformative columns
        del(df['REF'], df['ALT'], df['QUAL'])

        # Normalize fields
        df['SVLEN'] = np.abs(df['SVLEN'])

        df['SVTYPE_ORG'] = df['SVTYPE']
        df['SVTYPE'] = df['SVTYPE'].apply(lambda val: val.replace('/', ''))

        df['END'] = df.apply(lambda row: row['POS'] + (row['SVLEN'] if row['SVTYPE'] != 'INS' else 1), axis=1)

        # Set ID and arrange columns
        df['ID'] = analib.variant.get_variant_id(df)

        df = analib.variant.order_variant_columns(df)

        # Sort
        # Sniffles sometimes writes duplicate entries for variant calls. Sorting by GT (descending so 1/1 comes first)
        # DV (descending so higher values of supporting reads come first) and DR (ascending so lower values reference
        # reads come first) allows duplicates to be removed so that keeping the first variant in the table gives the
        # strongest version of the variant.
        df.sort_values(
            ['#CHROM', 'POS', 'SVTYPE', 'SVLEN', 'GT', 'DV', 'DR'],
            ascending=(True, True, True, True, False, False, True),
            inplace=True
        )

        # Write full caller annotations
        df.to_csv(output.bed_anno, sep='\t', index=False)

        # Filter by VCF column
        df = df.loc[df['FILTER'] == 'PASS']

        # Drop duplicates (Sniffles sometimes writes duplicate variant calls).
        df.drop_duplicates(['ID'], inplace=True)

        # Remove annotation columns (recorded in the annotation table already written)
        for rm_field in (
            'PRECISE', 'IMPRECISE',
            'Kurtosis_quant_start', 'Kurtosis_quant_stop', 'STD_quant_start', 'STD_quant_stop',
            'MAPQ', 'RE', 'REF_strand', 'STRANDS', 'SUPTYPE', 'AF', 'SVMETHOD', 'ZMW', 'SAMPLE',
            'FILTER', 'SVTYPE_ORG', 'UNRESOLVED'
        ):
            del(df[rm_field])


        # Rename fields
        prepend_caller_set = {'GT', 'DV', 'DR', 'CHR2'}

        df.columns = [('SNIF_' + col if col in prepend_caller_set else col) for col in df.columns]

        # Get SV tables
        df_ins = df.loc[df['SVTYPE'] == 'INS'].copy()
        df_del = df.loc[df['SVTYPE'] == 'DEL'].copy()
        df_inv = df.loc[df['SVTYPE'] == 'INV'].copy()
        df_dup = df.loc[(df['SVTYPE'] == 'DUP') | (df['SVTYPE'] == 'INVDUP')].copy()
        df_bnd = df.loc[df['SVTYPE'] == 'BND'].copy()

        # Annotate inverted duplications
        df_dup['SNIF_INVDUP'] = df_dup['SVTYPE'] == 'INVDUP'
        df_dup['SVTYPE'] = 'DUP'

        # CHR2 only has information for BND and INS
        del(df_del['SNIF_CHR2'])
        del(df_inv['SNIF_CHR2'])
        del(df_dup['SNIF_CHR2'])

        # Write SV
        df_ins.loc[df_ins['SVLEN'] >= 50].to_csv(output.sv_ins, sep='\t', index=False)
        df_del.loc[df_del['SVLEN'] >= 50].to_csv(output.sv_del, sep='\t', index=False)
        df_inv.loc[df_inv['SVLEN'] >= 50].to_csv(output.sv_inv, sep='\t', index=False)
        df_dup.loc[df_dup['SVLEN'] >= 50].to_csv(output.sv_dup, sep='\t', index=False)
        df_bnd.to_csv(output.bpt_bnd, sep='\t', index=False)


#
# VCF to TAB
#

# variant_sniffles_bed_vcf_to_tab
#
# Convert variant VCF to a tab-separated file with columns for annotations.
rule variant_sniffles_bed_vcf_to_tab:
    input:
        vcf=lambda wildcards: variant_global_sample_table_entry('sniffles', wildcards.sample)['DATA']
    output:
        tab=temp('temp/variant/caller/sniffles/bed/{sample}/all/all/byref/noqc/sv_all.tab')
    shell:
        """vcf2tsv -n 'NA' -g {input.vcf} """
        """> {output.tab}"""
