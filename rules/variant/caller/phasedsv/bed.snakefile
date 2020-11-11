"""
Extract calls from Phased-SV
"""

# variant_phasedsv_bed_empty
#
# Create an empty DUP and INV bed.
rule variant_phasedsv_bed_empty:
    output:
        bed='results/variant/caller/phasedsv/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        fa='results/variant/caller/phasedsv/fasta/{sample}/all/{filter}/{vartype}_{svtype}.fa.gz',
        fai='results/variant/caller/phasedsv/fasta/{sample}/all/{filter}/{vartype}_{svtype}.fa.gz.fai'
    wildcard_constraints:
        svtype='inv|dup'
    shell:
        """echo -e "#CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN" > {output.bed}; """
        """touch {output.fa} {output.fai}"""

# variant_phasedsv_bed_fa
#
# Make final BED and FASTA of ins/del/inv sequences.
rule variant_phasedsv_bed_fa:
    input:
        bed='temp/variant/caller/phasedsv/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed'
    output:
        bed='results/variant/caller/phasedsv/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        fa='results/variant/caller/phasedsv/fasta/{sample}/all/{filter}/{vartype}_{svtype}.fa.gz',
        fai='results/variant/caller/phasedsv/fasta/{sample}/all/{filter}/{vartype}_{svtype}.fa.gz.fai'
    wildcard_constraints:
        svtype='ins|del'
    run:
        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        if df.shape[0] > 0:
            if 'SEQ' not in df or 'ID' not in df:
                raise RuntimeError('Variant BED is missing SEQ or ID field.')

            # Write FASTA
            with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
                SeqIO.write(analib.seq.bed_to_seqrecord_iter(df), out_file, 'fasta')

            shell("""samtools faidx {output.fa}""")
        else:
            shell("""touch {output.fa} {output.fai}""")

        # Remove SEQ column
        del(df['SEQ'])

        # Write
        df.to_csv(output.bed, sep='\t', index=False)

# variant_phasedsv_bed_filter_region_sv
#
# Apply a BED filter to SVs.
rule variant_phasedsv_bed_filter_region_sv:
    input:
        sv_ins='temp/variant/caller/phasedsv/bed/{sample}/all/all/byref/{vartype}_ins.bed',
        sv_del='temp/variant/caller/phasedsv/bed/{sample}/all/all/byref/{vartype}_del.bed',
        filter=lambda wildcards: analib.variant.get_filter_bed(wildcards.filter, UCSC_REF_NAME, config, SVPOP_DIR)
    output:
        sv_ins=temp('temp/variant/caller/phasedsv/bed/{sample}/all/{filter}/byref/{vartype}_ins.bed'),
        sv_del=temp('temp/variant/caller/phasedsv/bed/{sample}/all/{filter}/byref/{vartype}_del.bed'),
        sv_filt='results/variant/caller/phasedsv/bed/{sample}/all/{filter}/byref/removed/{vartype}_in_filter.bed'
    wildcard_constraints:
        filter='((?!all)[^\/]*)|all.+'
    run:

        if wildcards.filter != 'all':
            shell(
                """bedtools intersect -wa -v -sorted -a {input.sv_ins} -b {input.filter} -header > {output.sv_ins}; """  # Insertions
                """bedtools intersect -wa -v -sorted -a {input.sv_del} -b {input.filter} -header > {output.sv_del}; """  # Deletions
                """{{"""
                """    bedtools intersect -wa -u -sorted -a {input.sv_ins} -b {input.filter}; """  # Capture filtered insertions
                """    bedtools intersect -wa -u -sorted -a {input.sv_del} -b {input.filter}; """  # Capture filtered deletions
                """}} | """
                """sort -k1,1 -k2,2n -o {output.sv_filt} """  # Sort filtered SVs
            )
        else:
            shell(
                """cp {input.sv_ins} {output.sv_ins}; """
                """cp {input.sv_del} {output.sv_del}; """
                """touch {output.sv_filt}; """
            )

# variant_phasedsv_bed_tab_to_bed
#
# Table to variant BED.
rule variant_phasedsv_bed_tab_to_bed:
    input:
        ref=config['reference'],
        tab='temp/variant/caller/phasedsv/bed/{sample}/all/all/byref/{vartype}_all.tab'
    output:
        sv_ins=temp('temp/variant/caller/phasedsv/bed/{sample}/all/all/byref/{vartype}_ins.bed'),
        sv_del=temp('temp/variant/caller/phasedsv/bed/{sample}/all/all/byref/{vartype}_del.bed')
    run:

        # Read variants
        var_df = pd.read_csv(input.tab, sep='\t', header=0)

        del(var_df['REF'], var_df['ALT'], var_df['QUAL'], var_df['FILTER'], var_df['SAMPLE'])

        # Adjust SVTYPE (is "insertion" and "deletion" in indels)
        var_df['SVTYPE'] = var_df['SVTYPE'].apply(lambda val: val[:3].upper())

        # ABS SVLEN
        var_df['SVLEN'] = np.abs(var_df['SVLEN'])

        # Add variant ID
        if var_df.shape[0] > 0:
            var_df['ID'] = analib.variant.get_variant_id(var_df)

            # Fix end position
            var_df['END'] = var_df.apply(
                lambda var_record:
                    var_record['POS'] + var_record['SVLEN'] if var_record['SVTYPE'] != 'INS' else var_record['POS'] + 1,
                    axis=1
            )

            # Sort (should already be sorted from phasedsv output)
            var_df.sort_values(['#CHROM', 'POS', 'END'], inplace=True)

        # Drop duplicates
        var_df.drop_duplicates('ID', inplace=True, keep='first')

        # Rearrange columns
        var_df = analib.variant.order_variant_columns(
            var_df,
            tail_cols=('GT', 'CONTIG', 'CONTIG_START', 'CONTIG_END', 'SEQ'),
            allow_missing=True
        )

        # Write
        var_df.iloc[
            list(var_df['SVTYPE'] == 'INS')
        ].to_csv(
            output.sv_ins, sep='\t', index=False, header=True, na_rep='.', float_format='%.0f'
        )

        var_df.iloc[
            list(var_df['SVTYPE'] == 'DEL')
        ].to_csv(
            output.sv_del, sep='\t', index=False, header=True, na_rep='.', float_format='%.0f'
        )


# variant_phasedsv_bed_sv_vcf_to_tab
#
# Convert SV VCF to a tab-separated file with columns for annotations.
rule variant_phasedsv_bed_indel_vcf_to_tab:
    input:
        vcf=lambda wildcards: variant_global_sample_table_entry('phasedsv', wildcards.sample)['DATA'].split(';')[1]
    output:
        tab=temp('temp/variant/caller/phasedsv/bed/{sample}/all/all/byref/indel_all.tab')
    shell:
        """vcfkeepinfo {input.vcf} END SVTYPE SVLEN SEQ | """
        """vcf2tsv -g -n 'NA' """
        """> {output.tab}"""

# variant_phasedsv_bed_sv_vcf_to_tab
#
# Convert SV VCF to a tab-separated file with columns for annotations.
rule variant_phasedsv_bed_sv_vcf_to_tab:
    input:
        vcf=lambda wildcards: variant_global_sample_table_entry('phasedsv', wildcards.sample)['DATA'].split(';')[0]
    output:
        tab=temp('temp/variant/caller/phasedsv/bed/{sample}/all/all/byref/sv_all.tab')
    shell:
        """vcfkeepinfo {input.vcf} END SVTYPE SVLEN CONTIG CONTIG_START CONTIG_END SEQ | """
        """vcf2tsv -g -n 'NA' """
        """> {output.tab}"""
