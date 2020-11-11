"""
Parse SMRT-SV variants to BED.
"""

###################
### Definitions ###
###################

def _variant_smrtsv_bed_get_sample_vcf(wildcards):
    """
    Locate variant VCF in SMRT-SV output.
    """

    # Get directory and SMRT-SV version from the sample table
    sample_entry = variant_global_sample_table_entry('smrtsv', wildcards.sample)

    sample_dir = sample_entry['DATA']
    version = sample_entry['VERSION']

    # Get VCF File
    if version == "2":
        vcf_file = os.path.join(sample_dir, 'variants.vcf.gz')

    if os.path.isfile(vcf_file):
        return vcf_file

    elif version == "1":
        vcf_file = os.path.join(sample_dir, 'variants.vcf')

    if not os.path.isfile(vcf_file):
        raise RuntimeError('Cannot locate input VCF for sample "{}" in "{}"'.format(wildcards.sample, sample_dir))

    return vcf_file

def _variant_bed_fix_indel_del_records(subdf, ref_file_name):
    """
    Fix indel DEL records from SMRT-SV (SEQ has extra sequence padded onto it).

    :subdf: A subset of a dataframe to correct.
    """

    subdf = subdf.copy()

    # Open and update DEL records
    with pysam.FastaFile(ref_file_name) as ref_file:
        for index, row in subdf.iterrows():
            if row['SVTYPE'] == 'DEL':
                subdf.loc[index, 'SEQ'] = ref_file.fetch(row['#CHROM'], row['POS'], row['POS'] + row['SVLEN'])

    # Return
    return subdf


#############
### Rules ###
#############

#
# Indel and SV (make final BED and FASTA)
#

# variant_smrtsv_bed_fa
#
# Make final BED and FASTA of ins/del/inv sequences.
rule variant_smrtsv_bed_fa:
    input:
        bed='temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed'
    output:
        bed='results/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        fa='results/variant/caller/smrtsv/fasta/{sample}/all/{filter}/{vartype}_{svtype}.fa.gz',
        fai='results/variant/caller/smrtsv/fasta/{sample}/all/{filter}/{vartype}_{svtype}.fa.gz.fai'
    wildcard_constraints:
        svtype='ins|del|inv|dup'
    run:
        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        if 'SEQ' not in df or 'ID' not in df:
            raise RuntimeError('Variant BED is missing SEQ or ID field.')

        # Write FASTA
        with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
            SeqIO.write(analib.seq.bed_to_seqrecord_iter(df), out_file, 'fasta')

        shell("""samtools faidx {output.fa}""")

        # Remove SEQ column
        del(df['SEQ'])

        # Write
        df.to_csv(output.bed, sep='\t', index=False)

#
# Indel
#

# variant_smrtsv_bed_filter_indel
#
# Remove indels inside deletions and inversions
rule variant_smrtsv_bed_filter_indel:
    input:
        sv_del='results/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/sv_del.bed',
        sv_inv='results/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/sv_inv.bed',
        indel_ins='temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/filtered/indel_ins.bed',
        indel_del='temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/filtered/indel_del.bed'
    output:
        indel_ins=temp('temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/indel_ins.bed'),
        indel_del=temp('temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/indel_del.bed'),
        indel_rm='results/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/removed/indel_in_sv.bed',
    shell:
        """head -n 1 {input.indel_ins} > {output.indel_ins}; """
        """head -n 1 {input.indel_del} > {output.indel_del}; """
        """head -n 1 {input.indel_ins} > {output.indel_rm}; """
        """bedtools intersect -wa -v -sorted -a {input.indel_ins} -b {input.sv_del} {input.sv_inv} >> {output.indel_ins}; """
        """bedtools intersect -wa -v -sorted -a {input.indel_del} -b {input.sv_del} {input.sv_inv} >> {output.indel_del}; """
        """{{"""
        """    bedtools intersect -wa -u -sorted -a {input.indel_ins} -b {input.sv_del} {input.sv_inv}; """
        """    bedtools intersect -wa -u -sorted -a {input.indel_del} -b {input.sv_del} {input.sv_inv}; """
        """}} | sort -k1,1 -k2,2n >> {output.indel_rm}"""

# variant_smrtsv_bed_filter_region_indel
#
# Apply a BED filter to SVs.
rule variant_smrtsv_bed_filter_region_indel:
    input:
        indel_ins='temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/indel_ins.bed',
        indel_del='temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/indel_del.bed',
        filter=lambda wildcards: analib.variant.get_filter_bed(wildcards.filter, UCSC_REF_NAME, config, SVPOP_DIR)
    output:
        indel_ins=temp('temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/filtered/indel_ins.bed'),
        indel_del=temp('temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/filtered/indel_del.bed'),
        indel_filt='results/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/removed/indel_in_filter.bed'
    wildcard_constraints:
        filter='((?!all)[^\/]*)|all.+'
    run:

        if wildcards.filter != 'all':
            shell(
                """head -n 1 {input.indel_ins} > {output.indel_ins}; """
                """head -n 1 {input.indel_del} > {output.indel_del}; """
                """head -n 1 {input.indel_ins} > {output.indel_filt}; """
                """bedtools intersect -wa -v -sorted -a {input.indel_ins} -b {input.filter} -header >> {output.indel_ins}; """  # Filter insertions
                """bedtools intersect -wa -v -sorted -a {input.indel_del} -b {input.filter} -header >> {output.indel_del}; """  # Filter deletions
                """{{"""
                """    bedtools intersect -wa -u -sorted -a {input.indel_ins} -b {input.filter}; """  # Capture filtered insertions
                """    bedtools intersect -wa -u -sorted -a {input.indel_del} -b {input.filter}; """  # Capture filtered deletions
                """}} | """
                """sort -k1,1 -k2,2n >> {output.indel_filt}"""
            )
        else:
            shell(
                """cp {input.indel_ins} {output.indel_ins}; """
                """cp {input.indel_del} {output.indel_del}; """
                """touch {output.indel_filt}"""
            )

# variant_smrtsv_bed_correct_indel_calls
#
# SMRT-SV makes indel DEL calls with flanking sequence, and this makes its way into the VCF. Correct here.
rule variant_smrtsv_bed_correct_indel_calls:
    input:
        bed='temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/indel_all_uncorrected.bed',
        ref=config['reference']
    output:
        indel_ins=temp('temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/indel_ins.bed'),
        indel_del=temp('temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/indel_del.bed')
    params:
        cores=4
    run:

        version = variant_global_sample_table_entry('smrtsv', wildcards.sample)['VERSION']

        # Get cores
        cores = np.int32(params.cores)

        if cores < 1:
            raise RuntimeError('Cores parameters must be at least 1: {}'.format(cores))

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Do fix in parallel (fixes for version 1 and 2)
        if version in {"1", "2"}:
            df = analib.pd.apply_parallel(
                df, _variant_bed_fix_indel_del_records, cores,
                kwds={
                    'ref_file_name': input.ref
                }
            )

        # Write
        df.loc[df['SVTYPE'] == 'INS'].to_csv(output.indel_ins, sep='\t', index=False)
        df.loc[df['SVTYPE'] == 'DEL'].to_csv(output.indel_del, sep='\t', index=False)


#
# SV
#

# variant_smrtsv_bed_filter_sv
#
# Remove SVs inside inversions.
rule variant_smrtsv_bed_filter_sv:
    input:
        sv_ins='temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/filtered/sv_ins.bed',
        sv_del='temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/filtered/sv_del.bed',
        sv_inv='temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/sv_inv.bed',
    output:
        sv_ins=temp('temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/sv_ins.bed'),
        sv_del=temp('temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/sv_del.bed'),
        sv_rm='results/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/removed/sv_in_inv.bed'
    shell:
        """head -n 1 {input.sv_ins} > {output.sv_ins}; """
        """head -n 1 {input.sv_del} > {output.sv_del}; """
        """head -n 1 {input.sv_ins} > {output.sv_rm}; """
        """bedtools intersect -wa -v -sorted -a {input.sv_ins} -b {input.sv_inv} >> {output.sv_ins}; """
        """bedtools intersect -wa -v -sorted -a {input.sv_del} -b {input.sv_inv} >> {output.sv_del}; """
        """{{"""
        """    bedtools intersect -wa -u -sorted -a {input.sv_ins} -b {input.sv_inv}; """
        """    bedtools intersect -wa -u -sorted -a {input.sv_del} -b {input.sv_inv}; """
        """}} | """
        """sort -k1,1 -k2,2n >> {output.sv_rm}"""

# variant_smrtsv_bed_filter_region_sv
#
# Apply a BED filter to SVs.
rule variant_smrtsv_bed_filter_region_sv:
    input:
        sv_ins='temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/sv_ins.bed',
        sv_del='temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/sv_del.bed',
        sv_inv='temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/sv_inv.bed',
        filter=lambda wildcards: analib.variant.get_filter_bed(wildcards.filter, UCSC_REF_NAME, config, SVPOP_DIR)
    output:
        sv_ins=temp('temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/filtered/sv_ins.bed'),
        sv_del=temp('temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/filtered/sv_del.bed'),
        sv_inv=temp('temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/sv_inv.bed'),
        sv_filt='results/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/removed/sv_in_filter.bed'
    wildcard_constraints:
        filter='((?!all)[^\/]*)|all.+'
    run:

        if wildcards.filter != 'all':
            shell(
                """head -n 1 {input.sv_ins} > {output.sv_ins}; """
                """head -n 1 {input.sv_del} > {output.sv_del}; """
                """head -n 1 {input.sv_inv} > {output.sv_inv}; """
                """head -n 1 {input.sv_ins} > {output.sv_filt}; """
                """bedtools intersect -wa -v -sorted -a {input.sv_ins} -b {input.filter} -header > {output.sv_ins}; """  # Insertions
                """bedtools intersect -wa -v -sorted -a {input.sv_del} -b {input.filter} -header > {output.sv_del}; """  # Deletions
                """bedtools intersect -wa -v -sorted -a {input.sv_inv} -b {input.filter} -header > {output.sv_inv}; """  # Inversions
                """{{"""
                """    bedtools intersect -wa -u -sorted -a {input.sv_ins} -b {input.filter}; """  # Capture filtered insertions
                """    bedtools intersect -wa -u -sorted -a {input.sv_del} -b {input.filter}; """  # Capture filtered deletions
                """    bedtools intersect -wa -u -sorted -a {input.sv_inv} -b {input.filter}; """  # Capture filtered inversions
                """}} | """
                """sort -k1,1 -k2,2n -o {output.sv_filt} """  # Sort filtered SVs
            )
        else:
            shell(
                """cp {input.sv_ins} {output.sv_ins}; """
                """cp {input.sv_del} {output.sv_del}; """
                """cp {input.sv_inv} {output.sv_inv}; """
                """touch {output.sv_filt}; """
            )

# variant_smrtsv_bed_tab_to_bed
#
# Format: bed3+12
# 1) Chromosome
# 2) Start position (0-based, inclusive)
# 3) End position (0-based, exclusive)
# 4) ID (CHROM-POS-SVTYPE-SVLEN)
# 5) Type (INS, DEL, INV)
# 6) Number of bases affected by the event
# 7) Quality score
# 8) Number of supporting contigs
# 9) Number of contigs at this locus
# 10) Contig variant was modeled on
# 11) Location in contig where variant starts (probably BED-style)
# 12) Location in contig where variant ends (probably BED-style)
# 13) Repeat annotation
# 14) Reference base
# 15) Sequence affected by the variant
rule variant_smrtsv_bed_tab_to_bed:
    input:
        ref=config['reference'],
        tab='temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/sv_all.tab'
    output:
        indel=temp('temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/indel_all_uncorrected.bed'),
        sv_ins=temp('temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/sv_ins.bed'),
        sv_del=temp('temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/sv_del.bed'),
        sv_inv=temp('temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/sv_inv.bed'),
        sv_no_support='results/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/sv_without_contig_support.bed'
    run:

        # Read variants
        var_df = pd.read_csv(input.tab, sep='\t', header=0)

        # Get abs version of SVLEN
        var_df['SVLEN_ABS'] = np.abs(var_df['SVLEN'])

        # Add variant ID
        var_df['ID'] = analib.variant.get_variant_id(var_df)

        # Fix end position
        var_df['END'] = var_df.apply(
            lambda var_record:
                var_record['POS'] + var_record['SVLEN_ABS'] if var_record['SVTYPE'] != 'INS' else var_record['POS'] + 1,
                axis=1
        )

        # Sort (should already be sorted from SMRTSV output)
        var_df.sort_values(['#CHROM', 'POS', 'END'], inplace=True)

        # Rearrange columns
        var_df = var_df.loc[:,
            ('#CHROM', 'POS', 'END', 'ID', 'FILTER',
                'SVTYPE', 'SVLEN',
                'QUAL', 'CONTIG_SUPPORT', 'CONTIG_DEPTH',
                'CONTIG', 'CONTIG_START', 'CONTIG_END',
                'REPEAT_TYPE', 'REF', 'SEQ'
            )
        ]

        # Filter negative length variants (Rare SMRTSV bug)
        #var_df = var_df[var_df.apply(lambda var_row: var_row['END'] >= var_row['POS'], axis=1)]

        del(var_df['FILTER'])

        # Separate variants
        df_indel = var_df.loc[abs(var_df['SVLEN']) < 50, :]
        df_sv_supported = var_df.loc[[a and b for a, b in zip(abs(var_df['SVLEN']) >= 50, var_df['CONTIG_SUPPORT'] > 1)], :]
        df_sv_no_support = var_df.loc[[a and b for a, b in zip(abs(var_df['SVLEN']) >= 50, var_df['CONTIG_SUPPORT'] == 1)], :]

        # Get
        df_ins = df_sv_supported.loc[df_sv_supported['SVTYPE'] == 'INS', :]
        df_del = df_sv_supported.loc[df_sv_supported['SVTYPE'] == 'DEL', :]
        df_inv = df_sv_supported.loc[df_sv_supported['SVTYPE'] == 'INV', :].copy()

        # Add SEQ for INV calls
        with pysam.FastaFile(input.ref) as in_file:
            df_inv['SEQ'] = df_inv.apply(lambda row: in_file.fetch(row['#CHROM'], row['POS'], row['END']), axis=1)

        # Write
        df_indel.to_csv(output.indel, sep='\t', index=False, header=True, na_rep='.', float_format='%.0f')

        df_ins.to_csv(output.sv_ins, sep='\t', index=False, header=True, na_rep='.', float_format='%.0f')
        df_del.to_csv(output.sv_del, sep='\t', index=False, header=True, na_rep='.', float_format='%.0f')
        df_inv.to_csv(output.sv_inv, sep='\t', index=False, header=True, na_rep='.', float_format='%.0f')

        df_sv_no_support.to_csv(output.sv_no_support, sep='\t', index=False, header=True, na_rep='.', float_format='%.0f')

# variant_smrtsv_bed_vcf_to_tab
#
# Convert variant VCF to a tab-separated file with columns for annotations.
rule variant_smrtsv_bed_vcf_to_tab:
    input:
        vcf=_variant_smrtsv_bed_get_sample_vcf
    output:
        tab=temp('temp/variant/caller/smrtsv/bed/{sample}/all/all/byref/noqc/sv_all.tab')
    shell:
        """vcfkeepinfo {input.vcf} END SVTYPE SVLEN CONTIG_SUPPORT CONTIG_DEPTH CONTIG CONTIG_START CONTIG_END REPEAT_TYPE SEQ | """
        """vcf2tsv -n 'NA' """
        """> {output.tab}"""
