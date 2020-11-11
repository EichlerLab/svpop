"""
1000 Genomes Phase 1 variants.

McVean, G. A., Altshuler, D. M., Durbin, R. M., Abecasis, G. R., Bentley, D. R., Chakravarti, A., … McVean, G. A. (2012). An integrated map of genetic variation from 1,092 human genomes. Nature, 491(7422), 56–65. http://doi.org/10.1038/nature11632
"""

# varset_set_1kgp1_annotate
#
# Annotate variants.
rule varset_set_1kgp1_annotate:
    input:
        bed='temp/varset/1kgp1/byref/sv_all.bed'
    output:
        bed_ins='temp/variant/varset/1kgp1/bed/all/all/all/byref/sv_ins.bed',
        bed_del='temp/variant/varset/1kgp1/bed/all/all/all/byref/sv_del.bed',
        bed_inv='temp/variant/varset/1kgp1/bed/all/all/all/byref/sv_inv.bed',
        bed_dup='temp/variant/varset/1kgp1/bed/all/all/all/byref/sv_dup.bed',
        bed_cnv='temp/variant/varset/1kgp1/bed/all/all/all/byref/sv_cnv.bed'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Separate BED by type and write
        for svtype in ('INS', 'DEL', 'INV', 'DUP', 'CNV'):
            df.loc[df['SVTYPE'] == svtype].to_csv(output['bed_{}'.format(svtype.lower())], sep='\t', index=False)


# varset_set_1kgp1_liftover
#
# Liftover variants to hg38.
rule varset_set_1kgp1_liftover:
    input:
        bed='temp/varset/set/1kgp1/all/hg19/sv_all.bed',
        chain='data/liftover/hg19_to_hg38.chain.gz',
    output:
        bed=temp('temp/varset/1kgp1/byref/sv_all.bed'),
        unmapped='results/varset/set/1kgp1/all/all/byref/unmapped.bed'
    shell:
        """{{\n"""
        """    head -n 1 {input.bed}; """
        """    liftOver -bedPlus=4 {input.bed} {input.chain} /dev/stdout {output.unmapped} | sort -k1,1 -k 2,2n; """
        """}} """
        """> {output.bed}"""

# varset_set_1kgp1_tab_to_bed
#
# Variant TAB to BED.
rule varset_set_1kgp1_tab_to_bed:
    input:
        tab='temp/varset/set/1kgp1/variants.tab.gz'
    output:
        snv='results/varset/set/1kgp1/all/all/byref/rm/snv.tab.gz',
        indel='results/varset/set/1kgp1/all/all/byref/rm/indel_all.tab.gz',
        sv_short='results/varset/set/1kgp1/all/all/byref/rm/sv_short.tab.gz',
        bed_sv=temp('temp/varset/set/1kgp1/all/hg19/sv_all.bed')
    run:

        # Read variants
        df = pd.read_csv(input.tab, sep='\t', header=0, dtype={'#CHROM': np.object})

        # GRCh37 to hg37 chromosome names
        df['#CHROM'] = analib.ref.grc_to_hg_chrom(df['#CHROM'], 'GRCh37')

        # Write and filter SNPs and Indels
        df.loc[df['VT'] == 'SNP'].to_csv(output.snv, sep='\t', index=False, compression='gzip')
        df.loc[df['VT'] == 'INDEL'].to_csv(output.indel, sep='\t', index=False, compression='gzip')

        df = df.loc[df['VT'].apply(lambda val: val not in {'SNP', 'INDEL'})]

        # Fix missing lengths (for uncertain sizes)
        df['POS'] = np.abs(df['POS'].astype(np.int16))
        df['SVLEN'] = df.apply(lambda row: np.int16(row['SVLEN']) if row['SVLEN'] != '.' else np.int16(row['END']) - row['POS'], axis=1).astype(np.int16)
        df['SVLEN'] = np.abs(df['SVLEN'])

        # Filter short SVs
        df.loc[df['SVLEN'] < 50].to_csv(output.sv_short, sep='\t', index=False, compression='gzip')

        df = df.loc[df['SVLEN'] >= 50]

        # Set ID
        df['ACC'] = df['ID']
        df['ID'] = analib.variant.get_variant_id(df)

        # Write
        df.loc[
            :, ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'ACC', 'AFR_AF', 'EUR_AF', 'QUAL', 'FILTER')
        ].to_csv(output.bed_sv, sep='\t', index=False)

# varset_set_1kgp1_sv_table
#
# VCF to variant table.
rule varset_set_1kgp1_sv_table:
    input:
        vcf='temp/varset/set/1kgp1/variants.vcf.gz'
    output:
        tab=temp('temp/varset/set/1kgp1/variants.tab.gz')
    shell:
        """zcat {input.vcf} | """
        """vcf2tsv | """
        """gzip > {output.tab}"""

# varset_set_1kgp1_dl_variants
#
# Get VCF.
rule varset_set_1kgp1_dl_variants:
    output:
        vcf=temp('temp/varset/set/1kgp1/variants.vcf.gz')
    params:
        url=config['varset']['set']['1kgp1']['url']
    shell:
        """wget {params.url} -O {output.vcf}"""
