"""
Genomes of the Netherlands (GoNL) variants.

Francioli, L. C., Menelaou, A., Pulit, S. L., Van Dijk, F., Palamara, P. F., Elbers, C. C., … Wijmenga, C. (2014). Whole-genome sequence variation, population structure and demographic history of the Dutch population. Nature Genetics, 46(8), 818–825. https://doi.org/10.1038/ng.3021

Francioli, L. C., Polak, P. P., Koren, A., Menelaou, A., Chun, S., Renkens, I., … Sunyaev, S. R. (2015). Genome-wide patterns and properties of de novo mutations in humans. Nature Genetics, 47(7), 822–826. https://doi.org/10.1038/ng.3292
"""

rule varset_set_gonl_make_bed:
    input:
        bed='temp/varset/gonl/byref/sv_all.bed'
    output:
        bed_ins='results/varset/set/gonl/all/all/byref/sv_ins.bed',
        bed_del='results/varset/set/gonl/all/all/byref/sv_del.bed',
        bed_inv='results/varset/set/gonl/all/all/byref/sv_inv.bed',
        bed_dup='results/varset/set/gonl/all/all/byref/sv_dup.bed',
        bed_cnv='results/varset/set/gonl/all/all/byref/sv_cnv.bed'
    run:
        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Separate BED by type and write
        for svtype in ('INS', 'DEL', 'INV', 'DUP', 'CNV'):
            df.loc[df['SVTYPE'] == svtype].to_csv(output['bed_{}'.format(svtype.lower())], sep='\t', index=False)

# varset_set_gonl_liftover
#
# Liftover variants to hg38.
rule varset_set_gonl_liftover:
    input:
        bed='temp/varset/set/gonl/all/hg37/sv_all.bed',
        chain='data/liftover/hg19_to_hg38.chain.gz'
    output:
        bed=temp('temp/varset/gonl/byref/sv_all.bed'),
        unmapped='results/varset/set/gonl/all/all/byref/rm/liftover_unmapped.bed'
    shell:
        """{{\n"""
        """    head -n 1 {input.bed}; """
        """    liftOver -bedPlus=4 {input.bed} {input.chain} /dev/stdout {output.unmapped} | sort -k1,1 -k 2,2n; """
        """}} """
        """> {output.bed}"""

# varset_set_gonl_tab_to_bed
#
# Parse variants into BED files.
rule varset_set_gonl_tab_to_bed:
    input:
        tab='temp/varset/set/gonl/variants.tab.gz'
    output:
        bed_rm_nosv='results/varset/set/gonl/all/all/byref/rm/sv_nosv_hg37.bed',
        bed_sv=temp('temp/varset/set/gonl/all/hg37/sv_all.bed')
    run:

        # Read variants
        df = pd.read_csv(input.tab, sep='\t', header=0, usecols=('#CHROM', 'POS', 'ID', 'ALT', 'SVTYPE', 'SVLEN', 'DISCOVERY', 'FILTER', 'QUAL'))

        # GRCh37 to hg37 chromosome names
        df['#CHROM'] = analib.ref.grc_to_hg_chrom(df['#CHROM'], 'GRCh37')

        # Remove variants without an SVLEN (removes 32 MEI insertions in GoNL 6.1)
        df_nosv = df.loc[df['SVLEN'] == '.']
        df = df.loc[df['SVLEN'] != '.']

        # Get length
        df['SVLEN'] = df['SVLEN'].apply(np.int16)
        df['SVLEN'] = df['SVLEN'].apply(np.abs)

        # Filter indel and 14 translocations; Write all non-SV calls
        df_nosv = pd.concat([df_nosv, df.loc[df['SVLEN'] < 50]], axis=0)

        df_nosv.to_csv(output.bed_rm_nosv, sep='\t', index=False)
        del(df_nosv)

        df = df.loc[df['SVLEN'] >= 50]

        # Set columns
        df['RS'] = df['ID']
        df['ID'] = analib.variant.get_variant_id(df)

        df['END'] = df['POS'] + df['SVLEN']

        # Write
        df = df.loc[:, ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'RS', 'DISCOVERY', 'FILTER', 'QUAL')]

        df.to_csv(output.bed_sv, sep='\t', index=False)

# varset_set_gonl_sv_table
#
# VCF to variant table.
rule varset_set_gonl_sv_table:
    input:
        vcf='temp/varset/set/gonl/variants.vcf.gz'
    output:
        tab=temp('temp/varset/set/gonl/variants.tab.gz')
    shell:
        """zcat {input.vcf} | """
        """vcf2tsv | """
        """gzip > {output.tab}"""


# varset_set_gonl_vcf_dl
#
# Get GoNL VCF.
rule varset_set_gonl_vcf_dl:
    output:
        vcf=temp('temp/varset/set/gonl/variants.vcf.gz')
    params:
        url=config['varset']['set']['gonl']['vcf']
    shell:
        """wget {params.url} -O {output.vcf}"""
