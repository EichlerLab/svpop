"""
SVs published by Mills (2011) from 1000 Genomes pilots 1 and 2.

Mills, R. E., Walter, K., Stewart, C., Handsaker, R. E., Chen, K., Alkan, C., … Korbel, J. O. (2011). Mapping copy number variation by population-scale genome sequencing. Nature, 470(7332), 59–65. http://doi.org/10.1038/nature09708
"""

#
# Bylen and empty sets
#

# varset_set_mills2011_ins_inv
#
# Write empty files for others.
rule varset_set_mills2011_inv:
    output:
        bed='results/varset/set/mills2011/all/all/byref/sv_inv.bed'
    shell:
        """echo -e "#CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN\tINFO" >{output.bed}; """


#
# Merge and Liftover
#

# varset_set_mills2011_liftover
#
# Liftover variants to hg38.
rule varset_set_mills2011_liftover:
    input:
        bed='temp/varset/set/mills2011/hg18/sv_{svtype}.bed',
        chain='data/liftover/hg18_to_hg38.chain.gz',
    output:
        bed='results/varset/set/mills2011/all/all/byref/sv_{svtype,ins|del}.bed',
        temp=temp('temp/varset/mills2011/byref/sv_{svtype}.bed'),
        unmapped='results/varset/set/mills2011/all/all/byref/unmapped_{svtype}.bed'
    shell:
        """liftOver {input.bed} {input.chain} {output.temp} {output.unmapped}; """
        """head -n 1 {input.bed} >{output.bed}; """
        """sort -k1,1 -k 2,2n {output.temp} >>{output.bed}"""

# varset_set_mills2011_merge_ins
#
# Merge MEI and novel INS.
rule varset_set_mills2011_merge_ins:
    input:
        bed_mei='temp/varset/set/mills2011/hg18/sv_ins_mei.bed',
        bed_nov='temp/varset/set/mills2011/hg18/sv_ins_novel.bed'
    output:
        bed=temp('temp/varset/set/mills2011/hg18/sv_ins.bed')
    run:

        # Read
        df = pd.concat(
            [
                pd.read_csv(input.bed_mei, sep='\t', header=0),
                pd.read_csv(input.bed_nov, sep='\t', header=0)
            ],
            axis=0
        )

        # Set ID
        df['POS1'] = df['POS'] + 1
        df['ID'] = df.apply(lambda row: '{#CHROM}-{POS1}-{SVTYPE}-{SVLEN}'.format(**row), axis=1)
        del(df['POS1'])

        # Sort
        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Arrange columns
        df = df.loc[:, ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'INFO')]

        # Write
        df.to_csv(output.bed, sep='\t', index=False)

#
# VCF to BED
#

rule varset_set_mills20100_tab_to_bed:
    input:
        tab='temp/varset/set/mills2011/hg18/sv_del.tab'
    output:
        bed=temp('temp/varset/set/mills2011/hg18/sv_del.bed')
    run:

        # Input
        df = pd.read_csv(input.tab, sep='\t', header=0)

        # Set BED fields
        df['SVLEN'] = np.abs(df['SVLEN'])
        df['END'] = df['POS'] + df['SVLEN']
        df['SVTYPE'] = 'DEL'

        # Set ID
        df['POS1'] = df['POS'] + 1
        df['ID'] = df.apply(lambda row: '{#CHROM}-{POS1}-{SVTYPE}-{SVLEN}'.format(**row), axis=1)
        del(df['POS1'])

        # Arrange columns
        df = df.loc[:, ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'INFO')]

        # Write
        df.to_csv(output.bed, sep='\t', index=False)


# varset_set_mills2011_del_to_tab
#
# DEL to BED.
rule varset_set_mills2011_del_to_tab:
    input:
        vcf='temp/varset/set/mills2011/variants_del.vcf.gz'
    output:
        tab=temp('temp/varset/set/mills2011/hg18/sv_del.tab')
    shell:
        """echo -e "#CHROM\tPOS\tSVLEN\tINFO" >{output.tab}; """
        """java -jar ${{SNPSIFT_JAR}} extractFields -e "." {input.vcf} CHROM POS SVLEN SVMETHOD | """
        """awk -vOFS="\\t" '"""
            """(NR > 1 && $3 != "." && $3 <= -50) {{$3 = -$3; print "chr"$1, $2, $3, $4}}"""
        """' | """
        """sort -k1,1 -k 2,2n """
        """>>{output.tab}"""

# varset_set_mills2011_ins_to_bed
#
# MEI and novel VCF to BED.
rule varset_set_mills2011_ins_to_bed:
    input:
        vcf='temp/varset/set/mills2011/variants_{svsubtype}.vcf.gz'
    output:
        bed=temp('temp/varset/set/mills2011/hg18/sv_ins_{svsubtype,mei|novel}.bed')
    shell:
        """echo -e "#CHROM\tPOS\tEND\tSVTYPE\tSVLEN\tINFO" >{output.bed}; """
        """java -jar ${{SNPSIFT_JAR}} extractFields -e "." {input.vcf} CHROM POS SVLEN ID ALT | """
        """awk -vOFS="\\t" '"""
            """(NR > 1 && $3 != ".") {{$5 = gensub("<|>", "", "g", $5); print "chr"$1, $2 - 1, $2, "INS", $3, $4 ";" $5}}"""
        """' | """
        """sort -k1,1 -k 2,2n """
        """>>{output}"""


#
# Download
#

# varset_set_mills2011_dl_del
#
# Get DEL variants.
rule varset_set_mills2011_dl_del:
    output:
        vcf=temp('temp/varset/set/mills2011/variants_del.vcf.gz')
    shell:
        """wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/paper_data_sets/companion_papers/mapping_structural_variation/union.2010_06.deletions.sites.vcf.gz -O {output.vcf}"""

# varset_set_mills2011_dl_novel
#
# Get novel sequence (non-reference) variants.
rule varset_set_mills2011_dl_novel:
    output:
        vcf=temp('temp/varset/set/mills2011/variants_novel.vcf.gz')
    shell:
        """wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/paper_data_sets/companion_papers/mapping_structural_variation/union.2010_06.novelsequences.sites.vcf.gz -O {output.vcf}"""

# varset_set_mills2011_dl_mei
#
# Get MEI variants.
rule varset_set_mills2011_dl_mei:
    output:
        vcf=temp('temp/varset/set/mills2011/variants_mei.vcf.gz')
    shell:
        """wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/paper_data_sets/companion_papers/mapping_structural_variation/union.2010_06.MobileElementInsertions.genotypes.vcf.gz -O {output.vcf}"""
