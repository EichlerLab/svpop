"""
SGDP variants.

Sudmant, P. H., Mallick, S., Nelson, B. J., Hormozdiari, F., Krumm, N., Huddleston, J., … Eichler, E. E. (2015). Global diversity, population stratification, and selection of human copy number variation. Science (New York, N.Y.), 1–16. http://doi.org/10.1126/science.aab3761
"""

def _varset_sudmant2015a_get_sv_type(row):

    if row['ALT'] == '<CN0>':
        return 'DEL'

    if row['ALT'].startswith('<CN'):
        return 'CNV'

    if row['ALT'] in ('A', 'C', 'G', 'T'):
        return row['SVTYPE_ORG'].split('_', 1)[0]

    if row['ALT'].startswith('<INS'):
        return 'INS'

    if row['ALT'] == '<INV>':
        return 'INV'

    raise ValueError('Cannot classify variant: ALT={ALT}, SVTYPE={SVTYPE_ORG}'.format(**row))


# varset_set_sudmant2015a_empty
#
# Write empty files for others.
rule varset_set_sudmant2015a_empty:
    output:
        bed=temp('temp/variant/varset/sudmant2015a/bed/all/all/all/byref/sv_{svtype}.bed')
    wildcard_constraints:
        svtype='ins|inv'
    shell:
        """echo -e "#CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN" >{output.bed}; """

# varset_set_sudmant2015a_annotate
#
# Annotate variants.
rule varset_set_sudmant2015a_annotate:
    input:
        bed='temp/variant/varset/sudmant2015a/byref/sv_{svtype}.bed'
    output:
        bed=temp('temp/variant/varset/sudmant2015a/bed/all/all/all/byref/sv_{svtype}.bed')
    wildcard_constraints:
        svtype='del|dup'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Set ID
        df['ID'] = analib.variant.get_variant_id(df)

        # Write
        df = df.loc[:, ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN')]

        df.to_csv(output[0], sep='\t', index=False, header=True)

# varset_set_sudmant2015a_liftover
#
# Liftover variants to hg38.
rule varset_set_sudmant2015a_liftover:
    input:
        bed='temp/variant/varset/sudmant2015a/hg19/sv_{svtype}.bed',
        chain='data/liftover/hg19_to_hg38.chain.gz',
    output:
        bed=temp('temp/variant/varset/sudmant2015a/byref/sv_{svtype}.bed'),
        unmapped='results/variant/varset/sudmant2015a/bed/all/all/all/dropped/liftover_unmapped_{svtype}.bed'
    shell:
        """{{\n"""
        """    head -n 1 {input.bed}; """
        """    liftOver -bedPlus=4 {input.bed} {input.chain} /dev/stdout {output.unmapped} | sort -k1,1 -k 2,2n; """
        """}} """
        """> {output.bed}"""

# varset_set_sudmant2015a_get_del
#
# Get deletions.
rule varset_set_sudmant2015a_get_del:
    input:
        tab='files/varset/sudmant2015a/Sudmant2015a.tab'
    output:
        bed=temp('temp/variant/varset/sudmant2015a/hg19/sv_del.bed'),
        dup=temp('temp/variant/varset/sudmant2015a/hg19/sv_dup.bed')
    run:

        # Read
        df = pd.read_csv(input.tab, sep='\t', header=0, usecols=('contig', 'start', 'end', 'size', 'type'))

        # Get columns and large DELs
        df = df.loc[:, ('contig', 'start', 'end', 'size', 'type')]

        df.columns = ['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE']

        # Write dups
        df.loc[df['SVTYPE'] != 'DEL'].to_csv(output.dup, sep='\t', index=False)

        # Get DELs
        df = df.loc[df['SVTYPE'] == 'DEL']

        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Write
        df.to_csv(output.bed, sep='\t', index=False)
