"""
AK1 assembly SV calls.

Seo, J.-S., Rhie, A., Kim, J., Lee, S., Sohn, M.-H., Kim, C.-U., … Kim, C. (2016). De novo assembly and phasing of a Korean human genome. Nature, 538(7624), 243–247. https://doi.org/10.1038/nature20098

Input file is Supplementary Table 16 (GRCh38 SV calls)
"""

# varset_set_ak1_svtype
#
# Get BED by SVTYPE.
rule varset_set_ak1_svtype:
    input:
        bed='temp/variant/varset/ak1/bed/all/all/all/byref/sv_all.bed'
    output:
        bed='temp/variant/varset/ak1/bed/all/all/all/byref/sv_{svtype}.bed'
    wildcard_constraints:
        svtype='ins|del|inv'
    run:
        df = pd.read_csv(input.bed, sep='\t', header=0)
        df = df.loc[df['SVTYPE'] == wildcards.svtype.upper()]
        df.to_csv(output.bed, sep='\t', index=False)


# varset_set_ak1_fixbed
#
# Fix AK1 variant BED file.
rule varset_set_ak1_all:
    input:
        bed=os.path.join(SVPOP_DIR, 'files/varset/ak1/ak1_sv_hg38.tab.gz')
    output:
        bed='temp/variant/varset/ak1/bed/all/all/all/byref/sv_all.bed',
        bed_complex='results/variant/varset/ak1/bed/all/all/all/byref/rm/sv_complex.bed'
    run:

        # Read
        df = pd.read_csv(
            input.bed,
            sep='\t',
            header=0,
            usecols=('Chromosome', 'Start', 'End', 'Type', 'Size', 'Source')
        )

        df.columns = ('#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'SOURCE')

        # Fix coordinates
        df['POS'] = df.apply(lambda row: row['POS'] - (2 if row['SVTYPE'] == 'INS' else 1), axis=1)
        df['END'] = df['END'] - 1

        # Set ID
        df['POS1'] = df['POS'] + 1
        df['ID'] = df.apply(lambda row: '{#CHROM}-{POS1}-{SVTYPE}-{SVLEN}'.format(**row), axis=1)
        del(df['POS1'])

        # Subset (cannot compare complex variants)
        df.loc[df['SVTYPE'] == 'COMPLEX'].to_csv(output.bed_complex, sep='\t', index=False)
        df = df.loc[df['SVTYPE'] != 'COMPLEX']

        # Collapse duplications (same variant with different annotations - e.g. intersects multiple genes)
        df_agg = df.groupby('ID').agg(lambda val: val.iloc[0])
        df_agg['ID'] = df_agg.index

        # Rearrange columns
        df_agg = df_agg.loc[:, ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'SOURCE')]

        # Sort
        df_agg.sort_values(['#CHROM', 'POS'], inplace=True)

        # Write
        df_agg.to_csv(output.bed, sep='\t', index=False)
