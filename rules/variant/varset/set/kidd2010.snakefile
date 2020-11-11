"""
SV Calls published by Kidd (2010). 1054 large SVs with capillary end-sequenced FOSMID clones.

Kidd, J. M., Graves, T., Newman, T. L., Fulton, R., Hayden, H. S., Malig, M., … Eichler, E. E. (2010). A human genome structural variation sequencing resource reveals insights into mutational mechanisms. Cell, 143(5), 837–847. http://doi.org/10.1016/j.cell.2010.10.027
"""


# varset_kidd2010_demux_svtype
#
# Separate into ins/del/inv.
rule varset_kidd2010_demux_svtype:
    input:
        bed='temp/varset/set/kidd2010/byref/sv_all.bed'
    output:
        bed_ins='results/varset/set/kidd2010/all/all/byref/sv_ins.bed',
        bed_del='results/varset/set/kidd2010/all/all/byref/sv_del.bed',
        bed_inv='results/varset/set/kidd2010/all/all/byref/sv_inv.bed'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Write
        df.loc[df['SVTYPE'] == 'INS'].to_csv(output.bed_ins, sep='\t', index=False)
        df.loc[df['SVTYPE'] == 'DEL'].to_csv(output.bed_del, sep='\t', index=False)
        df.loc[df['SVTYPE'] == 'INV'].to_csv(output.bed_inv, sep='\t', index=False)


# varset_kidd2010_liftover
#
# Liftover variants to hg38.
rule varset_kidd2010_liftover:
    input:
        bed='temp/varset/set/kidd2010/hg18/sv_all.bed',
        chain='data/liftover/hg18_to_hg38.chain.gz'
    output:
        bed=temp('temp/varset/set/kidd2010/byref/sv_all.bed'),
        temp=temp('temp/varset/kidd2010/byref/hg18/sv_all_temp.bed'),
        unmapped='results/varset/set/kidd2010/all/all/byref/hg18/unmapped_all.bed'
    shell:
        """liftOver -bedPlus=3 {input.bed} {input.chain} {output.temp} {output.unmapped}; """
        """head -n 1 {input.bed} >{output.bed}; """
        """sort -k1,1 -k 2,2n {output.temp} >>{output.bed}"""

# varset_kidd2010_merge_calls
#
# Merge calls. Indel, VNTR, and INV were given on different tabs of the supplemental spreadsheet S2.
rule varset_kidd2010_merge_calls:
    input:
        indel=config['varset']['set']['kidd2010']['indel'],
        vntr=config['varset']['set']['kidd2010']['vntr'],
        inv=config['varset']['set']['kidd2010']['inv']
    output:
        bed=temp('temp/varset/set/kidd2010/hg18/sv_all.bed')
    run:

        # Declarations
        kidd_type_to_svtype = {
            'I': 'INS',
            'D': 'DEL',
            'V': 'INV'
        }

        kidd_cols = (
            'clone name', 'variant accession', 'chrm', 'variant type', 'clone break 1', 'clone break 2',
            'chrm break 1', 'chrm break 2', 'classification', 'junction size', 'part of nonredundant set'
        )

        # Read
        df_indel = pd.read_csv(input.indel, sep='\t', header=0, usecols=kidd_cols)
        df_vntr = pd.read_csv(input.vntr, sep='\t', header=0)
        df_inv = pd.read_csv(input.inv, sep='\t', header=0, usecols=kidd_cols)

        df_vntr['classification'] = np.nan
        df_vntr['junction size'] = np.nan
        df_vntr = df_vntr.loc[:, kidd_cols]

        df_indel['DISCO_SET'] = 'INDEL'
        df_vntr['DISCO_SET'] = 'VNTR'
        df_inv['DISCO_SET'] = 'INV'

        # Merge
        df = pd.concat([df_indel, df_vntr, df_inv], axis=0)
        df.reset_index(inplace=True)

        # Filter for non-redundant variants
        df = df.loc[df['part of nonredundant set'] == 'part of nonredundant set']

        # Breakpoints to BED-like cooridnates (0-based, half-open)
        df['chrm break 1'] = df['chrm break 1'] - 1
        df['clone break 1'] = df['clone break 1'] - 1

        # Set SVTYPE
        df['SVTYPE'] = df['variant type'].apply(lambda val: kidd_type_to_svtype[val])

        # Get SVLEN
        df['SVLEN'] = df.apply(
            lambda row: row['clone break 2'] - row['clone break 1'] if row['SVTYPE'] == 'INS' else row['chrm break 2'] - row['chrm break 1'],
            axis=1
        )

        # Make BED columns
        df['#CHROM'] = df['chrm']
        df['POS'] = df['chrm break 1']
        df['END'] = df.apply(lambda row: row['POS'] + (1 if row['SVTYPE'] == 'INS' else row['SVLEN']), axis=1)

        # Set ID
        df['POS1'] = df['POS'] + 1
        df['ID'] = df.apply(lambda row: '{#CHROM}-{POS1}-{SVTYPE}-{SVLEN}'.format(**row), axis=1)
        del(df['POS1'])

        # Set info
        df['INFO'] = df.apply(lambda row: 'set={DISCO_SET};acc={variant accession};class={classification};junction={junction size};clone={clone name}'.format(**row), axis=1)

        df = df.loc[:, ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'INFO')]

        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Write
        df.to_csv(output.bed, sep='\t', index=False, na_rep='NA')
