"""
Comparing variants to a non-public set of SVs from Ira Hall. Do not share raw data
or results without expressed permission from Evan Eichler and/or Ira Hall.
"""

ruleorder: varset_set_hallsv_merge_all > global_variant_merge_all


#
# Definitions
#

def _data_varset_hallsv_bnd_cols(row):
    """
    Parse ALT field for BND records and return fields for this record.

    Fields:
        * BND_CHROM: Chromosome of other breakpoint.
        * BND_POS: Position of other breakpoint.
        * BND_STRAND: Strand of joined piece.
        * BND_REF: Reference bases at this breakpoint.
        * BND_ORIENT: Breakpoint orientation (POST: Piece joined after this one; PRE: piece joined before).
    """

    match_up = re.match('(.+)([\[\]])(\w+):(\d+)[\[\]]', row['ALT'])

    if match_up:
        return pd.Series(
            [match_up.group(3), match_up.group(4), ('+' if match_up.group(2) == '[' else '-'), match_up.group(1), 'POST'],
            index=['BND_CHROM', 'BND_POS', 'BND_STRAND', 'BND_REF', 'BND_ORIENT']
        )

    match_dn = re.match('([\[\]])(\w+):(\d+)[\[\]](.+)', row['ALT'])

    if match_dn:
        return pd.Series(
            [match_dn.group(2), match_dn.group(3), ('+' if match_dn.group(1) == ']' else '-'), match_dn.group(4), 'PRE'],
            index=['BND_CHROM', 'BND_POS', 'BND_STRAND', 'BND_REF', 'BND_ORIENT']
        )

    raise RuntimeError('Cannot match BND="{}"'.format(row['ALT']))


#
# Rules
#

# varset_set_hallsv_sv_types
#
# Get SV by type.
rule varset_set_hallsv_sv_types:
    input:
        bed='results/varset/set/hallsv_{sample}/all/all/byref/sv_all.bed'
    output:
        bed_ins='results/varset/set/hallsv_{sample}/all/all/byref/sv_ins.bed',
        bed_del='results/varset/set/hallsv_{sample}/all/all/byref/sv_del.bed',
        bed_inv='results/varset/set/hallsv_{sample}/all/all/byref/sv_inv.bed'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Write
        df.loc[df['SVTYPE'] == 'INS'].to_csv(output.bed_ins, sep='\t', index=False)
        df.loc[df['SVTYPE'] == 'DEL'].to_csv(output.bed_del, sep='\t', index=False)
        df.loc[df['SVTYPE'] == 'INV'].to_csv(output.bed_inv, sep='\t', index=False)

# varset_set_hallsv_merge_all
#
# Merge SVs into one.
#
# NOTE: BND DELs are not included (commented out).
rule varset_set_hallsv_merge_all:
    input:
        sv_other='temp/varset/set/hallsv/sv_other_{sample}.bed'
#        sv_bnd_del='temp/varset/set/hallsv/bnd_del_{sample}.bed'
    output:
        bed=temp('results/varset/set/hallsv_{sample}/all/all/byref/sv_all.bed')
    shell:
        """head -n 1 {input.sv_other} """
        """> {output.bed}; """
        """awk 'FNR > 1' {input} | """
        """sort -k1,1 -k2,2n """
        """>> {output.bed}"""

rule varset_set_hallsv_get_mei:
    input:
        tab='temp/varset/set/hallsv/variants_{sample}.tab'
    output:
        bed='results/varset/set/hallsv_{sample}/all/all/byref/sv_dup.bed'
    run:

        # Read
        df = pd.read_csv(input.tab, sep='\t', header=0, low_memory=False, na_values=['', '.'])
        df = df.loc[df['SVTYPE'] == 'DUP']

        df = df.loc[df['GT'].apply(lambda val: val not in {'./.', '0/0'})]
        df = df.loc[df['FILTER'] == 'PASS']

        # Set fields
        df['POS1'] = df['POS']
        df['POS'] -= 1

        df['SVLEN'] = np.abs(df['SVLEN']).astype(np.int)

        df['END'] = df['END'].astype(np.int)

        df['END'] = df.apply(lambda row: row['POS'] + (1 if row['SVTYPE'] == 'INS' else row['SVLEN']), axis=1)

        df['ID'] = df.apply(lambda row: '{#CHROM}-{POS1}-{SVTYPE}-{SVLEN}'.format(**row), axis=1)

        # Filter for SVs
        df = df.loc[df['SVLEN'] >= 50]

        # Sort
        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Write
        df.loc[:, ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN')].to_csv(output.bed, sep='\t', index=False)

# varset_set_hallsv_other
#
# Get INS (MEI), DEL, and INV records.
rule varset_set_hallsv_other:
    input:
        tab='temp/varset/set/hallsv/variants_{sample}.tab'
    output:
        bed=temp('temp/varset/set/hallsv/sv_other_{sample}.bed')
    run:

        # Read
        df = pd.read_csv(input.tab, sep='\t', header=0, low_memory=False, na_values=['', '.'])
        #df = df.loc[df['SVTYPE'].apply(lambda svtype: svtype in {'DEL', 'INV', 'MEI'})]
        df = df.loc[df['SVTYPE'] != 'BND']

        df = df.loc[df['GT'].apply(lambda val: val not in {'./.', '0/0'})]
        df = df.loc[df['FILTER'] == 'PASS']

        # Set fields
        #df['SVTYPE'] = df['SVTYPE'].apply(lambda svtype: svtype if svtype != 'MEI' else 'INS')
        df['ALT'].apply(lambda val: re.sub('<(\w+)[:>].*', '\\1', val))

        df['POS1'] = df['POS']
        df['POS'] -= 1

        df['SVLEN'] = np.abs(df['SVLEN']).astype(np.int)

        df['END'] = df['END'].astype(np.int)

        df['END'] = df.apply(lambda row: row['POS'] + (1 if row['SVTYPE'] == 'INS' else row['SVLEN']), axis=1)

        df['ID'] = df.apply(lambda row: '{#CHROM}-{POS1}-{SVTYPE}-{SVLEN}'.format(**row), axis=1)

        # Filter for SVs
        df = df.loc[df['SVLEN'] >= 50]

        # Write
        df.loc[:, ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN')].to_csv(output.bed, sep='\t', index=False)


# varset_set_hallsv_bnd_del
#
# Get deletions using breakend notation (BND VCF records).
rule varset_set_hallsv_bnd_del:
    input:
        tab=temp('temp/varset/set/hallsv/variants_{sample}.tab')
    output:
        bed=temp('temp/varset/set/hallsv/bnd_del_{sample}.bed')
    run:

        # Read BND records
        df = pd.read_csv(input.tab, sep='\t', header=0, low_memory=False, na_values=['', '.'])
        df = df.loc[df['SVTYPE'] == 'BND']

        df = df.loc[df['GT'].apply(lambda val: val not in {'./.', '0/0'})]

        # Get paired BND records (ignore singleton, triplicate, etc)
        bnd2_key = df.groupby('EVENT').apply(lambda subdf: subdf.shape[0] == 2)
        df = df.loc[df['EVENT'].apply(lambda val: val in bnd2_key.index[bnd2_key])]

        del(bnd2_key)

        # Parse breakpoint information into fields
        df = pd.concat([df, df.apply(_data_varset_hallsv_bnd_cols, axis=1)], axis=1)

        df['BND_POS'] = df['BND_POS'].astype(np.int)


        # Filter for deletions (left-most record, right-most record is not needed)
        df = df.loc[(df['BND_ORIENT'] == 'POST') & (df['BND_STRAND'] == '+')]
        df = df.loc[(df['#CHROM'] == df['BND_CHROM']) & (df['POS'] < df['BND_POS'])]

        # Make SV BED
        df['POS1'] = df['POS']
        df['POS'] -= 1
        df['END'] = df['BND_POS']

        df['SVLEN'] = df['END'] - df['POS']
        df['SVTYPE'] = 'DEL'

        df['ID'] = df.apply(lambda row: '{#CHROM}-{POS1}-DEL-{SVLEN}'.format(**row), axis=1)

        # Filter for SVs
        df = df.loc[df['SVLEN'] >= 50]

        # Write
        df.loc[:, ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN')].to_csv(output.bed, sep='\t', index=False)

# varset_set_hallsv_extract_sample
#
# Extract sample table from VCF.
rule varset_set_hallsv_extract_sample:
    input:
        vcf=config['varset']['set']['hallsv']['vcf']
    output:
        tab=temp('temp/varset/set/hallsv/variants_{sample}.tab')
    shell:
        """vcfkeepsamples {input.vcf} {wildcards.sample} | """
        """vcf2tsv -n "" -g """
        """>{output.tab}"""
