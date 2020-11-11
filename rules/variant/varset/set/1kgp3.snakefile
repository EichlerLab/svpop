"""
Get SVs from the 1000 Genomes Phase 3 project using the augmented Sudmant SV calls.

Citation:
Sudmant, P. H., Rausch, T., Gardner, E. J., Handsaker, R. E., Abyzov, A., Huddleston, J., … Korbel, J. O. (2015). An integrated map of structural variation in 2,504 human genomes. Nature, 526(7571), 75–81. http://doi.org/10.1038/nature15394
"""


#
# Definitions
#

def _varset_1kgp3_get_sv_type(row):
    """
    Determine the SVTYPE per row.

    :param row: Dataframe row (Series).
    """

    # DELs annotated as CN0
    if row['ALT'] == '<CN0>':
        return 'DEL'

    # CN2
    if row['ALT'] == '<CN2>':

        # CN2-DEL? SVTYPE is unclear
        if row['SVTYPE'] == 'DEL':
            return 'UNKNOWN'

        if row['SVTYPE'] == 'DUP':
            return 'DUP'

    # CN3-9
    if row['ALT'].startswith('<CN'):
        return 'DUP'

    # Pindel DELs
    if row['ALT'] in ('A', 'C', 'G', 'T'):
        return row['SVTYPE'].split('_', 1)[0]

    # MEI insertions
    if row['ALT'].startswith('<INS'):
        return 'INS'

    # Inversions
    if row['ALT'] == '<INV>':
        return 'INV'

    # SVTYPE could not be found
    raise ValueError('Cannot classify variant: ALT={ALT}, SVTYPE={SVTYPE}'.format(**row))


#
# Rules
#

# varset_set_1kgp3_table_to_bed
#
# Convert the variant table to BED and fix chromosome names for hg38 (from GRCh38).
rule varset_set_1kgp3_table_to_bed:
    input:
        tab='temp/variant/varset/1kgp3/variants.tab.gz'
    output:
        bed_ins=temp('temp/variant/varset/1kgp3/bed/all/all/all/byref/sv_ins.bed'),
        bed_del=temp('temp/variant/varset/1kgp3/bed/all/all/all/byref/sv_del.bed'),
        bed_inv=temp('temp/variant/varset/1kgp3/bed/all/all/all/byref/sv_inv.bed'),
        bed_cnv=temp('temp/variant/varset/1kgp3/bed/all/all/all/byref/sv_cnv.bed'),
        bed_dup=temp('temp/variant/varset/1kgp3/bed/all/all/all/byref/sv_dup.bed'),
        bed_uk='results/variant/varset/1kgp3/bed/all/all/all/byref/rm/sv_unknown.bed'
    run:

        # Read
        df = pd.read_csv(input.tab, sep='\t', header=0, dtype={'#CHROM': np.object})

        # GRCh38 to hg38 chromosome names
        df['#CHROM'] = analib.ref.grc_to_hg_chrom(df['#CHROM'], 'GRCh38')

        # Save extended SV type (ALU, CNV, etc)
        df['SVTYPE_ORG'] = df['SVTYPE']

        # Get SVTYPE
        df['SVTYPE'] = df.apply(_varset_1kgp3_get_sv_type, axis=1)

        df.loc[df['SVTYPE'] == 'UNKNOWN'].to_csv(output.bed_uk, sep='\t', index=False)

        df = df.loc[df['SVTYPE'] != 'UNKNOWN']


        # Set end
        df['SVLEN'] = df.apply(lambda row: abs(int(row['SVLEN'])) if row['SVLEN'] != '.' else (int(row['END']) - int(row['POS'])), axis=1)

        df['END'] = df.apply(
            lambda row: row['POS'] + (1 if row['SVTYPE'] == 'INS' else abs(int(row['SVLEN']))),
            axis=1
        )

        # Set ID
        df['ID'] = analib.variant.get_variant_id(df)

        # De-duplicate by ID
        df = df.groupby('ID').apply(lambda subdf: pd.Series(
            [
                subdf.iloc[0]['#CHROM'], subdf.iloc[0]['POS'], subdf.iloc[0]['END'],
                subdf.iloc[0]['ID'], subdf.iloc[0]['SVTYPE'], subdf.iloc[0]['SVLEN'],
                ','.join(subdf['ALT']), ','.join(subdf['SVTYPE_ORG']), ','.join(subdf['CS']),
                ','.join(subdf['QUAL'].apply(str)), ','.join(subdf['FILTER']),
                subdf.shape[0]
            ],
            index=['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'ALT', 'SVTYPE_ORG', 'CS', 'QUAL', 'FILTER', 'N_CALL']
        ))

        # Rearrange columns
        df = df.loc[:, ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'ALT', 'SVTYPE_ORG', 'CS', 'QUAL', 'FILTER')]

        # Sort
        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Write
        for svtype in ('INS', 'DEL', 'INV', 'DUP', 'CNV'):
            df.loc[df['SVTYPE'] == svtype].to_csv(output['bed_{}'.format(svtype.lower())], sep='\t', index=False)

# varset_set_1kgp3_variant_table
#
# Get a table of variants.
rule varset_set_1kgp3_variant_table:
    input:
        vcf='temp/variant/varset/1kgp3/variants.vcf.gz'
    output:
        tab=temp('temp/variant/varset/1kgp3/variants.tab.gz')
    shell:
        """zcat {input.vcf} | """
        """vcf2tsv | """
        """gzip > {output.tab}"""

# varset_set_1kgp3_dl_variants
#
# Get VCF.
rule varset_set_1kgp3_dl_variants:
    output:
        vcf=temp('temp/variant/varset/1kgp3/variants.vcf.gz')
    params:
        url=config['varset']['set']['1kgp3']['url']
    shell:
        """wget {params.url} -O {output.vcf}"""


#
# SVTYPE Translation
#
# The table below shows the ALT, SVTYPE, and CS (Caller) fields obtained from the VCF and the number of variant
# calls with each combination. The "NEW_SVTYPE" column shows how these are translated.
#
# Note: NEW_SVTYPE == DUP means Tandem Duplication.
#

# ALT	SVTYPE	CS	COUNT	NEW_SVTYPE
# <CN0>	CNV	DUP_gs	2306	DEL
# <CN0>	CNV	DUP_uwash	466	DEL
# <CN0>	DEL	DEL_union	31981	DEL
# <CN0>	DEL	DUP_gs	15	DEL
# <CN0>	DEL	DUP_uwash	1092	DEL
# <CN0>	DEL_ALU	DEL_union	2	DEL
# <CN0>	DEL_HERV	DEL_union	1	DEL
# <CN0>	DEL_LINE1	DEL_union	15	DEL
# <CN0>	DEL_SVA	DEL_union	1	DEL
# <CN0>	DUP	DUP_gs	306	DEL
# <CN2>	CNV	DUP_gs	2411	CNV
# <CN2>	CNV	DUP_uwash	488	CNV
# <CN2>	DEL	DUP_gs	15	UNKNOWN
# <CN2>	DEL	DUP_uwash	8	UNKNOWN
# <CN2>	DUP	DUP_delly	364	DUP
# <CN2>	DUP	DUP_gs	4413	CNV
# <CN2>	DUP	DUP_uwash	1229	CNV
# <CN3>	CNV	DUP_gs	250	CNV
# <CN3>	CNV	DUP_uwash	27	CNV
# <CN3>	DUP	DUP_gs	3	CNV
# <CN4>	CNV	DUP_gs	96	CNV
# <CN4>	CNV	DUP_uwash	3	CNV
# <CN5>	CNV	DUP_gs	46	CNV
# <CN6>	CNV	DUP_gs	16	CNV
# <CN7>	CNV	DUP_gs	6	CNV
# <CN8>	CNV	DUP_gs	2	CNV
# <CN9>	CNV	DUP_gs	1	CNV
# <INS:ME:ALU>	ALU	ALU_umary	12743	INS
# <INS:ME:LINE1>	LINE1	L1_umary	3047	INS
# <INS:ME:SVA>	SVA	SVA_umary	835	INS
# <INS:MT>	INS	NUMT_umich	162	INS
# <INV>	INV	CINV_delly	698	INV
# <INV>	INV	INV_delly	88	INV
# A	DEL	DEL_pindel	2258	DEL
# A	DEL_ALU	DEL_pindel	199	DEL
# A	DEL_LINE1	DEL_pindel	8	DEL
# C	DEL	DEL_pindel	1684	DEL
# C	DEL_ALU	DEL_pindel	132	DEL
# C	DEL_LINE1	DEL_pindel	6	DEL
# C	DEL_SVA	DEL_pindel	2	DEL
# G	DEL	DEL_pindel	1744	DEL
# G	DEL_ALU	DEL_pindel	118	DEL
# G	DEL_LINE1	DEL_pindel	1	DEL
# T	DEL	DEL_pindel	2148	DEL
# T	DEL_ALU	DEL_pindel	780	DEL
# T	DEL_LINE1	DEL_pindel	26	DEL
# T	DEL_SVA	DEL_pindel	6	DEL