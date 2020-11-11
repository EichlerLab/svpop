"""
Get SVs from dbVAR.

https://www.ncbi.nlm.nih.gov/dbvar/
"""


#
# Definitions
#

# varset_set_dbvar_sv_bed
#
# Extract variants and do a non-redundant merge within dbVAR variants.
rule varset_set_dbvar_sv_bed:
    input:
        tab='temp/variant/varset/dbvar/variants.tab.gz',
        fai=config['reference_fai']
    output:
        bed=temp('temp/variant/varset/dbvar/bed/all/all/all/byref/sv_{svtype}.bed')
    wildcard_constraints:
        svtype='ins|del|inv|dup|cnv'
    run:

        # Column types for input data
        column_dtype = {
            '#CHROM': np.object,
            'POS': np.int32,
            'END': np.int32,
            'ID': np.object,
            'SVTYPE': np.object,
            'SVLEN': np.object,
            'ALT': np.object,
            'SAMPLE': np.object
        }

        # Contig name translation
        ref_contigs = list(pd.read_csv(input.fai, sep='\t', usecols=(0, ), header=None, squeeze=True))

        ref_contigs = set([contig for contig in ref_contigs if contig not in {'chrEBV', 'chrM'}])

        # Read
        df = pd.read_csv(input.tab, sep='\t',
            usecols=('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'ALT', 'SAMPLE'),
            dtype=column_dtype
        )

        # GRCh37 to hg38 and subset for know chromosomes
        df['#CHROM'] = analib.ref.grc_to_hg_chrom(df['#CHROM'], 'GRCh37', missing=True)
        df = df.loc[~ pd.isnull(df['#CHROM'])]  # Remove unknown chromosome names

        df['#CHROM'] = df['#CHROM'].apply(str)

        df = df.loc[df['#CHROM'].apply(lambda val: val in ref_contigs)]

        # Subset type and records
        df = df.loc[df['SVTYPE'].apply(lambda val: val == wildcards.svtype.upper())]

        if wildcards.svtype in {'inv', 'cnv'}:
            # Use END to calculate SVLEN
            df['SVLEN'] = df['END'] - df['POS']

        else:

            # Filter by integer SVLEN
            df = df.loc[df['SVLEN'].apply(analib.util.is_int)]

            # Set SVLEN and END
            df['SVLEN'] = np.abs(df['SVLEN'].astype(np.int))
            df['END'] = df['POS'] + df['SVLEN']

        # Subset SVs
        df = df.loc[df['SVLEN'] >= 50]

        # Set ID
        df['ORG_ID'] = df['ID']

        df['ID'] = analib.variant.get_variant_id(df)

        # Rearrange columns
        df = df.loc[: , ('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'ORG_ID', 'ALT', 'SAMPLE')]

        # Get a list of non-redundant dataframes (one per chromosome)
        df_nr_list = list()

        for chrom in sorted(set(df['#CHROM'])):
            df_nr_list.append(
                analib.variant.nr_interval_merge(df.loc[df['#CHROM'] == chrom], 0.50)
            )

        # Concatenate non-redundant variants from each chromosome to one dataframe
        df = pd.concat(df_nr_list, axis=0)

        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # To byref
        if wildcards.svtype == 'ins':
            df['END'] = df['POS'] + 1

        # Write
        df.to_csv(output.bed, sep='\t', index=False)


# varset_set_dbvar_sv_table
#
# VCF to variant table.
rule varset_set_dbvar_sv_table:
    input:
        vcf='temp/variant/varset/dbvar/variants.vcf.gz'
    output:
        tab=temp('temp/variant/varset/dbvar/variants.tab.gz')
    shell:
        """vcf2tsv {input.vcf} | """
        """gzip > {output.tab}"""

# varset_set_dbvar_dl
#
# Get variant VCF files.
rule varset_set_dbvar_dl:
    output:
        vcf=temp('temp/variant/varset/dbvar/variants.vcf.gz'),
        tbi=temp('temp/variant/varset/dbvar/variants.vcf.gz.tbi')
    params:
        url=config['varset']['set']['dbvar']['vcf']
    shell:
        """wget {params.url} -O {output.vcf}; """
        """wget {params.url}.tbi -O {output.tbi} """
