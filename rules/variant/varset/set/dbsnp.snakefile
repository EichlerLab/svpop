"""
Get variants from dbSNP.
"""

# varset_set_dbsnp_sv_bed
#
# Extract variants and do a non-redundant merge within dbVAR variants.
rule varset_set_dbsnp_sv_bed:
    input:
        tab='temp/variant/varset/dbsnp151/variants.tab.gz',
        ref=config['reference']
    output:
        sv_ins=temp('temp/variant/varset/dbsnp151/bed/all/all/all/byref/sv_ins.bed')
    run:

        # Column types for input data
        column_dtype = {
            '#CHROM': np.object,
            'POS': np.int32,
            'ID': np.object,
            'REF': np.object,
            'ALT': np.object,
            'SAMPLE': np.object,
            'AC': np.float32,
            'AF': np.float32,
            'AN': np.float32,
            'CIPOS': np.object,
            'CIEND': np.object,
            'SVTYPE': np.object
        }

        # Contig name translation
        ref_contigs = list(get_df_fai(input.ref + '.fai').index)

        ref_contigs = set([contig for contig in ref_contigs if contig not in {'chrEBV', 'chrM'}])

        # Read
        df = pd.read_csv(input.tab, sep='\t',
            usecols=('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'SAMPLE', 'AC', 'AF', 'AN', 'CIPOS', 'CIEND', 'SVTYPE'),
            dtype=column_dtype
        )

        df['DBSNPID'] = df['ID']
        del(df['ID'])

        # GRCh37 to hg38 and subset for know chromosomes
        df['#CHROM'] = analib.ref.grc_to_hg_chrom(df['#CHROM'], 'GRCh37', missing=True)
        df = df.loc[~ pd.isnull(df['#CHROM'])]  # Remove unknown chromosome names

        df['#CHROM'] = df['#CHROM'].apply(str)

        df = df.loc[df['#CHROM'].apply(lambda val: val in ref_contigs)]

        # Get indel record




        # Subset type and records


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

# varset_set_dbsnp_sv_table
#
# VCF to variant table.
rule varset_set_dbsnp_sv_table:
    input:
        vcf='temp/variant/varset/dbsnp151/variants.vcf.gz'
    output:
        tab=temp('temp/variant/varset/dbsnp151/variants.tab.gz')
    shell:
        """vcf2tsv {input.vcf} | """
        """gzip > {output.tab}"""

# varset_set_dbsnp_dl
#
# Get variant VCF files.
rule varset_set_dbsnp_dl:
    output:
        vcf=temp('temp/variant/varset/dbsnp151/variants.vcf.gz'),
        tbi=temp('temp/variant/varset/dbsnp151/variants.vcf.gz.tbi')
    params:
        url=config['varset']['set']['dbsnp151']['vcf']
    shell:
        """wget {params.url} -O {output.vcf}; """
        """wget {params.url}.tbi -O {output.tbi} """
