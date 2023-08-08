"""
Global rules for variants.
"""

# # variant_global_filter_fa
# rule variant_global_filter_fa:
#     input:
#         bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
#         fa='temp/variant/caller/{sourcename}/{sample}/all/all/bed/pre_filter/fa/{vartype}_{svtype}.fa.gz'
#     output:
#         fa='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz'
#     run:
#
#         # Read variant IDs
#         id_set = set(
#             pd.read_csv(input.bed, sep='\t', usecols=('ID', ), squeeze=True)
#         )
#
#         # Filter
#         with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
#             SeqIO.write(svpoplib.seq.fa_to_record_iter(input.fa, id_set), out_file, 'fasta')

def variant_global_filter_input_bed(wildcards):
    """
    Get input file for BED. Pull from caller if "filter" is all (pre-filter), pull from all for all filters.

    :param wildcards: Wildcards.

    :return: Path.
    """

    if wildcards.filter == 'all':
        return svpoplib.rules.parse_wildcards(
            'temp/variant/caller/{callertype}/{sourcename}/{sample}/all/all/bed/pre_filter/{vartype}_{svtype}.bed.gz',
            wildcards.sourcename, SAMPLE_TABLE, wildcards.sample, wildcards
        )

    return 'results/variant/caller/{sourcename}/{sample}/all/all/bed/{vartype}_{svtype}.bed.gz'.format(**wildcards)

def variant_global_filter_input_fa(wildcards):
    """
    Get input file for FASTA. Pull from caller if "filter" is all (pre-filter), pull from all for all filters.

    :param wildcards: Wildcards.

    :return: Path.
    """

    if wildcards.filter == 'all':
        return svpoplib.rules.parse_wildcards(
            'temp/variant/caller/{callertype}/{sourcename}/{sample}/all/all/bed/pre_filter/fa/{vartype}_{svtype}.fa.gz',
            wildcards.sourcename, SAMPLE_TABLE, wildcards.sample, wildcards
        )

    return 'results/variant/caller/{sourcename}/{sample}/all/all/bed/fa/{vartype}_{svtype}.fa.gz'.format(**wildcards)


# variant_global_filter_region
#
# Apply a BED filter to SVs.
rule variant_global_filter_region:
    input:
        bed=variant_global_filter_input_bed,
        fa=variant_global_filter_input_fa,
        filter=lambda wildcards: svpoplib.variant.get_filter_bed(wildcards.filter, UCSC_REF_NAME, config, SVPOP_DIR)
            if wildcards.filter != 'all' else []
    output:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        fa='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz',
        fai='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz.fai',
        bed_filt='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/filter_dropped/{vartype}_{svtype}_dropped.bed.gz'
    wildcard_constraints:
        vartype='sv|indel|snv|rgn|sub',
        svtype='ins|del|inv|snv|dup|rgn|sub'
    params:
        batch_size=10000
    run:

        # Read filter file
        empty_filter = False

        if wildcards.filter != 'all':

            with svpoplib.seq.PlainOrGzReader(input.bed_filt, 'rt') as in_file:
                try:
                    has_header = next(in_file).strip().startswith('#')
                except StopIteration:
                    empty_filter = True

            if has_header and not empty_filter:
                df_filter = pd.read_csv(input.bed_filt, sep='\t')

                missing_cols = [col for col in ('#CHROM', 'POS', 'END') if col not in set(df_filter.columns)]

                if missing_cols:
                    raise RuntimeError(f'Filter file is missing required columns: {", ".join(missing_cols)}')

                df_filter = df_filter[['#CHROM', 'POS', 'END']]

            elif not empty_filter:
                df_filter = pd.read_csv(input.bed_filt, sep='\t', header=None)

                if df_filter.shape[1] < 3:
                    raise RuntimeError(f'Filter file must have at least three columns: #CHROM, POS, and END (header line optional): {input.bed_filt}')

                df_filter = df_filter[[0, 1, 2]]
                df_filter.columns = ['#CHROM', 'POS', 'END']

            df_filter['#CHROM'] = df_filter['#CHROM'].astype(str)

        else:
            empty_filter = True

        # Get filter interval
        filter_tree = collections.defaultdict(intervaltree.intervaltree.IntervalTree)

        if not empty_filter:
            for index, row in df_filter.iterrows():
                filter_tree[row['#CHROM']][row['POS']:row['END']] = True

        # Process BED
        df_iter = pd.read_csv(
            input.bed, sep='\t',
            iterator=True,
            chunksize=params.batch_size
        )

        # Do filter
        chrom_set = set(svpoplib.ref.get_df_fai(config['reference_fai']).index)

        col_names = None

        write_header_pass = True
        write_header_filt = True

        id_set = set()  # Set of variant IDs passing filters

        with gzip.open(output.bed, 'wt') as out_file, gzip.open(output.bed_filt, 'wt') as out_file_filt:
            for df in df_iter:

                if col_names is None and df.shape[1] > 0:
                    col_names = list(df.columns)

                if df.shape[0] == 0:
                    continue

                filter_pass = df['#CHROM'].isin(chrom_set)
                filter_pass &= df.apply(lambda row: len(filter_tree[row['#CHROM']][row['POS']:row['END']]) == 0, axis=1)

                df_filt = df[~ filter_pass]
                df = df[filter_pass]

                if df.shape[0] > 0:
                    df.to_csv(out_file, sep='\t', index=False, header=write_header_pass)
                    write_header_pass = False

                    id_set |= set(df['ID'])

                if df_filt.shape[0] > 0:
                    df_filt.to_csv(out_file_filt, sep='\t', index=False, header=write_header_filt)
                    write_header_filt = False

            # Write empty BED files
            if col_names is None:
                col_names = ['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN'] + (['REF', 'ALT'] if wildcards.vartype == 'snv' else [])

            if write_header_pass:
                pd.DataFrame([], columns=col_names).to_csv(out_file, sep='\t', index=False, header=True)

            if write_header_filt:
                pd.DataFrame([], columns=col_names).to_csv(out_file_filt, sep='\t', index=False, header=True)

        # Create output FASTA
        if os.path.isfile(input.fa) and os.stat(input.fa).st_size > 0 and len(id_set) > 0:

            # Subset FASTA
            with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
                SeqIO.write(svpoplib.seq.fa_to_record_iter(input.fa, id_set, require_all=False), out_file, 'fasta')

            pysam.faidx(output.fa)

        else:
            # Write empty FASTA
            with open(output.fa, 'w') as out_file:
                pass

            with open(output.fai, 'w') as out_file:
                pass


# # variant_global_filter_region
# #
# # Apply a BED filter to SVs.
# rule variant_global_filter_region:
#     input:
#         bed='results/variant/caller/{sourcename}/{sample}/all/all/bed/{vartype}_{svtype}.bed.gz',
#         fa='results/variant/caller/{sourcename}/{sample}/all/all/bed/fa/{vartype}_{svtype}.fa.gz',
#         filter=lambda wildcards: svpoplib.variant.get_filter_bed(wildcards.filter, UCSC_REF_NAME, config, SVPOP_DIR)
#     output:
#         bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
#         fa='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz',
#         bed_filt='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/filter_dropped/{vartype}_{svtype}_dropped.bed.gz'
#     wildcard_constraints:
#         svtype='ins|del|inv|snv|dup|rgn|sub',
#         filter='^(?!all).+'
#     run:
#
#         # Filter
#         shell(
#             """bedtools intersect -wa -v -sorted -a {input.bed} -b {input.filter} -header | gzip > {output.bed}; """
#             """bedtools intersect -wa -u -sorted -a {input.bed} -b {input.filter} -header | gzip > {output.bed_filt}; """
#         )
#
#         # Create output FASTA
#         if os.path.isfile(input.fa) and os.stat(input.fa).st_size > 0:
#
#             # Subset FASTA
#             id_set = set(pd.read_csv(output.bed, sep='\t', usecols=('ID', ))['ID'])
#
#             with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
#                 SeqIO.write(svpoplib.seq.fa_to_record_iter(input.fa, id_set, require_all=False), out_file, 'fasta')
#
#         else:
#             # Write empty FASTA
#             with open(output.fa, 'w') as out_file:
#                 pass
#
#
# # variant_global_filter_region_all
# #
# # All variants before filter
# rule variant_global_filter_region_all:
#     input:
#         bed=lambda wildcards: svpoplib.rules.parse_wildcards(
#             'temp/variant/caller/{callertype}/{sourcename}/{sample}/all/all/bed/pre_filter/{vartype}_{svtype}.bed.gz',
#             wildcards.sourcename, SAMPLE_TABLE, wildcards.sample, wildcards
#         ),
#         fa=lambda wildcards: svpoplib.rules.parse_wildcards(
#             'temp/variant/caller/{callertype}/{sourcename}/{sample}/all/all/bed/pre_filter/fa/{vartype}_{svtype}.fa.gz',
#             wildcards.sourcename, SAMPLE_TABLE, wildcards.sample, wildcards
#         )
#     output:
#         bed='results/variant/caller/{sourcename}/{sample}/all/all/bed/{vartype}_{svtype}.bed.gz',
#         fa='results/variant/caller/{sourcename}/{sample}/all/all/bed/fa/{vartype}_{svtype}.fa.gz',
#         bed_filt='results/variant/caller/{sourcename}/{sample}/all/all/bed/filter_dropped/{vartype}_{svtype}_dropped.bed.gz'
#     wildcard_constraints:
#         svtype='ins|del|inv|snv|dup|rgn|sub'
#     run:
#
#         shell(
#             """cp -l {input.bed} {output.bed}; """
#             """touch {output.bed_filt}; """
#         )
#
#         # Create output FASTA
#         if os.path.isfile(input.fa) and os.stat(input.fa).st_size > 0:
#             shell(
#                 """cp -l {input.fa} {output.fa}"""
#             )
#
#         else:
#             # Write empty FASTA
#             with open(output.fa, 'w') as out_file:
#                 pass

# variant_global_vcf_gz
#
# Compress VCF files
# rule variant_global_vcf_gz:
#     input:
#         vcf='temp/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/vcf/{ref}_{fmt}/{vartype}_{svtype}.vcf'
#     output:
#         vcf='results/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/vcf/{ref}_{fmt}/{vartype}_{svtype}.vcf.gz',
#         tbi='results/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/vcf/{ref}_{fmt}/{vartype}_{svtype}.vcf.gz.tbi'
#     wildcard_constraints:
#         ref='grc|hc',
#         fmt='sv|alt'
#     shell:
#         """bgzip -c {input.vcf} > {output.vcf}; """
#         """sleep 5; """
#         """tabix {output.vcf}"""

# variant_global_uncompress_fa
#
# Uncompress for tools that cannot read a gzipped FASTA.
rule variant_global_uncompress_fa:
    input:
        fa='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz'
    output:
        fa=temp('temp/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa')
    run:

        if os.stat(input.fa).st_size > 0:
            shell("""zcat {input.fa} > {output.fa}""")
        else:
            shell("""> {output.fa}""")
