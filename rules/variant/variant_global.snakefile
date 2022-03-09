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
        bed_filt='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/filter_dropped/{vartype}_{svtype}_dropped.bed.gz'
    wildcard_constraints:
        svtype='ins|del|inv|snv|dup|rgn|sub'
    run:

        if wildcards.filter == 'all':
            # Pull from base caller temp

            shell(
                """cp {input.bed} {output.bed}; """
                """touch {output.bed_filt}; """
            )

            # Create output FASTA
            if os.path.isfile(input.fa) and os.stat(input.fa).st_size > 0:
                shell(
                    """cp {input.fa} {output.fa}"""
                )

            else:
                # Write empty FASTA
                with open(output.fa, 'w') as out_file:
                    pass

        else:
            # Filter from all

            shell(
                """bedtools intersect -wa -v -sorted -a {input.bed} -b {input.filter} -header | gzip > {output.bed}; """
                """bedtools intersect -wa -u -sorted -a {input.bed} -b {input.filter} -header | gzip > {output.bed_filt}; """
            )

            # Create output FASTA
            if os.path.isfile(input.fa) and os.stat(input.fa).st_size > 0:

                # Subset FASTA
                id_set = set(pd.read_csv(output.bed, sep='\t', usecols=('ID', ))['ID'])

                with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
                    SeqIO.write(svpoplib.seq.fa_to_record_iter(input.fa, id_set, require_all=False), out_file, 'fasta')

            else:
                # Write empty FASTA
                with open(output.fa, 'w') as out_file:
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
