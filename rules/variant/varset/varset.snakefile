"""
Rules for processing and filtering varsets.
"""

###################
### Definitions ###
###################


#############
### Rules ###
#############

# # variant_varset_filter_region
# #
# # Apply a BED filter to SVs.
# rule variant_varset_filter_region:
#     input:
#         bed='temp/variant/varset/{varset}/bed_merged/{sample}/all/all/byref/{vartype}_{svtype}.bed',
#         filter=lambda wildcards: analib.variant.get_filter_bed(wildcards.filter, UCSC_REF_NAME, config, SVPOP_DIR)
#     output:
#         bed='results/variant/varset/{varset}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed'
#     wildcard_constraints:
#         vartype='sv|indel|snv|sub|rgn',
#         svtype='ins|del|inv|dup|snv|sub|rgn',
#         filter='((?!all)[^\/]*)|all.+'
#     run:
#
#         if wildcards.filter != 'all':
#             shell(
#                 """bedtools intersect -wa -v -sorted -a {input.bed} -b {input.filter} -header > {output.bed}; """
#             )
#         else:
#             shell(
#                 """cp {input.bed} {output.bed}; """
#             )

# variant_varset_bed
#
# Copy variant BED or merge SVs from from a merged variant set.
rule variant_varset_bed:
    input:
        bed=lambda wildcards: analib.varset.get_config_entry(wildcards.varset, wildcards, config)['input_bed'],
        filter=lambda wildcards: analib.variant.get_filter_bed(wildcards.filter, UCSC_REF_NAME, config, SVPOP_DIR)
    output:
        bed='results/variant/varset/{varset}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed'
    run:

        # Get config entry
        varset_entry = analib.varset.get_config_entry(wildcards.varset, wildcards, config)

        # Copy single source if not merged.
        if len(input.bed) == 1:

            if wildcards.filter != 'all':
                shell(
                    """bedtools intersect -wa -v -sorted -a {input.bed} -b {input.filter} -header > {output.bed}; """
                )
            else:
                shutil.copyfile(input.bed[0], output.bed)

            return

        # Merge
        df = analib.svmerge.merge_variants(
            varset_entry['input_bed'],
            varset_entry['varset'],
            strategy=varset_entry['merge_strategy'],
            threads=6
        )

        # Bylen to byref
        df['END'] = df.apply(lambda row: row['END'] if row['SVTYPE'] != 'INS' else row['POS'] + 1, axis=1)

        # Write
        df.to_csv(output.bed, sep='\t', index=False, na_rep='.')
