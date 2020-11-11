"""
Merge annotations from each caller for a caller set.
"""

###################
### Definitions ###
###################


#############
### Rules ###
#############

# variant_callerset_anno_merge
#
# Merge annotations for a caller set.
rule variant_callerset_anno_merge:
    input:
        bed='results/variant/callerset/{callerset}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        tab=lambda wildcards: analib.callerset.get_caller_set_input(
            wildcards.callerset,
            'results/variant/{sourcetype}/{sourcename}/anno/{sample}/all/{filter}/{annodir}/{annotype}_{vartype}_{svtype}.{ext}',
            config,
            wildcards
        )
    output:
        tab='results/variant/callerset/{callerset}/anno/{sample}/all/{filter}/{annodir}/{annotype}_{vartype}_{svtype}.{ext}'
    params:
        mem=lambda wildcards: analib.callerset.cluster_param_anno_mem(wildcards, config),
    wildcard_constraints:
        filter='\\w+',
        annodir='[a-zA-Z0-9\\.\\-]+',
        sample='[a-zA-Z0-9\\.-]+',
        svtype='ins|del|inv|dup|sub|rgn|snv',
        ext='tab|bed'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Get sources.
        callerset_input = analib.callerset.get_caller_set_input(
            wildcards.callerset,
            'results/variant/{sourcetype}/{sourcename}/anno/{sample}/all/{filter}/{annodir}/{annotype}_{vartype}_{svtype}.{ext}',
            config,
            wildcards
        )

        # Get entry
        callerset_entry = analib.callerset.get_config_entry(wildcards.callerset, config)

        # Merge annotations
        df_merge = analib.callerset.merge_annotations(
            df, callerset_input, callerset_entry
        )

        # Write
        df_merge.to_csv(output.tab, sep='\t', index=False)
