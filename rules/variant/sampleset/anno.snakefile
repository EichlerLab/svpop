"""
Merge annotations from a merged set of samples.
"""

###################
### Definitions ###
###################


#############
### Rules ###
#############

# variant_sampleset_anno_merge
#
# Merge annotations for a sample set.
rule variant_sampleset_anno_merge:
    input:
        bed='results/variant/sampleset/{sampleset}/bed/{samplelist}/all/{filter}/byref/{vartype}_{svtype}.bed',
        tab=lambda wildcards: analib.sampleset.get_sample_set_input(
            wildcards.sampleset,
            wildcards.samplelist,
            'results/variant/{sourcetype}/{sourcename}/anno/{sample}/all/{filter}/{annodir}/{annotype}_{vartype}_{svtype}.{ext}',
            config,
            wildcards
        )
    output:
        tab='results/variant/sampleset/{sampleset}/anno/{samplelist}/all/{filter}/{annodir}/{annotype}_{vartype}_{svtype}.{ext}'
    params:
        mem=lambda wildcards: analib.sampleset.cluster_param_anno_mem(wildcards, config),
    wildcard_constraints:
        filter='\\w+',
        svtype='ins|del|inv|dup|sub|rgn|snv',
        sample='[a-zA-Z0-9\\.]+',
        ext='tab|bed'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Get sources.
        sampleset_input = analib.sampleset.get_sample_set_input(
            wildcards.sampleset,
            wildcards.samplelist,
            'results/variant/{sourcetype}/{sourcename}/anno/{sample}/all/{filter}/{annodir}/{annotype}_{vartype}_{svtype}.{ext}',
            config,
            wildcards
        )

        # Get entry
        sampleset_entry = analib.sampleset.get_config_entry(wildcards.sampleset, wildcards.samplelist, config)

        # Merge annotations
        df_merge = analib.sampleset.merge_annotations(
            df, sampleset_input, sampleset_entry
        )

        # Write
        df_merge.to_csv(output.tab, sep='\t', index=False)
