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
        bed='results/variant/sampleset/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        tab=lambda wildcards: svpoplib.sampleset.get_sample_set_input(
            wildcards.sourcename,
            wildcards.sample,
            'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/{annodir}/{annotype}_{vartype}_{svtype}.{ext}.gz',
            config,
            wildcards
        )
    output:
        tab='results/variant/sampleset/{sourcename}/{sample}/{filter}/all/anno/{annodir}/{annotype}_{vartype}_{svtype}.{ext}.gz'
    params:
        mem=lambda wildcards: svpoplib.sampleset.cluster_param_anno_mem(wildcards, config)
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|rgn|snv',
        ext='tsv|bed'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Get sources.
        sampleset_input = svpoplib.sampleset.get_sample_set_input(
            wildcards.sourcename,
            wildcards.sample,
            'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/{annodir}/{annotype}_{vartype}_{svtype}.{ext}.gz',
            config,
            wildcards
        )

        # Get entry
        sampleset_entry = svpoplib.sampleset.get_config_entry(wildcards.sourcename, wildcards.sample, config)

        # Merge annotations
        df_merge = svpoplib.sampleset.merge_annotations(
            df, sampleset_input, sampleset_entry
        )

        # Write
        df_merge.to_csv(output.tab, sep='\t', index=False, compression='gzip')
