"""
Merge annotations from a merged set of samples.
"""

import pandas as pd

import svpoplib

#
# Rules
#

# Merge annotations for a sample set.
rule variant_sampleset_anno_merge:
    input:
        bed='results/variant/sampleset/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        tsv=lambda wildcards: svpoplib.sampleset.get_sample_set_input(
            wildcards.sourcename,
            wildcards.sample,
            'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/{annodir}/{annotype}_{vartype}_{svtype}.{ext}.gz',
            config,
            wildcards
        )
    output:
        tsv='results/variant/sampleset/{sourcename}/{sample}/{filter}/all/anno/{annodir}/{annotype}_{vartype}_{svtype}.{ext}.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|rgn|snv',
        vartype='sv|indel|snv|sub|rgn',
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
        df_merge.to_csv(output.tsv, sep='\t', index=False, compression='gzip')
