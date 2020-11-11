"""
Merge variants from multiple callers for a single sample.
"""

###################
### Definitions ###
###################



#############
### Rules ###
#############

# variant_callerset_bed_merge
#
# Merge callerset variants from multiple sources for one sample.
rule variant_callerset_bed_merge:
    input:
        bed=expand(
            'temp/variant/callerset/{{callerset}}/bed/{{sample}}/all/{{filter}}/byref/{{vartype}}_{{svtype}}/chrom_{chrom}.bed',
            chrom=analib.ref.get_df_fai(config['reference_fai']).index
        )
    output:
        bed='results/variant/callerset/{callerset}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed'
    wildcard_constraints:
        vartype='sv|indel|snv',
        svtype='ins|del|inv|dup|snv',
        sample='[^/]+'
    run:

        df = pd.concat(
            [pd.read_csv(bed_file_name, sep='\t') for bed_file_name in input.bed],
            sort=False,
            axis=0
        ).sort_values(
            ['#CHROM', 'POS']
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False)


# variant_callerset_bed_merge_chrom
#
# Merge one chromosome.
rule variant_callerset_bed_merge_chrom:
    input:
        bed=lambda wildcards: analib.callerset.get_caller_set_input(
            wildcards.callerset,
            'results/variant/{sourcetype}/{sourcename}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
            config,
            wildcards
        )
    output:
        bed=temp('temp/variant/callerset/{callerset}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}/chrom_{chrom}.bed')
    params:
        cpu=lambda wildcards: analib.callerset.cluster_param_cpu(wildcards, config),
        mem=lambda wildcards: analib.callerset.cluster_param_mem(wildcards, config),
        rt=lambda wildcards: analib.callerset.cluster_param_rt(wildcards, config)
    run:

        # Get entry
        callerset_entry = analib.callerset.get_config_entry(wildcards.callerset, config)

        merge_strategy = analib.sampleset.get_merge_strategy(callerset_entry, wildcards.vartype, wildcards.svtype)

        # Merge
        df = analib.svmerge.merge_variants(
            input.bed,
            callerset_entry['name_list'],
            merge_strategy['strategy'],
            subset_chrom=wildcards.chrom,
            threads=params.cpu
        )

        # Rename MERGE for callerset
        rename_dict = {
            'MERGE_SRC': 'CALLERSET_SRC',
            'MERGE_SRC_ID': 'CALLERSET_SRC_ID',
            'MERGE_AC': 'CALLERSET_N',
            'MERGE_AF': 'CALLERSET_PROP',
            'MERGE_SAMPLES': 'CALLERSET_LIST',
            'MERGE_VARIANTS': 'CALLERSET_VARIANTS',
            'MERGE_RO': 'CALLERSET_RO',
            'MERGE_DIST': 'CALLERSET_DIST',
            'MERGE_SZRO': 'CALLERSET_SZRO',
            'MERGE_OFFSZ': 'CALLERSET_OFFSZ'
        }

        df.columns = [rename_dict.get(col, col) for col in df.columns]

        if 'DISC_CLASS' in df.columns:
            del(df['DISC_CLASS'])

        # Bylen to byref
        if df.shape[0] > 0:
            df['END'] = df.apply(lambda row: (row['POS'] + 1) if row['SVTYPE'] == 'INS' else row['END'], axis=1)

        # Write
        df.to_csv(output.bed, sep='\t', index=False)