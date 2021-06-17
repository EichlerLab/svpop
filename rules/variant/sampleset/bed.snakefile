"""
Merge variants from multiple samples.
"""

#############
### Rules ###
#############

# variant_sampleset_bed_merge
#
# Merge callerset variants from multiple sources for one sample.
rule variant_sampleset_bed_merge:
    input:
        bed=expand(
            'temp/variant/sampleset/{{sourcename}}/{{sample}}/{{filter}}/all/bed/{{vartype}}_{{svtype}}/chrom_{chrom}.bed.gz',
            chrom=svpoplib.ref.get_df_fai(config['reference_fai']).index
        )
    output:
        bed='results/variant/sampleset/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz'
    wildcard_constraints:
        vartype='sv|indel|snv',
        svtype='ins|del|inv|dup|snv',
        samplelist='[^/]+'
    params:
        mem=lambda wildcards: svpoplib.sampleset.cluster_param_anno_mem(wildcards, config),
    run:

        df = pd.concat(
            [pd.read_csv(bed_file_name, sep='\t') for bed_file_name in input.bed],
            axis=0,
            sort=False
        ).sort_values(
            ['#CHROM', 'POS']
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# variant_sampleset_bed_merge_chrom
#
# Get a set of variants merged from multiple samples.
rule variant_sampleset_bed_merge_chrom:
    input:
        bed=lambda wildcards: svpoplib.sampleset.get_sample_set_input(
            wildcards.sourcename,
            wildcards.sample,
            'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
            config,
            wildcards
        )
    output:
        bed=temp(
            'temp/variant/sampleset/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}/chrom_{chrom}.bed.gz'
        )
    params:
        cpu=lambda wildcards: svpoplib.sampleset.cluster_param_cpu(wildcards, config),
        mem=lambda wildcards: svpoplib.sampleset.cluster_param_mem(wildcards, config),
        rt=lambda wildcards: svpoplib.sampleset.cluster_param_rt(wildcards, config)
    run:

        # Get sample and input list
        sampleset_entry = svpoplib.sampleset.get_config_entry(wildcards.sourcename, wildcards.sample, config)
        merge_strategy = svpoplib.sampleset.get_merge_strategy(sampleset_entry, wildcards.vartype, wildcards.svtype)

        # Merge
        df = svpoplib.svmerge.merge_variants(
            input.bed,
            sampleset_entry['samples'],
            merge_strategy['strategy'],
            subset_chrom=wildcards.chrom,
            threads=params.cpu
        )

        # Bylen to byref
        if df.shape[0] > 0:
            df['END'] = df.apply(lambda row: (row['POS'] + 1) if row['SVTYPE'] == 'INS' else row['END'], axis=1)

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# variant_sampleset_fa_merge
#
# Merge FASTA files for sampleset.
rule variant_sampleset_fa_merge:
    input:
        bed='results/variant/sampleset/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        fa=lambda wildcards: svpoplib.sampleset.get_sample_set_input(
            wildcards.sourcename,
            wildcards.sample,
            'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz',
            config,
            wildcards
        )
    output:
        fa='results/variant/sampleset/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz'
    wildcard_constraints:
        svtype='ins|del|inv|sub|rgn'
    run:

        fa_input_pattern = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz'

        df = pd.read_csv(input.bed, sep='\t', usecols=('ID', 'MERGE_SAMPLES', 'MERGE_VARIANTS'))

        df['SAMPLE'] = df['MERGE_SAMPLES'].apply(lambda val: val.split(',')[0])
        df['ID_SAMPLE'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[0])

        del df['MERGE_SAMPLES']
        del df['MERGE_VARIANTS']

        sampleset_entry = svpoplib.sampleset.get_config_entry(wildcards.sourcename, wildcards.sample, config)

        with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
            SeqIO.write(
                svpoplib.sampleset.fa_write_func(
                    df, wildcards, sampleset_entry, fa_input_pattern, config
                ),
                out_file, 'fasta'
            )
