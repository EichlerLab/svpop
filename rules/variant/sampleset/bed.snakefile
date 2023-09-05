"""
Merge variants from multiple samples.
"""

#############
### Rules ###
#############

def _variant_sampleset_bed_get_partition_count(wildcards):
        sampleset_entry = svpoplib.sampleset.get_config_entry(wildcards.sourcename, None, config)

        partition_count = sampleset_entry.get('params', dict()).get('partitions', svpoplib.constants.DEFAULT_PARTITIONS)

        try:
            return int(partition_count)
        except ValueError:
            raise RuntimeError(f'Error getting "config/partitions" from sampleset config for {wildcards.sourcename}: Partition count is not an integer: {partition_count}')

# variant_sampleset_bed_merge_svindel
#
# Merge callerset variants from multiple sources for one sample.
rule variant_sampleset_bed_merge_svindel:
    input:
        bed=lambda wildcards: expand(
            'temp/variant/sampleset/{{sourcename}}/{{sample}}/{{filter}}/all/bed/svindel_{{svtype}}/part_{part}.bed.gz',
            part=range(_variant_sampleset_bed_get_partition_count(wildcards))
        )
    output:
        bed_sv='results/variant/sampleset/{sourcename}/{sample}/{filter}/all/bed/sv_{svtype}.bed.gz',
        bed_indel='results/variant/sampleset/{sourcename}/{sample}/{filter}/all/bed/indel_{svtype}.bed.gz'
    wildcard_constraints:
        svtype='ins|del'
    params:
        mem=lambda wildcards: svpoplib.sampleset.cluster_param_anno_mem(wildcards, config, 'svindel')
    run:

        df = pd.concat(
            [pd.read_csv(bed_file_name, sep='\t') for bed_file_name in input.bed if os.stat(bed_file_name).st_size > 0],
            axis=0,
            sort=False
        ).sort_values(
            ['#CHROM', 'POS']
        )

        # Write
        df.loc[df['SVLEN'] >= 50].to_csv(output.bed_sv, sep='\t', index=False, compression='gzip')
        df.loc[df['SVLEN'] < 50].to_csv(output.bed_indel, sep='\t', index=False, compression='gzip')

# variant_sampleset_bed_merge
#
# Merge callerset variants from multiple sources for one sample.
rule variant_sampleset_bed_merge:
    input:
        bed=lambda wildcards: expand(
            'temp/variant/sampleset/{{sourcename}}/{{sample}}/{{filter}}/all/bed/{{varsvtype}}/part_{part}.bed.gz',
            part=range(_variant_sampleset_bed_get_partition_count(wildcards))
        )
    output:
        bed='results/variant/sampleset/{sourcename}/{sample}/{filter}/all/bed/{varsvtype}.bed.gz'
    wildcard_constraints:
        varsvtype='sv_inv|sv_dup|snv_snv'
    params:
        mem=lambda wildcards: svpoplib.sampleset.cluster_param_anno_mem(wildcards, config, wildcards.varsvtype.split('_')[0], wildcards.varsvtype.split('_')[1])
    run:

        df = pd.concat(
            [pd.read_csv(bed_file_name, sep='\t') for bed_file_name in input.bed if os.stat(bed_file_name).st_size > 0],
            axis=0,
            sort=False
        ).sort_values(
            ['#CHROM', 'POS']
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# variant_sampleset_bed_merge_part
#
# Get a set of variants merged from multiple samples.
rule variant_sampleset_bed_merge_part:
    input:
        bed=lambda wildcards: svpoplib.sampleset.get_sample_set_input(
            wildcards.sourcename,
            wildcards.sample,
            'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
            config,
            wildcards
        ),
        fa=lambda wildcards: svpoplib.sampleset.get_sample_set_input(
            wildcards.sourcename,
            wildcards.sample,
            'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz',
            config,
            wildcards
        ) if svpoplib.sampleset.is_read_seq(wildcards, config) else [],
        tsv_part='results/variant/sampleset/{sourcename}/partition_table.tsv.gz'

    output:
        bed=temp(
            'temp/variant/sampleset/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}/part_{part}.bed.gz'
        )
    params:
        mem=lambda wildcards: svpoplib.sampleset.cluster_param_mem(wildcards, config),
        rt=lambda wildcards: svpoplib.sampleset.cluster_param_rt(wildcards, config)
    threads: lambda wildcards: svpoplib.sampleset.cluster_param_cpu(wildcards, config)
    run:

        # Get sample and input list
        sampleset_entry = svpoplib.sampleset.get_config_entry(wildcards.sourcename, wildcards.sample, config)
        merge_strategy = svpoplib.sampleset.get_merge_strategy(sampleset_entry, wildcards.vartype, wildcards.svtype, config)

        # Get chromosome list
        df_part = pd.read_csv(input.tsv_part, sep='\t')
        df_part = df_part.loc[df_part['PARTITION'] == int(wildcards.part)]

        if df_part.shape[0] == 0:
            with open(output.bed, 'wt') as out_file:
                pass  # Empty output file if no chromosomes in this partition

        # Merge
        df = svpoplib.svmerge.merge_variants(
            input.bed,
            sampleset_entry['samples'],
            merge_strategy['strategy'],
            fa_list=input.fa if input.fa else None,
            subset_chrom=set(df_part['CHROM']),
            threads=threads
        )

        # Bylen to byref
        if df.shape[0] > 0:
            df['END'] = df.apply(lambda row: (row['POS'] + 1) if row['SVTYPE'] == 'INS' else row['END'], axis=1)

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# variant_sampleset_batch_table
#
# Assign each chromosome to a batch.
rule variant_sampleset_batch_table:
    output:
        tsv_part='results/variant/sampleset/{sourcename}/partition_table.tsv.gz'
    run:

        partition_count = _variant_sampleset_bed_get_partition_count(wildcards)

        # Read (Series of chrom -> len)
        df = svpoplib.ref.get_df_fai(config['reference_fai'])

        # Get a list of partitions (each element is a tuple of chrom names)
        partitions = svpoplib.partition.chrom.partition(df, partition_count)

        # Build table
        df = pd.DataFrame(df)

        df['PARTITION'] = -1

        for index in range(len(partitions)):
            for chrom in partitions[index]:
                if df.loc[chrom, 'PARTITION'] != -1:
                    raise RuntimeError(f'Chromosome {chrom} assigned to multiple partitions')
                df.loc[chrom, 'PARTITION'] = index

        missing_list = list(df.loc[df['PARTITION'] == -1].index)

        if missing_list:
            raise RuntimeError(f'Failed assigning {len(missing_list)} chromosomes to partitions: {", ".join(missing_list)}')

        # Write
        df.to_csv(output.tsv_part, sep='\t', index=True, compression='gzip')


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
        fa='results/variant/sampleset/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz',
        fai='results/variant/sampleset/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz.fai'
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

        shell("""samtools faidx {output.fa}""")
