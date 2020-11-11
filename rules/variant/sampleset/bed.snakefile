"""
Merge variants from multiple samples.
"""


###################
### Definitions ###
###################




#############
### Rules ###
#############

# variant_sampleset_bed_merge
#
# Merge callerset variants from multiple sources for one sample.
rule variant_sampleset_bed_merge:
    input:
        bed=expand(
            'temp/variant/sampleset/{{sampleset}}/bed/{{samplelist}}/all/{{filter}}/byref/{{vartype}}_{{svtype}}/chrom_{chrom}.bed',
            chrom=analib.ref.get_df_fai(config['reference_fai']).index
        )
    output:
        bed='results/variant/sampleset/{sampleset}/bed/{samplelist}/all/{filter}/byref/{vartype}_{svtype}.bed'
    wildcard_constraints:
        vartype='sv|indel|snv',
        svtype='ins|del|inv|dup|snv',
        samplelist='[^/]+'
    params:
        mem=lambda wildcards: analib.sampleset.cluster_param_anno_mem(wildcards, config),
    run:

        df = pd.concat(
            [pd.read_csv(bed_file_name, sep='\t') for bed_file_name in input.bed],
            axis=0,
            sort=False
        ).sort_values(
            ['#CHROM', 'POS']
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False)


# variant_sampleset_bed_merge_chrom
#
# Get a set of variants merged from multiple samples.
rule variant_sampleset_bed_merge_chrom:
    input:
        bed=lambda wildcards: analib.sampleset.get_sample_set_input(
            wildcards.sampleset,
            wildcards.samplelist,
            'results/variant/{sourcetype}/{sourcename}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
            config,
            wildcards
        )
    output:
        bed=temp(
            'temp/variant/sampleset/{sampleset}/bed/{samplelist}/all/{filter}/byref/{vartype}_{svtype}/chrom_{chrom}.bed'
        )
    params:
        cpu=lambda wildcards: analib.sampleset.cluster_param_cpu(wildcards, config),
        mem=lambda wildcards: analib.sampleset.cluster_param_mem(wildcards, config),
        rt=lambda wildcards: analib.sampleset.cluster_param_rt(wildcards, config)
    run:

        # Get sample and input list
        sampleset_entry = analib.sampleset.get_config_entry(wildcards.sampleset, wildcards.samplelist, config)
        merge_strategy = analib.sampleset.get_merge_strategy(sampleset_entry, wildcards.vartype, wildcards.svtype)

        # Merge
        df = analib.svmerge.merge_variants(
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
        df.to_csv(output.bed, sep='\t', index=False)

# variant_sampleset_fa_merge
#
# Merge FASTA files for sampleset.
rule variant_sampleset_fa_merge:
    input:
        bed='results/variant/sampleset/{sampleset}/bed/{samplelist}/all/{filter}/byref/{vartype}_{svtype}.bed',
        fa=lambda wildcards: analib.sampleset.get_sample_set_input(
            wildcards.sampleset,
            wildcards.samplelist,
            'results/variant/{sourcetype}/{sourcename}/fasta/{sample}/all/{filter}/{vartype}_{svtype}.fa.gz',
            config,
            wildcards
        )
    output:
        fa='results/variant/sampleset/{sampleset}/fasta/{samplelist}/all/{filter}/{vartype}_{svtype}.fa.gz'
    wildcard_constraints:
        svtype='ins|del|inv|sub|rgn'
    run:

        df = pd.read_csv(input.bed, sep='\t', usecols=('ID', 'MERGE_SAMPLES', 'MERGE_VARIANTS'))

        df['SAMPLE'] = df['MERGE_SAMPLES'].apply(lambda val: val.split(',')[0])
        df['ID_SAMPLE'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[0])

        del df['MERGE_SAMPLES']
        del df['MERGE_VARIANTS']

        sampleset_entry = analib.sampleset.get_config_entry(wildcards.sampleset, wildcards.samplelist, config)

        with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
            SeqIO.write(fa_write_func(df, wildcards, sampleset_entry), out_file, 'fasta')

def fa_write_func(df, wildcards, sampleset_entry):
    """
    Function to yield a sequence record iterator for rule variant_sampleset_fa_merge.

    :param df: DataFrame with SAMPLE, ID_SAMPLE (original pre-merged ID), and ID.
    :param wildcards: Rule wildcards.
    :param sampleset_entry: Sampleset entry.

    :return: SeqRecord iterator.
    """

    fa_dict = {sample: fa_file for sample, fa_file in analib.sampleset.get_sample_set_input(
        wildcards.sampleset,
        wildcards.samplelist,
        'results/variant/{sourcetype}/{sourcename}/fasta/{sample}/all/{filter}/{vartype}_{svtype}.fa.gz',
        config,
        wildcards,
        as_tuple=True
    )}

    # Process samples
    for sample in sampleset_entry['samples']:

        id_dict = dict(df.loc[df['SAMPLE'] == sample, ['ID_SAMPLE', 'ID']].set_index('ID_SAMPLE').squeeze())
        id_set = set(id_dict.keys())

        # Open file
        in_file_name = fa_dict[sample]

        found_ids = set()

        with gzip.open(in_file_name, 'rt') as in_file:
            for record in SeqIO.parse(in_file, 'fasta'):

                found_ids.add(record.id)

                if record.id in id_set:
                    record.id = id_dict[record.id]
                    yield record

        # Check for missing IDs
        missing_ids = id_set - found_ids

        if missing_ids:
            raise RuntimeError('Missing {} variant ID(s) for sample {} when merging FASTAs: {}{}: {}'.format(
                len(missing_ids), sample,
                ', '.join(sorted(missing_ids)[:3]), '...' if len(missing_ids) > 3 else '',
                in_file_name
            ))
