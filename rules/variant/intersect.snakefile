"""
Variant set intersections.
"""

ruleorder: var_intersect_bymerge_svset_diff > var_intersect_by_merge

###################
### Definitions ###
###################

def intersect_is_read_seq(wildcards, config):
    """
    Determine if merge requires input sequence.

    :param wildcards: Rule wildcards.
    :param config: Configuration.

    :return: `True` if sequences should be read.
    """

    config_def = svpoplib.svmerge.get_merge_def(wildcards.merge_def, config)

    if config_def is None:
        config_def = wildcards.merge_def

    return svpoplib.svmergeconfig.params.get_merge_config(config_def).read_seq

MERGE_INFO_FIELD_LIST = [
    'MERGE_OFFSET', 'MERGE_RO', 'MERGE_SZRO', 'MERGE_OFFSZ', 'MERGE_MATCH'
]


#############
### Rules ###
#############

#
# Variant comparison Venn by sample
#

# var_intersect_venn
#
# Make Venn PDFs.
rule var_intersect_venn:
    input:
        a='results/variant/{sourcetype_a}/{sourcename_a}/{sample_a}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz',
        b='results/variant/{sourcetype_b}/{sourcename_b}/{sample_b}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz',
        tsv='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/intersect.tsv.gz'
    output:
        img='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/venn/variant_venn.{ext}'
    run:

        # Matplotlib imports
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        # Read intersect table
        df = pd.read_csv(input.tsv, sep='\t')

        df['SOURCE_SET'] = df['SOURCE_SET'].apply(lambda val: set(val.split(',')))

        # Read original variant table
        len_a = pd.read_csv(input.a, sep='\t', usecols=('ID', 'SVLEN'), index_col='ID', squeeze=True)
        len_b = pd.read_csv(input.b, sep='\t', usecols=('ID', 'SVLEN'), index_col='ID', squeeze=True)

        df['ID'] = df.apply(lambda row:
            'A:{}'.format(row['ID_A']) if 'A' in row['SOURCE_SET'] else
                'B:{}'.format(row['ID_B']),
            axis=1
        )

        # Get a series of lengths per variant (average if intersecting two different size variants)
        if wildcards.svtype != 'snv':
            def val_len(row):
                if {'A', 'B'} == row['SOURCE_SET']:
                    return int(
                        (len_a[row['ID_A']] + len_b[row['ID_B']]) / 2
                    )

                if {'A'} == row['SOURCE_SET']:
                    return len_a[row['ID_A']]

                if {'B'} == row['SOURCE_SET']:
                    return len_b[row['ID_B']]

                raise RuntimeError('Cannot get variant lengths: SOURCE_SET is not A, B, or A and B: "{}"'.format(','.join(row['SOURCE_SET'])))

            df['LEN'] = df.apply(val_len, axis=1)

            len_series = df.set_index('ID')['LEN']
        else:
            len_series = None

        # Get sets of variants.
        # For set A, choose names from A
        # For set B, choose names from A for variants that intersected, and choose names from B for variants that did
        # not intersect.

        set_a = set(df.loc[df['SOURCE_SET'].apply(lambda source_set: 'A' in source_set), 'ID'])
        set_b = set(df.loc[df['SOURCE_SET'].apply(lambda source_set: 'B' in source_set), 'ID'])

        # Get repeat class label
        svset_label = get_svset_label(wildcards.svset)

        # Get variant type label
        svtype_label = get_svtype_label(wildcards.svtype)

        # Set sample name
        name_a = get_sample_name(wildcards.sourcetype_a, wildcards.sourcename_a, wildcards.sample_a, wildcards.sourcename_a != wildcards.sourcename_b)
        name_b = get_sample_name(wildcards.sourcetype_b, wildcards.sourcename_b, wildcards.sample_b, wildcards.sourcename_a != wildcards.sourcename_b)

        # Make Venn
        fig = svpoplib.plot.venn.get_venn_fig(
            set_a, set_b,
            name_a, name_b,
            '{} - {}'.format(svtype_label, svset_label),
            len_series
        )

        # Write
        fig.savefig(output.img, bbox_inches='tight')
        plt.close(fig)


#
# Merge combined svtypes.
#

# var_intersect_combined_insdel
#
# Merge INS/DEL.
rule var_intersect_combined_insdel:
    input:
        tsv_ins='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_ins/intersect.tsv.gz',
        tsv_del='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_del/intersect.tsv.gz',
    output:
        tsv='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_insdel/intersect.tsv.gz'
    wildcard_constraints:
        vartype='sv|indel'
    run:

        # Read and merge
        df = pd.concat(
            [
                pd.read_csv(input.tsv_ins, sep='\t', header=0),
                pd.read_csv(input.tsv_del, sep='\t', header=0)
            ]
        )

        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

# var_intersect_combined_insdelinv
#
# Merge INS/DEL/INV.
rule var_intersect_combined_insdelinv:
    input:
        tsv_ins='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/sv_ins/intersect.tsv.gz',
        tsv_del='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/sv_del/intersect.tsv.gz',
        tsv_inv='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/sv_inv/intersect.tsv.gz'
    output:
        tsv='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/sv_insdelinv/intersect.tsv.gz'
    run:

        # Read and merge
        df = pd.concat(
            [
                pd.read_csv(input.tsv_ins, sep='\t', header=0),
                pd.read_csv(input.tsv_del, sep='\t', header=0),
                pd.read_csv(input.tsv_inv, sep='\t', header=0)
            ]
        )

        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

#
# Do intersect
#

# var_intersect_by_merge
#
# Get sets of variants with reciprocal overlaps using the same svset filter.
rule var_intersect_by_merge:
    input:
        a='results/variant/{sourcetype_a}/{sourcename_a}/{sample_a}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz',
        b='results/variant/{sourcetype_b}/{sourcename_b}/{sample_b}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz',
        fa=lambda wildcards:
            [
                'results/variant/{sourcetype_a}/{sourcename_a}/{sample_a}/{filter}/{svset}/bed/fa/{vartype}_{svtype}.fa.gz'.format(**wildcards),
                'results/variant/{sourcetype_b}/{sourcename_b}/{sample_b}/{filter}/{svset}/bed/fa/{vartype}_{svtype}.fa.gz'.format(**wildcards)
            ] if intersect_is_read_seq(wildcards, config) else []
    output:
        tsv='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/intersect.tsv.gz'
    params:
        cpu=8
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|snv|rgn'
    run:

        # Get configured merge definition
        config_def = svpoplib.svmerge.get_merge_def(wildcards.merge_def, config)

        # Merge
        df = svpoplib.svmerge.merge_variants(
            bed_list=[input.a, input.b],
            sample_names=['A', 'B'],
            strategy=config_def,
            threads=params.cpu,
            fa_list=input.fa if input.fa else None
        )

        support_col_list = [col for col in df.columns if col in MERGE_INFO_FIELD_LIST]

        # Subset columns
        df = df.loc[:, ['ID', 'MERGE_SAMPLES', 'MERGE_SRC', 'MERGE_VARIANTS'] + support_col_list]

        # Create an ID column for sample (empty string if the variant was not in that sample)
        df['MERGE_SAMPLES'] = df['MERGE_SAMPLES'].apply(lambda val:
            val + ',' if val == 'A' else (',' + val if val == 'B' else val)
        )

        df['MERGE_VARIANTS'] = df.apply(lambda row:
            row['MERGE_VARIANTS'] + ',' if row['MERGE_SAMPLES'] == 'A,' else (',' + row['MERGE_VARIANTS'] if row['MERGE_SAMPLES'] == ',B' else row['MERGE_VARIANTS']),
            axis=1
        )

        df['ID_A'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[0])
        df['ID_B'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[1])

        # Set support columns
        new_col_list = list()

        for col in support_col_list:
            new_col = col[len('MERGE_'):]
            new_col_list.append(new_col)

            split_list = df[col].apply(lambda val: val.split(',') if not pd.isnull(val) else '')

            df[new_col] = split_list.apply(lambda val: val[1] if len(val) > 1 else np.nan)

        # Subset and write
        df['SOURCE_SET'] = df['MERGE_SAMPLES']
        df = df.loc[:, ['ID_A', 'ID_B', 'SOURCE_SET'] + new_col_list]

        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')


# var_intersect_bymerge_svset_diff
#
# Get sets of variants with reciprocal overlaps using a different svset filter for each sample.
rule var_intersect_bymerge_svset_diff:
    input:
        a='results/variant/{sourcetype_a}/{sourcename_a}/{sample_a}/{filter}/{svset_a}/bed/{vartype}_{svtype}.bed.gz',
        b='results/variant/{sourcetype_b}/{sourcename_b}/{sample_b}/{filter}/{svset_b}/bed/{vartype}_{svtype}.bed.gz',
        fa=lambda wildcards:
            [
                'results/variant/{sourcetype_a}/{sourcename_a}/{sample_a}/{filter}/{svset}/bed/fa/{vartype}_{svtype}.fa.gz'.format(**wildcards),
                'results/variant/{sourcetype_b}/{sourcename_b}/{sample_b}/{filter}/{svset}/bed/fa/{vartype}_{svtype}.fa.gz'.format(**wildcards)
            ] if intersect_is_read_seq(wildcards, config) else []
    output:
        tsv='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset_a}_vs_{svset_b}/{vartype}_{svtype}/intersect.tsv.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|snv|rgn'
    run:

        # Get configured merge definition
        config_def = svpoplib.svmerge.get_merge_def(wildcards.merge_def, config)

        # Merge
        df = svpoplib.svmerge.merge_variants(
            bed_list=[input.a, input.b],
            sample_names=['A', 'B'],
            strategy=config_def,
            threads=params.cpu,
            fa_list=input.fa if input.fa else None
        )

        support_col_list = [col for col in df.columns if col in MERGE_INFO_FIELD_LIST]

        # Subset columns
        df = df.loc[:, ['ID', 'MERGE_SAMPLES', 'MERGE_SRC', 'MERGE_VARIANTS'] + support_col_list]

        # Create an ID column for sample (empty string if the variant was not in that sample)
        df.loc[df['MERGE_SAMPLES'] == 'A', 'MERGE_VARIANTS'] = df.loc[df['MERGE_SRC'] == 'A', 'MERGE_VARIANTS'] + ','
        df.loc[df['MERGE_SAMPLES'] == 'B', 'MERGE_VARIANTS'] = ',' + df.loc[df['MERGE_SRC'] == 'B', 'MERGE_VARIANTS']

        df['ID_A'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[0])
        df['ID_B'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[1])

        # Set support columns
        new_col_list = list()

        for col in support_col_list:
            new_col = col[len('MERGE_'):]
            new_col_list.append(new_col)

            split_list = df[col].apply(lambda val: val.split(',') if not pd.isnull(val) else '')

            df[new_col] = split_list.apply(lambda val: val[1] if len(val) > 1 else np.nan)

        # Subset and write
        df['SOURCE_SET'] = df['MERGE_SAMPLES']
        df = df.loc[:, ['ID_A', 'ID_B', 'SOURCE_SET'] + new_col_list]

        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

