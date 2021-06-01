"""
Variant set intersections.
"""

ruleorder: var_intersect_bymerge_svset_diff > var_intersect_by_merge

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
        tsv='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/intersect.tsv.gz'
    output:
        pdf='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/venn/variant_venn.pdf',
        png='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/venn/variant_venn.png'
    run:

        # Read intersect table

        df = pd.read_csv(input.tsv, sep='\t')

        # Get sets of variants.
        # For set A, choose names from A
        # For set B, choose names from A for variants that intersected, and choose names from B for variants that did
        # not intersect.

        set_a = set(df.loc[~ pd.isnull(df['ID_A']), 'ID_A'])
        set_b = set(df.loc[~ pd.isnull(df['ID_B'])].apply(lambda row: row['ID_A'] if not pd.isnull(row['ID_A']) else row['ID_B'], axis=1))

        # Get repeat class label
        svset_label = get_svset_label(wildcards.svset)

        # Get variant type label
        svtype_label = get_svtype_label(wildcards.svtype)

        # Set sample name
        name_a = get_sample_name(wildcards.sourcetype_a, wildcards.sourcename_a, wildcards.sample_a, wildcards.sourcename_a != wildcards.sourcename_b)
        name_b = get_sample_name(wildcards.sourcetype_b, wildcards.sourcename_b, wildcards.sample_b, wildcards.sourcename_a != wildcards.sourcename_b)

        # Set flag to get length distributions (mean & median)
        len_stat = wildcards.svtype != 'snv'

        # Make Venn
        svpoplib.plot.venn.make_venn(
            set_a, set_b, [output.pdf, output.png],
            name_a, name_b,
            '{} - {}'.format(svtype_label, svset_label),
            len_stat = wildcards.svtype != 'snv'
        )


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
        b='results/variant/{sourcetype_b}/{sourcename_b}/{sample_b}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz'
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
            threads=params.cpu
        )

        # Subset columns
        df = df.loc[:, ['ID', 'MERGE_SAMPLES', 'MERGE_SRC', 'MERGE_VARIANTS']]

        # Create an ID column for sample (empty string if the variant was not in that sample)
        df.loc[df['MERGE_SRC'] == 'A', 'MERGE_VARIANTS'] =  df.loc[df['MERGE_SRC'] == 'A', 'MERGE_VARIANTS'] + ','
        df.loc[df['MERGE_SRC'] == 'B', 'MERGE_VARIANTS'] =  ',' + df.loc[df['MERGE_SRC'] == 'B', 'MERGE_VARIANTS']

        df['ID_A'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[0])
        df['ID_B'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[1])

        # Subset and write
        df['SOURCE_SET'] = df['MERGE_SAMPLES']
        df = df.loc[:, ['ID_A', 'ID_B', 'SOURCE_SET']]

        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')


# var_intersect_bymerge_svset_diff
#
# Get sets of variants with reciprocal overlaps using a different svset filter for each sample.
rule var_intersect_bymerge_svset_diff:
    input:
        a='results/variant/{sourcetype_a}/{sourcename_a}/{sample_a}/{svset_a}/{filter}/bed/{vartype}_{svtype}.bed.gz',
        b='results/variant/{sourcetype_b}/{sourcename_b}/{sample_b}/{svset_b}/{filter}/bed/{vartype}_{svtype}.bed.gz'
    output:
        tsv='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset_a}_vs_{svset_b}/{vartype}_{svtype}/intersect.tsv.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|snv|rgn'
    run:

        # Get configured merge definition
        config_def = svpoplib.svmerge.get_merge_def(wildcards.merge_def, config)

        if config_def is None:
            config_def = wildcards.merge_def

        # Merge
        df = svpoplib.svmerge.merge_variants(
            bed_list=[input.a, input.b],
            sample_names=['A', 'B'],
            strategy=config_def,
            threads=6
        )

        # Subset columns
        df = df.loc[:, ['ID', 'MERGE_SAMPLES', 'MERGE_SRC', 'MERGE_VARIANTS']]

        # Create an ID column for sample (empty string if the variant was not in that sample)
        df.loc[df['MERGE_SRC'] == 'A', 'MERGE_VARIANTS'] =  df.loc[df['MERGE_SRC'] == 'A', 'MERGE_VARIANTS'] + ','
        df.loc[df['MERGE_SRC'] == 'B', 'MERGE_VARIANTS'] =  ',' + df.loc[df['MERGE_SRC'] == 'B', 'MERGE_VARIANTS']

        df['ID_A'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[0])
        df['ID_B'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[1])

        # Subset and write
        df['SOURCE_SET'] = df['MERGE_SAMPLES']
        df = df.loc[:, ['ID_A', 'ID_B', 'SOURCE_SET']]

        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')
