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
        tab='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{merge_def}/{vartype}_{svtype}/intersect.tab'
    output:
        pdf='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{merge_def}/{vartype}_{svtype}/venn/variant_venn.pdf',
        png='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{merge_def}/{vartype}_{svtype}/venn/variant_venn.png'
    run:

        # Read intersect table

        df = pd.read_csv(input.tab, sep='\t')

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
        analib.plot.venn.make_venn(
            set_a, set_b, [output.pdf, output.png],
            name_a, name_b,
            '{} - {}'.format(svtype_label, svset_label),
            len_stat = wildcards.svtype != 'snv'
        )


#
# Merge "all" and "insdel" sets.
#

# var_intersect_ro_insdel
#
# Merge INS/DEL.
rule var_intersect_ro_insdel:
    input:
        tab_ins='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{merge_def}/{vartype}_ins/intersect.tab',
        tab_del='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{merge_def}/{vartype}_del/intersect.tab',
    output:
        tab='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{merge_def}/{vartype}_insdel/intersect.tab'
    run:

        # Read and merge
        df = pd.concat(
            [
                pd.read_csv(input.tab_ins, sep='\t', header=0),
                pd.read_csv(input.tab_del, sep='\t', header=0)
            ]
        )

        # Write
        df.to_csv(output.tab, sep='\t', index=False)

# var_intersect_ro_all
#
# Merge INS/DEL/INV.
rule var_intersect_ro_all:
    input:
        tab_ins='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{merge_def}/{vartype}_ins/intersect.tab',
        tab_del='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{merge_def}/{vartype}_del/intersect.tab',
        tab_inv='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{merge_def}/{vartype}_inv/intersect.tab'
    output:
        tab='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{merge_def}/{vartype}_all/intersect.tab'
    run:

        # Read and merge
        df = pd.concat(
            [
                pd.read_csv(input.tab_ins, sep='\t', header=0),
                pd.read_csv(input.tab_del, sep='\t', header=0),
                pd.read_csv(input.tab_inv, sep='\t', header=0)
            ]
        )

        # Write
        df.to_csv(output.tab, sep='\t', index=False)

#
# Do intersect
#

# var_intersect_by_merge
#
# Get sets of variants with reciprocal overlaps using the same svset filter.
rule var_intersect_by_merge:
    input:
        a='results/variant/{sourcetype_a}/{sourcename_a}/bed/{sample_a}/{svset}/{filter}/byref/{vartype}_{svtype}.bed',
        b='results/variant/{sourcetype_b}/{sourcename_b}/bed/{sample_b}/{svset}/{filter}/byref/{vartype}_{svtype}.bed'
    output:
        tab='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{merge_def}/{vartype}_{svtype}/intersect.tab'
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|snv'
#        svset='^[^(_vs_)]+$'
    params:
        cpu=8,
        mem='8G',
        rt='96:00:00'
    run:

        # Get configured merge definition
        config_def = analib.svmerge.get_merge_def(wildcards.merge_def, config)

        # Merge
        df = analib.svmerge.merge_variants(
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

        df.to_csv(output.tab, sep='\t', index=False)


# var_intersect_bymerge_svset_diff
#
# Get sets of variants with reciprocal overlaps using a different svset filter for each sample.
rule var_intersect_bymerge_svset_diff:
    input:
        a='results/variant/{sourcetype_a}/{sourcename_a}/bed/{sample_a}/{svset_a}/{filter}/byref/{vartype}_{svtype}.bed',
        b='results/variant/{sourcetype_b}/{sourcename_b}/bed/{sample_b}/{svset_b}/{filter}/byref/{vartype}_{svtype}.bed'
    output:
        tab='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset_a}_vs_{svset_b}/{filter}/{merge_def}/{vartype}_{svtype}/intersect.tab'
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|snv'
    run:

        # Get configured merge definition
        config_def = analib.svmerge.get_merge_def(wildcards.merge_def, config)

        if config_def is None:
            config_def = wildcards.merge_def

        # Merge
        df = analib.svmerge.merge_variants(
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

        df.to_csv(output.tab, sep='\t', index=False)
