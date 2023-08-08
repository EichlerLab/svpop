"""
Variant set intersections.
"""


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

        df['SOURCE_SET'] = df['SOURCE_SET'].apply(lambda val: set(val.split(',')) - {''})

        # Read original variant table
        len_a = pd.read_csv(input.a, sep='\t', usecols=('ID', 'SVLEN'), index_col='ID').squeeze()
        len_b = pd.read_csv(input.b, sep='\t', usecols=('ID', 'SVLEN'), index_col='ID').squeeze()

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

ruleorder: var_intersect_by_merge_svindel > var_intersect_by_merge

# var_intersect_by_merge_svindel
#
# Overlap variants: svindel
rule var_intersect_by_merge_svindel:
    input:
        bed=lambda wildcards: svpoplib.intersect.intersect_get_input(
            [
                'results/variant/{sourcetype_a}/{sourcename_a}/{sample_a}/{filter_a}/{svset_a}/bed/svindel_{svtype}.bed.gz',
                'results/variant/{sourcetype_b}/{sourcename_b}/{sample_b}/{filter_b}/{svset_b}/bed/svindel_{svtype}.bed.gz'
            ], wildcards, config
        ),
        fa=lambda wildcards: svpoplib.intersect.intersect_get_input(
            [
                'results/variant/{sourcetype_a}/{sourcename_a}/{sample_a}/{filter_a}/{svset_a}/bed/fa/svindel_{svtype}.fa.gz',
                'results/variant/{sourcetype_b}/{sourcename_b}/{sample_b}/{filter_b}/{svset_b}/bed/fa/svindel_{svtype}.fa.gz'
            ], wildcards, config
        )
    output:
        tsv='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/svindel_{svtype}/intersect.tsv.gz',
        tsv_sv='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/sv_{svtype}/intersect.tsv.gz',
        tsv_indel='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/indel_{svtype}/intersect.tsv.gz'
    threads: 8
    wildcard_constraints:
        svtype='ins|del'
    run:


        # Intersect
        df = svpoplib.intersect.run_intersect(
            input.bed,
            svpoplib.svmerge.get_merge_def(wildcards.merge_def, config),
            input.fa if len(input.fa) > 0 else None,
            threads=threads
        )

        # Read SV set
        df_a = pd.read_csv(input.bed[0], sep='\t', usecols=('ID', 'SVLEN'))
        df_b = pd.read_csv(input.bed[1], sep='\t', usecols=('ID', 'SVLEN'))

        sv_set_a = set(df_a.loc[df_a['SVLEN'] >= 50, 'ID'])
        sv_set_b = set(df_b.loc[df_b['SVLEN'] >= 50, 'ID'])

        del(df_a)
        del(df_b)

        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

        df.loc[
            df['ID_A'].isin(sv_set_a) | df['ID_B'].isin(sv_set_a)
        ].to_csv(
            output.tsv_sv, sep='\t', index=False, compression='gzip'
        )

        df.loc[
            (~ df['ID_A'].isin(sv_set_a)) | (~ df['ID_B'].isin(sv_set_a))
        ].to_csv(
            output.tsv_indel, sep='\t', index=False, compression='gzip'
        )

# var_intersect_by_merge
#
# Overlap variants: non-svindel.
# svindel: Rule "var_intersect_by_merge_svindel" has higher priority.
rule var_intersect_by_merge:
    input:
        bed=lambda wildcards: svpoplib.intersect.intersect_get_input(
            [
                'results/variant/{sourcetype_a}/{sourcename_a}/{sample_a}/{filter_a}/{svset_a}/bed/{vartype}_{svtype}.bed.gz',
                'results/variant/{sourcetype_b}/{sourcename_b}/{sample_b}/{filter_b}/{svset_b}/bed/{vartype}_{svtype}.bed.gz'
            ], wildcards, config
        ),
        fa=lambda wildcards: svpoplib.intersect.intersect_get_input(
            [
                'results/variant/{sourcetype_a}/{sourcename_a}/{sample_a}/{filter_a}/{svset_a}/bed/fa/{vartype}_{svtype}.fa.gz',
                'results/variant/{sourcetype_b}/{sourcename_b}/{sample_b}/{filter_b}/{svset_b}/bed/fa/{vartype}_{svtype}.fa.gz'
            ], wildcards, config
        )
    output:
        tsv='results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/intersect.tsv.gz'
    threads: 8
    wildcard_constraints:
        vartype='sv|indel|snv|sub|rgn',
        svtype='ins|del|inv|dup|sub|snv|rgn'
    run:


        # Intersect
        df = svpoplib.intersect.run_intersect(
            input.bed,
            svpoplib.svmerge.get_merge_def(wildcards.merge_def, config),
            input.fa if len(input.fa) > 0 else None,
            threads=threads
        )

        # Read SV set
        df_a = pd.read_csv(input.bed[0], sep='\t', usecols=('ID', 'SVLEN'))
        df_b = pd.read_csv(input.bed[1], sep='\t', usecols=('ID', 'SVLEN'))

        sv_set_a = set(df_a.loc[df_a['SVLEN'] >= 50, 'ID'])
        sv_set_b = set(df_b.loc[df_b['SVLEN'] >= 50, 'ID'])

        del(df_a)
        del(df_b)

        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

        df.to_csv(
            output.tsv, sep='\t', index=False, compression='gzip'
        )
