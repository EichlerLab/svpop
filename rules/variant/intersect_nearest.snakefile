"""
Find distance to nearest SV.
"""


#############
### Rules ###
#############

### Intersect ###

# intersect_nearest_distance_insdel
#
# Merge INS/DEL.
rule intersect_nearest_distance_insdel:
    input:
        tsv_ins='results/variant/intersect_nearest/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{filter}/{svset}/{vartype}_ins/intersect.tsv.gz',
        tsv_del='results/variant/intersect_nearest/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{filter}/{svset}/{vartype}_del/intersect.tsv.gz'
    output:
        tsv='results/variant/intersect_nearest/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{filter}/{svset}/{vartype}_insdel/intersect.tsv.gz'
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

# intersect_nearest_distance_table
#
# Get table of distance between variants in two call sets.
rule intersect_nearest_distance_table:
    input:
        a='results/variant/{sourcetype_a}/{sourcename_a}/{sample_a}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz',
        b='results/variant/{sourcetype_b}/{sourcename_b}/{sample_b}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz'
    output:
        tsv='results/variant/intersect_nearest/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{filter}/{svset}/{vartype}_{svtype}/intersect.tsv.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|rgn|snv'
    run:

        # Match variants
        df = svpoplib.variant.var_nearest(
            pd.read_csv(input.a, sep='\t', header=0),
            pd.read_csv(input.b, sep='\t', header=0),
            ref_alt=False
        )

        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')


### Nearest - Match Ref/Alt ###

# intersect_nearest_ra_distance_insdel
#
# Merge INS/DEL.
rule intersect_nearest_ra_distance_insdel:
    input:
        tsv_ins='results/variant/intersect_nearest_ra/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{filter}/{svset}/{vartype}_ins/intersect.tsv.gz',
        tsv_del='results/variant/intersect_nearest_ra/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{filter}/{svset}/{vartype}_del/intersect.tsv.gz'
    output:
        tsv='results/variant/intersect_nearest_ra/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{filter}/{svset}/{vartype}_insdel/intersect.tsv.gz'
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

# intersect_nearest_distance_table
#
# Get table of distance between variants in two call sets. Variants are matched on REF/ALT columns (must be the same to
# be compared).
rule intersect_nearest_ra_distance_table:
    input:
        a='results/variant/{sourcetype_a}/{sourcename_a}/{sample_a}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz',
        b='results/variant/{sourcetype_b}/{sourcename_b}/{sample_b}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz'
    output:
        tsv='results/variant/intersect_nearest_ra/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{filter}/{svset}/{vartype}_{svtype}/intersect.tsv.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|rgn|snv'
    run:

        # Match variants
        df = svpoplib.variant.var_nearest(
            pd.read_csv(input.a, sep='\t', header=0),
            pd.read_csv(input.b, sep='\t', header=0),
            ref_alt=True
        )

        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')
