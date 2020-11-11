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
        tab_ins='results/variant/intersect_nearest/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{vartype}_ins/intersect.tab',
        tab_del='results/variant/intersect_nearest/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{vartype}_del/intersect.tab'
    output:
        tab='results/variant/intersect_nearest/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{vartype}_insdel/intersect.tab'
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

# intersect_nearest_distance_table
#
# Get table of distance between variants in two call sets.
rule intersect_nearest_distance_table:
    input:
        a='results/variant/{sourcetype_a}/{sourcename_a}/bed/{sample_a}/{svset}/{filter}/byref/{vartype}_{svtype}.bed',
        b='results/variant/{sourcetype_b}/{sourcename_b}/bed/{sample_b}/{svset}/{filter}/byref/{vartype}_{svtype}.bed'
    output:
        tab='results/variant/intersect_nearest/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{vartype}_{svtype}/intersect.tab'
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|rgn|snv'
    run:

        # Match variants
        df = analib.variant.var_nearest(
            pd.read_csv(input.a, sep='\t', header=0),
            pd.read_csv(input.b, sep='\t', header=0),
            ref_alt=False
        )

        # Write
        df.to_csv(output.tab, sep='\t', index=False)


### Nearest - Match Ref/Alt ###

# intersect_nearest_ra_distance_insdel
#
# Merge INS/DEL.
rule intersect_nearest_ra_distance_insdel:
    input:
        tab_ins='results/variant/intersect_nearest_ra/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{vartype}_ins/intersect.tab',
        tab_del='results/variant/intersect_nearest_ra/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{vartype}_del/intersect.tab'
    output:
        tab='results/variant/intersect_nearest_ra/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{vartype}_insdel/intersect.tab'
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

# intersect_nearest_distance_table
#
# Get table of distance between variants in two call sets. Variants are matched on REF/ALT columns (must be the same to
# be compared).
rule intersect_nearest_ra_distance_table:
    input:
        a='results/variant/{sourcetype_a}/{sourcename_a}/bed/{sample_a}/{svset}/{filter}/byref/{vartype}_{svtype}.bed',
        b='results/variant/{sourcetype_b}/{sourcename_b}/bed/{sample_b}/{svset}/{filter}/byref/{vartype}_{svtype}.bed'
    output:
        tab='results/variant/intersect_nearest_ra/{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}/{svset}/{filter}/{vartype}_{svtype}/intersect.tab'
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|rgn|snv'
    run:

        # Match variants
        df = analib.variant.var_nearest(
            pd.read_csv(input.a, sep='\t', header=0),
            pd.read_csv(input.b, sep='\t', header=0),
            ref_alt=True
        )

        # Write
        df.to_csv(output.tab, sep='\t', index=False)
