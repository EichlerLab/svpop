"""
Annotate variants by homopolymers.
"""

# variant_anno_homopoly_intersect
#
# Intersect with homopolymers.
# * all: All homopolymer regions (could be a minimum size at some time).
# * any: Any intersect
rule variant_anno_homopoly_intersect:
    input:
        var_bed='results/variant/{sourcetype}/{caller}/bed4/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        hpy_bed='data/anno/homopolymer/{rep_len}_regions_all.bed'
    output:
        tab='results/variant/{sourcetype}/{caller}/anno/{sample}/all/{filter}/homopolymer/{rep_len}_intersect_all_any_{vartype}_{svtype}.tab'
    wildcard_constraints:
        sourcetype='caller|varset',
        svtype='ins|del|inv|dup|snv|rgn|sub'
    shell:
        """echo -e "ID\\tBASE\\tLEN" > {output.tab}; """
        """bedtools intersect -a {input.var_bed} -b {input.hpy_bed} -loj -sorted | """
        """awk -vOFS="\\t" '($5 != ".") {{print $4, $8, $9}}' """
        """>> {output.tab}"""


# variant_anno_homopolymer_nearest
#
# Find nearest non-intersecting homopolymer.
rule variant_anno_homopolymer_nearest:
    input:
        var_bed='results/variant/{sourcetype}/{caller}/bed4/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        hpy_bed='data/anno/homopolymer/{rep_len}_regions_all.bed'
    output:
        tab='results/variant/{sourcetype}/{caller}/anno/{sample}/all/{filter}/homopolymer/{rep_len}_nearest_all_{flank}_{vartype}_{svtype}.tab'
    wildcard_constraints:
        sourcetype='caller|varset',
        svtype='ins|del|inv|dup|snv|rgn|sub',
        flank='up|dn'
    run:

        # Get ignored direction
        if wildcards.flank == 'dn':
            ignore_dir = 'u'
        elif wildcards.flank == 'up':
            ignore_dir = 'd'
        else:
            raise RuntimeError('Unknown flank: {}'.format(wildcards.flank))

        # Intersect
        shell(
            """echo -e "ID\\tBASE\\tLEN\\tDISTANCE" > {output.tab}; """
            """bedtools closest -io -D ref -i{ignore_dir} -a {input.var_bed} -b {input.hpy_bed} -sorted | """
            """awk -vOFS="\\t" '"""
                """function abs(v) {{return v < 0 ? -v : v}} """
                """($4 != ".") {{print $4, $8, $9, abs($10)}}"""
            """' """
            """>> {output.tab}"""
        )

        # abs() from: https://unix.stackexchange.com/questions/220588/how-to-take-the-absolute-value-using-awk
