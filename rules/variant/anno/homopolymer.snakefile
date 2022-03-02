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
        var_bed='temp/variant/caller/{sourcename}/{sample}/{filter}/all/bed4/{vartype}_{svtype}.bed.gz',
        hpy_bed='data/anno/homopolymer/{rep_len}_regions_all.bed.gz'
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/homopolymer/{rep_len}_intersect_all_any_{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    shell:
        """{{\n"""
        """    echo -e "ID\\tBASE\\tLEN";\n"""
        """    bedtools intersect -a {input.var_bed} -b {input.hpy_bed} -loj -sorted | """
            """awk -vOFS="\\t" '($5 != ".") {{print $4, $8, $9}}';\n"""
        """}} | gzip """
        """> {output.tsv}"""


# variant_anno_homopolymer_nearest
#
# Find nearest non-intersecting homopolymer.
rule variant_anno_homopolymer_nearest:
    input:
        var_bed='temp/variant/caller/{sourcename}/{sample}/{filter}/all/bed4/{vartype}_{svtype}.bed.gz',
        hpy_bed='data/anno/homopolymer/{rep_len}_regions_all.bed.gz'
    output:
        tab='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/homopolymer/{rep_len}_nearest_all_{flank}_{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
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
            """{{\n"""
            """    echo -e "ID\\tBASE\\tLEN\\tDISTANCE"; """
            """    bedtools closest -io -D ref -i{ignore_dir} -a {input.var_bed} -b {input.hpy_bed} -sorted | """
            """    awk -vOFS="\\t" '"""
                """function abs(v) {{return v < 0 ? -v : v}} """
                """($4 != ".") {{print $4, $8, $9, abs($10)}}"""
                """';\n """
            """}} | gzip """
            """> {output.tab}"""
        )

        # abs() from: https://unix.stackexchange.com/questions/220588/how-to-take-the-absolute-value-using-awk
