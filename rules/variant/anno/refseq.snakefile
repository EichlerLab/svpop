"""
Intersect variants with RefSeq annotations.
"""

###################
### Definitions ###
###################

def _variant_anno_refseq_prox_get_distance(row, direction):
    """
    Get distance to the edge of a window. If the SV is upstream RefSeq annotation (relative to the RefSeq record,
    +/up or -/dn), then calculate the distance using the end-points of the SV and window upstream of the RefSeq record.
    Use start-points to calculate distance if the SV is downstream (-/up, +/dn).

    :param row: Row of the intersect table with "POS" and "END" representing the coordinates of the SV, and "txStart"
        and "txEnd" representing the window upstream or downstream of RefSeq annotations. Must also include "strand" to
        indicated if the RefSeq annotation orientation.
    :param direction: Direction of interest. Indicates whether windows are upstream (up) or downstream (dn) of the
        RefSeq record.
    """

    # Check arguments
    if direction not in ('up', 'dn'):
        raise RuntimeError('Unknown direction (up/dn): {}'.format(direction))

    if row['strand'] not in ('+', '-'):
        raise RuntimeError('Unknown strand (+/-) on line {}: {}'.format(row.name, row['strand']))

    # Get distance
    if (row['strand'] == '+' and direction == 'up') or (row['strand'] == '-' and direction == 'dn'):
        return row['txEnd'] - row['END'] + 1

    else:
        return row['POS'] - row['txStart'] + 1


#############
### Rules ###
#############

#
# RefSeq count - Count number of bases by RefSeq annotations (CDS, UTR5, etc)
#

# variant_anno_refseq_intersect_count
#
# Count affected bases for different gene regions (intron, CDS, etc) for each variant/gene intersection.
rule variant_anno_refseq_intersect_count:
    input:
        bed='temp/variant/caller/{sourcename}/{sample}/{filter}/all/anno/refseq/{vartype}_{svtype}_merged.bed.gz'
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/refseq/refseq-count_{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
        vartype='sv|indel|snv|rgn|sub',
        svtype='ins|del|inv|dup|snv|rgn|sub'
    run:

        df_sv = pd.read_csv(input.bed, sep='\t', header=0)

        df = svpoplib.refseq.get_refgene_counts_df(df_sv)

        df.to_csv(output.tsv, sep='\t', index=False)

# variant_anno_refseq_intersect_bed
#
# Merge variant BED with refgene annotations. Each variant is concatenated with each refseq gene it
# overlaps with.
rule variant_anno_refseq_intersect_bed:
    input:
        bed='temp/variant/caller/{sourcename}/{sample}/{filter}/all/anno/refseq/{vartype}_{svtype}_subcol.bed.gz',
        refseq='data/anno/refseq/refseq.bed.gz'
    output:
        bed=temp('temp/variant/caller/{sourcename}/{sample}/{filter}/all/anno/refseq/{vartype}_{svtype}_merged.bed.gz')
    shell:
        """{{\n"""
        """    {{\n"""
        """        head -n 1 <(zcat {input.bed});\n"""
        """        head -n 1 <(zcat {input.refseq});\n"""
        """    }} | \n"""
        """    awk -vORS="\\t" {{print}} | \n"""
        """    sed -re 's/\\s+/\\t/g' | \n"""
        """    sed -re 's/\\s*$/\\n/';\n"""
        """    bedtools intersect -a {input.bed} -b {input.refseq} -wa -wb;\n"""
        """}} | gzip > {output.bed}"""

# variant_anno_refseq_intersect_bed_subcol
#
# Subset BED columns for RefSeq intersections.
rule variant_anno_refseq_intersect_bed_subcol:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz'
    output:
        bed=temp('temp/variant/caller/{sourcename}/{sample}/{filter}/all/anno/refseq/{vartype}_{svtype}_subcol.bed.gz')
    run:

        cols = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN']

        df = pd.read_csv(input.bed, sep='\t', header=0, usecols=cols)

        df.loc[:, cols].to_csv(output.bed, sep='\t', index=False, compression='gzip')


#
# Refseq proximal (upstream and downstream)
#

# variant_anno_refseq_proximal_table
#
# Make table of proximal hits for RefSeq annotations.
rule variant_anno_refseq_proximal_table:
    input:
        bed='temp/variant/caller/{sourcename}/{sample}/{filter}/all/anno/refseq-prox/refseq-prox-{direction}-{flank}_{vartype}_{svtype}.bed.gz',
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/refseq-prox/refseq-prox-{direction}-{flank}_{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
        direction='up|dn',
        vartype='sv|indel|snv|rgn|sub',
        svtype='ins|del|inv|dup|snv|rgn|sub'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0, na_values='.')

        # Filter non-intersected records
        df = df.loc[~ pd.isnull(df['name'])]

        # Calculate distance to RefSeq annotation
        if df.shape[0] > 0:
            df['DISTANCE'] = df.apply(lambda row: _variant_anno_refseq_prox_get_distance(row, wildcards.direction), axis=1)
            df = df.loc[df['DISTANCE'] > 0]
        else:
            df['DISTANCE'] = []

        # Subset table
        df = df.loc[: , ('ID', 'comName', 'name', 'strand', 'DISTANCE')]
        df.columns = ('ID', 'GENE', 'GENE_ID', 'STRAND', 'DISTANCE')

        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

# variant_anno_refseq_proximal_intersect
#
# Intersect variants with proximal RefSeq annotations.
rule variant_anno_refseq_proximal_intersect:
    input:
        sv_bed='temp/variant/caller/{sourcename}/{sample}/{filter}/all/bed4/{vartype}_{svtype}.bed.gz',
        rs_bed='data/anno/refseq/region/flank_{flank}/refseq_{direction}.bed.gz'
    output:
        bed=temp('temp/variant/caller/{sourcename}/{sample}/{filter}/all/anno/refseq-prox/refseq-prox-{direction}-{flank}_{vartype}_{svtype}.bed.gz')
    wildcard_constraints:
        direction='up|dn'
    shell:
        """{{\n"""
        """    {{\n"""
        """        head -n 1 <(zcat {input.sv_bed});\n"""
        """        head -n 1 <(zcat {input.rs_bed});\n"""
        """    }} | \n"""
        """    awk -vORS="\\t" {{print}} | \n"""
        """    sed -re 's/\\s+/\\t/g' | \n"""
        """    sed -re 's/\\s*$/\\n/';\n"""
        """    bedtools intersect -a {input.sv_bed} -b {input.rs_bed} -loj;\n"""
        """}} | gzip > {output.bed}"""
