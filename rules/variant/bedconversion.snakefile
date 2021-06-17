"""
BED file conversion rules.
"""

#
# BED4
#

# variant_bedconversion_bed4
#
# Make bed4 format (position and ID).
rule variant_bedconversion_bed4:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz'
    output:
        bed=temp('temp/variant/caller/{sourcename}/{sample}/{filter}/all/bed4/{vartype}_{svtype}.bed.gz')
    run:
        pd.read_csv(input.bed, sep='\t').loc[: , ['#CHROM', 'POS', 'END', 'ID']].to_csv(output.bed, sep='\t', index=False, compression='gzip')


#
# By-length conversion
#

# variant_bedconversion_bylen
#
# Adds SVLEN to POS for insertions to create a pseudo-region for intersections. Used for
# bedtools intesects when comparing callsets, but SV-Pop doesn't need it for its own intersects.
rule variant_bedconversion_bylen:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz'
    output:
        bed=temp('temp/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/bylen/{vartype}_{svtype}.bed.gz')
    run:

        if wildcards.svtype == 'ins':
            # Read
            df = pd.read_csv(input.bed, sep='\t', header=0)

            df['END'] = df['POS'] + df['SVLEN']
            df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

        else:
            shell(
                """ln -sf ../byref/{wildcards.vartype}_{wildcards.svtype}.bed {output.bed}; """
            )

# variant_bedconversion_bylen
#
# Adds SVLEN to POS for insertions to create a pseudo-region for intersections. Used for
# bedtools intesects when comparing callsets, but SV-Pop doesn't need it for its own intersects.
rule variant_bedconversion_uncompress:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz'
    output:
        bed=temp('temp/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed')
    run:

        if wildcards.svtype == 'ins':
            # Read
            df = pd.read_csv(input.bed, sep='\t', header=0)

            df['END'] = df['POS'] + df['SVLEN']
            df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

        else:
            shell(
                """ln -sf ../byref/{wildcards.vartype}_{wildcards.svtype}.bed {output.bed}; """
            )

#
# BED4 and by-length
#

# variant_bedconversion_bed4_bylen
#
# Variant byref to bylen.
rule variant_bedconversion_bed4_bylen:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz'
    output:
        bed=temp('temp/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed4/bylen/{vartype}_{svtype}.bed.gz')
    run:

        df = pd.read_csv(input.bed, sep='\t', header=0)

        if wildcards.svtype == 'ins':
            df['END'] = df['POS'] + df['SVLEN']

        df = df[['#CHROM', 'POS', 'END', 'ID']]

        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')
