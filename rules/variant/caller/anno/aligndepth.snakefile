"""
Get alignment depth for each variant
"""

# variant_caller_anno_depth_merge
#
# Merge depths for all chromosomes.
rule variant_caller_anno_depth_merge:
    input:
        tab=expand(
            'temp/variant/caller/{{sourcename}}/anno/{{sample}}/all/{{filter}}/aligndepth/aligndepth-{{alignset}}-{{alignsample}}_{{vartype}}_{{svtype}}/{chrom}.tab',
            chrom=analib.ref.get_df_fai(config['reference_fai']).index
        )
    output:
        tab='results/variant/caller/{sourcename}/anno/{sample}/all/{filter}/aligndepth/aligndepth-{alignset}-{alignsample}_{vartype}_{svtype}.tab'
    wildcard_constraints:
        vartype='sv|indel|snv',
        svtype='ins|del|inv|dup|snv|rgn|sub',
    run:

        df = pd.concat(
            [pd.read_csv(bed_file_name, sep='\t') for bed_file_name in input.tab],
            axis=0,
            sort=False
        ).to_csv(
            output.tab, sep='\t', index=False
        )


# variant_caller_anno_depth_chrom
#
# Get alignment depth per chromosome.
rule variant_caller_anno_depth_chrom:
    input:
        bed='results/variant/caller/{sourcename}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed'
    output:
        tab=temp(
            'temp/variant/caller/{sourcename}/anno/{sample}/all/{filter}/aligndepth/aligndepth-{alignset}-{alignsample}_{vartype}_{svtype}/{chrom}.tab'
        )
    wildcard_constraints:
        vartype='sv|indel|snv',
        svtype='ins|del|inv|dup|snv|rgn|sub',
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', low_memory=False)

        # Subset to chromosome
        df = df.loc[df['#CHROM'] == wildcards.chrom]

        # Get depth
        df_depth = analib.anno.align.get_depth(
            df.loc[df['#CHROM'] == wildcards.chrom],
            wildcards.alignset,
            wildcards.alignsample,
            config
        )

        # Write
        df_depth.to_csv(output.tab, sep='\t', index=False)
