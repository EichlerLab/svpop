"""
HGSVC variants.

Chaisson, M. J. P., Sanders, A. D., Zhao, X., Malhotra, A., Porubsky, D., Rausch, T., … Lee, C. (2017). Multi-platform discovery of haplotype-resolved structural variation in human genomes. BioRxiv. https://doi.org/10.1101/193144
"""



# varset_set_hgsvc1_bed_all_samples
#
# Merge samples using default merge parameters.
rule varset_set_hgsvc1_bed_all_samples:
    input:
        bed_HG00514='temp/variant/varset/hgsvc1/bed/HG00514/all/all/byref/sv_{svtype}.bed',
        bed_HG00733='temp/variant/varset/hgsvc1/bed/HG00733/all/all/byref/sv_{svtype}.bed',
        bed_NA19420='temp/variant/varset/hgsvc1/bed/NA19240/all/all/byref/sv_{svtype}.bed'
    output:
        bed=temp('temp/variant/varset/hgsvc1/bed/all/all/all/byref/sv_{svtype}.bed')
    wildcard_constraints:
        svtype='ins|del|inv'
    run:

        # Get configured merge definition
        config_def = 'nr:szro=50:offset=200'

        # Merge
        df = analib.svmerge.merge_variants(
            bed_list=[input.bed_HG00514, input.bed_HG00733, input.bed_NA19420],
            sample_names=['HG00514', 'HG00733', 'NA19240'],
            strategy=config_def,
            threads=6
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False)


# varset_set_hgsvcpilot_unified_bed
#
# Unified SV BED.
rule varset_set_hgsvc1_bed:
    input:
        bed=config['varset']['set']['hgsvc1']['bed']
    output:
        bed=temp('temp/variant/varset/hgsvc1/bed/{sample}/all/all/byref/sv_{svtype}.bed')
    wildcard_constraints:
        sample='HG00514|HG00733|NA19240',
        svtype='ins|del|inv'
    shell:
         """zcat {input.bed} > {output.bed}"""

