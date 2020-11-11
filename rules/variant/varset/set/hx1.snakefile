"""
HX1 PacBio SVs

Shi, L., Guo, Y., Dong, C., Huddleston, J., Yang, H., Han, X., … Wang, K. (2016). Long-read sequencing and de novo assembly of a Chinese genome. Nature Communications, 7(May), 12065. https://doi.org/10.1038/ncomms12065
"""

# varset_set_hx1_bed
#
# Make BED files
rule varset_set_hx1_bed:
    input:
        bed='temp/variant/varset/hx1/variants.bed.gz'
    output:
        bed_ins='temp/variant/varset/hx1/bed/all/all/all/byref/sv_ins.bed',
        bed_del='temp/variant/varset/hx1/bed/all/all/all/byref/sv_del.bed',
        bed_inv='temp/variant/varset/hx1/bed/all/all/all/byref/sv_inv.bed'
    run:

        # Read
        df = pd.read_csv(
            input.bed, sep='\t', header=None,
            names=('#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'SEQ', 'REPEAT', 'TSD', 'CONTIG', 'CONTIG_START', 'CONTIG_END', 'UK_1')
        )

        # Byref, fix SVTYPE and ID
        df['END'] = df.apply(lambda row: row['POS'] + (1 if row['SVTYPE'] == 'INS' else row['SVLEN']), axis=1)
        df['SVTYPE'] = df['SVTYPE'].apply(lambda val: val[:3].upper())
        df.insert(3, 'ID', analib.variant.get_variant_id(df))

        # Write
        df.loc[df['SVTYPE'] == 'INS'].to_csv(output.bed_ins, sep='\t', index=False)
        df.loc[df['SVTYPE'] == 'DEL'].to_csv(output.bed_del, sep='\t', index=False)
        df.loc[df['SVTYPE'] == 'INV'].to_csv(output.bed_inv, sep='\t', index=False)


# varset_set_hx1_dl
#
# Get variant BED file for SV calls from the PacBio assembly.
rule varset_set_hx1_dl:
    output:
        bed=temp('temp/variant/varset/hx1/variants.bed.gz')
    params:
        url=config['varset']['set']['hx1']['sv_bed']
    shell:
        """wget {params.url} -O {output.bed}; """
