"""
Variants for Huddleston 2017 (GR). Retrieved from dbVar.

Huddleston, J., Chaisson, M. J. P., Steinberg, K. M., Warren, W., Hoekzema, K., Gordon, D., … Eichler, E. E. (2017). Discovery and genotyping of structural variation from long-read haploid genome sequence data. Genome Research, 27(5), 677–685. https://doi.org/10.1101/gr.214007.116
"""

# varset_set_huddleston2017_variants_to_bed
#
# Write BED files
rule varset_set_huddleston2017_variants_to_bed:
    input:
        bed=os.path.join(SVPOP_DIR, 'files/varset/huddleston2017/sv_all.bed.gz')
    output:
        bed_ins='temp/variant/varset/huddleston2017/bed/all/all/all/byref/sv_ins.bed',
        bed_del='temp/variant/varset/huddleston2017/bed/all/all/all/byref/sv_del.bed',
        bed_inv='temp/variant/varset/huddleston2017/bed/all/all/all/byref/sv_inv.bed'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Write
        df.loc[df['SVTYPE'] == 'INS'].to_csv(output.bed_ins, sep='\t', index=False)
        df.loc[df['SVTYPE'] == 'DEL'].to_csv(output.bed_del, sep='\t', index=False)
        df.loc[df['SVTYPE'] == 'INV'].to_csv(output.bed_inv, sep='\t', index=False)
