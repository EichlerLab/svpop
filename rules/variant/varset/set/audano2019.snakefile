"""
Variants for Audano 2019.

Audano, P.A., Sulovari, A., Graves-Lindsay, T.A., Cantsilieris, S., Sorensen, M., Welch, A.E., Dougherty, M.L., Nelson, B.J., Shah, A., Dutcher, S.K., et al. (2019). Characterizing the Major Structural Variant Alleles of the Human Genome. Cell 176, 1–13.
"""

# varset_set_audano2018_bed_sv
#
# Get SV BED
rule varset_set_audano2019_bed_sv:
    input:
        bed=os.path.join(SVPOP_DIR, config['varset']['set']['audano2019']['sv'])
    output:
        bed_ins=temp('temp/variant/varset/audano2019/bed/{sample}/all/all/byref/sv_ins.bed'),
        bed_del=temp('temp/variant/varset/audano2019/bed/{sample}/all/all/byref/sv_del.bed'),
        bed_inv=temp('temp/variant/varset/audano2019/bed/{sample}/all/all/byref/sv_inv.bed')
    run:

        bed_cols = [
            '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'REPEAT_TYPE',
            'MERGE_SOURCE', 'MERGE_AC', 'MERGE_AF', 'MERGE_SAMPLES', 'DISC_CLASS',
            'TANDEM', 'GT_PROP_CALL', 'GT_AN', 'GT_AF', 'VALIDATION'
        ]

        # Read
        df = pd.read_csv(input.bed, sep='\t', usecols=bed_cols)

        df = df.loc[:, bed_cols]

        df = df.fillna('.')

        # Subset
        if wildcards.sample != 'all':
            df = df.loc[
                df['MERGE_SAMPLES'].apply(
                    lambda sample_list: wildcards.sample in set(sample_list.split(','))
                )
            ]

            if df.shape[0] == 0:
                raise RuntimeError('No record for sample: ' + sample)

        # Set ID
        df['ID'] = df['ID'].apply(lambda val: val.split('_', 1)[1])

        # Write
        df.loc[df['SVTYPE'] == 'INS'].to_csv(output.bed_ins, sep='\t', index=False, float_format='%.4f')
        df.loc[df['SVTYPE'] == 'DEL'].to_csv(output.bed_del, sep='\t', index=False, float_format='%.4f')
        df.loc[df['SVTYPE'] == 'INV'].to_csv(output.bed_inv, sep='\t', index=False, float_format='%.4f')


# varset_set_audano2019_bed_dup
#
# Get SV BED
rule varset_set_audano2019_bed_dup:
    input:
        bed=os.path.join(SVPOP_DIR, config['varset']['set']['audano2019']['dup'])
    output:
        bed='temp/variant/varset/audano2019/bed/all/all/all/byref/sv_dup.bed'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t')

        df['ID'] = df['ID'].apply(lambda val: val.split('_', 1)[1])

        df = analib.variant.order_variant_columns(df)

        df = df.fillna('.')

        # Write
        df.to_csv(output.bed, sep='\t', index=False, float_format='%.4f')

rule varset_set_audano2019_bed_indel:
    input:
        bed=os.path.join(SVPOP_DIR, config['varset']['set']['audano2019']['indel'])
    output:
        bed_ins=temp('temp/variant/varset/audano2019/bed/all/all/all/byref/indel_ins.bed'),
        bed_del=temp('temp/variant/varset/audano2019/bed/all/all/all/byref/indel_del.bed')
    run:

        bed_cols = [
            '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
            'MERGE_SRC', 'MERGE_SAMPLES', 'MERGE_AC', 'MERGE_AF',
            'DISC_CLASS'
        ]

        # Read
        df = pd.read_csv(input.bed, sep='\t', usecols=bed_cols)

        # Write
        df.loc[df['SVTYPE'] == 'INS'].to_csv(output.bed_ins, sep='\t', index=False, float_format='%.4f')
        df.loc[df['SVTYPE'] == 'DEL'].to_csv(output.bed_del, sep='\t', index=False, float_format='%.4f')