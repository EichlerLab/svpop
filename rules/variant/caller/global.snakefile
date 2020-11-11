"""
Global rules for callers.
"""

###################
### Definitions ###
###################


#############
### Rules ###
#############

#
# BED4
#

# variant_caller_global_bed4
#
# Make bed4 format (position and ID).
rule variant_caller_global_bed4:
    input:
        bed='results/variant/caller/{caller}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed'
    output:
        bed='results/variant/caller/{caller}/bed4/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed'
    run:
        pd.read_csv(input.bed, sep='\t').loc[: , ['#CHROM', 'POS', 'END', 'ID']].to_csv(output.bed, sep='\t', index=False)


#
# Group variant calls and annotations - sv_all (ins/del/inv)
#

# variant_caller_global_concat_anno_sv_all
#
# Concatenate annotations.
rule variant_caller_global_concat_anno_sv_all:
    input:
        tab_ins='results/variant/{sourcetype}/{sourcename}/anno/{sample}/{varset}/{filter}/{annodir}/{annotype}_sv_ins.{ext}',
        tab_del='results/variant/{sourcetype}/{sourcename}/anno/{sample}/{varset}/{filter}/{annodir}/{annotype}_sv_del.{ext}',
        tab_inv='results/variant/{sourcetype}/{sourcename}/anno/{sample}/{varset}/{filter}/{annodir}/{annotype}_sv_inv.{ext}'
    output:
        tab_all='results/variant/{sourcetype}/{sourcename}/anno/{sample}/{varset}/{filter}/{annodir}/{annotype}_sv_all.{ext,tab|bed}'
    run:

        # Concat
        df = analib.pd.concat_frames(
            [pd.read_csv(file_name, sep='\t', header=0) for file_name in input]
        )

        # Sort
        sort_cols = [col for col in df.columns if col in {'#CHROM', 'POS', 'END', 'ID'}]

        if len(sort_cols) > 0:
            df = df.sort_values(sort_cols)

        # Write
        df.to_csv(output.tab_all, sep='\t', index=False)

# variant_caller_global_concat_sv_all
#
# Concatenate all SV types (ins, del, and inv).
rule variant_caller_global_concat_sv_all:
    input:
        sv_ins='results/variant/{sourcetype}/{sourcename}/bed/{sample}/all/{filter}/byref/sv_ins.bed',
        sv_del='results/variant/{sourcetype}/{sourcename}/bed/{sample}/all/{filter}/byref/sv_del.bed',
        sv_inv='results/variant/{sourcetype}/{sourcename}/bed/{sample}/all/{filter}/byref/sv_inv.bed'
    output:
        all_bed='results/variant/{sourcetype}/{sourcename}/bed/{sample}/all/{filter}/byref/sv_all.bed'
    run:

        analib.pd.concat_frames(
            [pd.read_csv(file_name, header=0, sep='\t') for file_name in input]
        ).sort_values(
            ['#CHROM', 'POS', 'END', 'ID']
        ).to_csv(
            output.all_bed, sep='\t', index=False
        )


#
# Group variant calls and annotations - insdel
#

# variant_caller_global_concat_anno_insdel
#
# Concatenate annotations.
rule variant_caller_global_concat_anno_insdel:
    input:
        tab_ins='results/variant/{sourcetype}/{sourcename}/anno/{sample}/{varset}/{filter}/{annodir}/{annotype}_{vartype}_ins.{ext}',
        tab_del='results/variant/{sourcetype}/{sourcename}/anno/{sample}/{varset}/{filter}/{annodir}/{annotype}_{vartype}_del.{ext}'
    output:
        tab_all='results/variant/{sourcetype}/{sourcename}/anno/{sample}/{varset}/{filter}/{annodir}/{annotype}_{vartype}_insdel.{ext,tab|bed}'
    run:

        # Concat
        df = analib.pd.concat_frames(
            [pd.read_csv(file_name, sep='\t', header=0) for file_name in input]
        )

        # Sort
        sort_cols = [col for col in df.columns if col in {'#CHROM', 'POS', 'END', 'ID'}]

        if len(sort_cols) > 0:
            df = df.sort_values(sort_cols)

        # Write
        df.to_csv(output.tab_all, sep='\t', index=False)

# variant_caller_global_concat_insdel
#
# Concatenate all SV types (ins, del, and inv).
rule variant_caller_global_concat_insdel:
    input:
        sv_ins='results/variant/{sourcetype}/{sourcename}/bed/{sample}/all/{filter}/byref/{vartype}_ins.bed',
        sv_del='results/variant/{sourcetype}/{sourcename}/bed/{sample}/all/{filter}/byref/{vartype}_del.bed'
    output:
        all_bed='results/variant/{sourcetype}/{sourcename}/bed/{sample}/all/{filter}/byref/{vartype}_insdel.bed'
    run:

        analib.pd.concat_frames(
            [pd.read_csv(file_name, sep='\t', header=0) for file_name in input]
        ).sort_values(
            ['#CHROM', 'POS', 'END', 'ID']
        ).to_csv(
            output.all_bed, sep='\t', index=False
        )

# variant_caller_global_concat_insdel_fa
#
# Concatenate insertion and deletion records into one FASTA.
rule variant_caller_global_concat_insdel_fa:
    input:
        fa_ins='results/variant/{sourcetype}/{sourcename}/fasta/{sample}/{varset}/{filter}/{vartype}_ins.fa.gz',
        fa_del='results/variant/{sourcetype}/{sourcename}/fasta/{sample}/{varset}/{filter}/{vartype}_del.fa.gz'
    output:
        fa='results/variant/{sourcetype}/{sourcename}/fasta/{sample}/{varset}/{filter}/{vartype}_insdel.fa.gz'
    run:

        def get_fa_iter():
            for in_file_name in (input.fa_ins, input.fa_del):
                with gzip.open(in_file_name, 'rt') as in_file:
                    for record in SeqIO.parse(in_file, 'fasta'):
                        yield record

        with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
            SeqIO.write(get_fa_iter(), out_file, 'fasta')

# variant_caller_global_concat_all_fa
#
# Concatenate insertion, deletion, and inversion records into one FASTA.
rule variant_caller_global_concat_all_fa:
    input:
        fa_ins='results/variant/{sourcetype}/{sourcename}/fasta/{sample}/{varset}/{filter}/{vartype}_ins.fa.gz',
        fa_del='results/variant/{sourcetype}/{sourcename}/fasta/{sample}/{varset}/{filter}/{vartype}_del.fa.gz',
        fa_inv='results/variant/{sourcetype}/{sourcename}/fasta/{sample}/{varset}/{filter}/{vartype}_inv.fa.gz'
    output:
        fa='results/variant/{sourcetype}/{sourcename}/fasta/{sample}/{varset}/{filter}/{vartype}_all.fa.gz'
    run:

        def get_fa_iter():
            for in_file_name in (input.fa_ins, input.fa_del, input.fa_inv):
                with gzip.open(in_file_name, 'rt') as in_file:
                    for record in SeqIO.parse(in_file, 'fasta'):
                        yield record

        with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
            SeqIO.write(get_fa_iter(), out_file, 'fasta')
