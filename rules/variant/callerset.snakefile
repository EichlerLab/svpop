"""
Merge variants from multiple callers for a single sample.
"""

import Bio.bgzf
import Bio.SeqIO
import pandas as pd

import svpoplib

global expand
global temp
global shell


#
# Rules
#

# Merge callerset variants from multiple sources for one sample.
rule variant_callerset_merge:
    input:
        bed=expand(
            'temp/variant/callerset/{{sourcename}}/{{sample}}/{{filter}}/all/bed/{{vartype}}_{{svtype}}/chrom_{chrom}.bed.gz',
            chrom=svpoplib.ref.get_df_fai(config['reference_fai']).index
        )
    output:
        bed='results/variant/callerset/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz'
    wildcard_constraints:
        vartype='sv|indel|snv',
        svtype='ins|del|inv|dup|snv',
        sample='[^/]+'
    run:

        df = pd.concat(
            [pd.read_csv(bed_file_name, sep='\t') for bed_file_name in input.bed],
            sort=False,
            axis=0
        ).sort_values(
            ['#CHROM', 'POS']
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# Merge one chromosome.
rule variant_callerset_merge_chrom:
    input:
        bed=lambda wildcards:
            svpoplib.callerset.get_caller_set_input(
                wildcards.sourcename,
                'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
                config,
                wildcards
            ),
        fa=lambda wildcards:
            svpoplib.callerset.get_caller_set_input(
                wildcards.sourcename,
                'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz',
                config,
                wildcards
            ) if svpoplib.callerset.is_read_seq(wildcards, config) else []
    output:
        bed=temp('temp/variant/callerset/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}/chrom_{chrom}.bed.gz')
    threads: 8
    run:

        # Get entry
        callerset_entry = svpoplib.callerset.get_config_entry(wildcards.sourcename, config)

        merge_strategy = svpoplib.sampleset.get_merge_strategy(callerset_entry, wildcards.vartype, wildcards.svtype, config)

        # Merge
        df = svpoplib.svmerge.merge_variants(
            input.bed,
            callerset_entry['name_list'],
            merge_strategy['strategy'],
            fa_list = input.fa if input.fa else None,
            subset_chrom=wildcards.chrom,
            threads=threads,
            ref_filename=config['reference']
        )

        # Rename MERGE for callerset
        rename_dict = {
            'MERGE_SRC': 'CALLERSET_SRC',
            'MERGE_SRC_ID': 'CALLERSET_SRC_ID',
            'MERGE_AC': 'CALLERSET_N',
            'MERGE_AF': 'CALLERSET_PROP',
            'MERGE_SAMPLES': 'CALLERSET_LIST',
            'MERGE_VARIANTS': 'CALLERSET_VARIANTS',
            'MERGE_RO': 'CALLERSET_RO',
            'MERGE_DIST': 'CALLERSET_DIST',
            'MERGE_SZRO': 'CALLERSET_SZRO',
            'MERGE_OFFSZ': 'CALLERSET_OFFSZ'
        }

        df.columns = [rename_dict.get(col, col) for col in df.columns]

        if 'DISC_CLASS' in df.columns:
            del(df['DISC_CLASS'])

        # Bylen to byref
        if df.shape[0] > 0:
            df['END'] = df.apply(lambda row: (row['POS'] + 1) if row['SVTYPE'] == 'INS' else row['END'], axis=1)

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# Merge annotations for a caller set.
rule variant_callerset_fa:
    input:
        bed='results/variant/callerset/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        tsv=lambda wildcards: svpoplib.callerset.get_caller_set_input(
            wildcards.sourcename,
            'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz',
            config,
            wildcards
        )
    output:
        fa='results/variant/callerset/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz',
        fai='results/variant/callerset/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz.fai'
    wildcard_constraints:
        filter='\\w+',
        annodir='[a-zA-Z0-9\\.\\-]+',
        sample='[a-zA-Z0-9\\.-]+',
        svtype='ins|del|inv|dup|sub|rgn|snv',
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0, usecols=('ID', 'CALLERSET_SRC_ID', 'CALLERSET_SRC'))

        # Get sources.
        callerset_input = svpoplib.callerset.get_caller_set_input(
            wildcards.sourcename,
            'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz',
            config,
            wildcards
        )

        # Get entry
        callerset_entry = svpoplib.callerset.get_config_entry(wildcards.sourcename, config)

        # Open FA and write
        with Bio.bgzf.BgzfWriter(output.fa, 'wt') as out_file:
            Bio.SeqIO.write(
                svpoplib.callerset.fa_iter(df, callerset_entry, callerset_input),
                out_file,
                'fasta'
            )

        shell("""samtools faidx {output.fa}""")


# Merge annotations for a caller set.
rule variant_callerset_anno_merge:
    input:
        bed='results/variant/callerset/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        anno=lambda wildcards: svpoplib.callerset.get_caller_set_input(
            wildcards.sourcename,
            'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/{annodir}/{annotype}_{vartype}_{svtype}.{ext}.gz',
            config,
            wildcards
        )
    output:
        anno='results/variant/callerset/{sourcename}/{sample}/{filter}/all/anno/{annodir}/{annotype}_{vartype}_{svtype}.{ext}.gz'
    wildcard_constraints:
        filter='\\w+',
        annodir='[a-zA-Z0-9\\.\\-]+',
        sample='[a-zA-Z0-9\\.-]+',
        svtype='ins|del|inv|dup|sub|rgn|snv',
        ext='tsv|bed'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Get sources.
        callerset_input = svpoplib.callerset.get_caller_set_input(
            wildcards.sourcename,
            'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/{annodir}/{annotype}_{vartype}_{svtype}.{ext}.gz',
            config,
            wildcards
        )

        # Get entry
        callerset_entry = svpoplib.callerset.get_config_entry(wildcards.sourcename, config)

        # Merge annotations
        df_merge = svpoplib.callerset.merge_annotations(
            df, callerset_input, callerset_entry
        )

        # Write
        df_merge.to_csv(output.anno, sep='\t', index=False, compression='gzip')
