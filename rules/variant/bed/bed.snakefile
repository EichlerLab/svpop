"""
Import external BED with variant calls for samples.
"""

import collections
import numpy as np
import os
import pandas as pd
import svpoplib

global SAMPLE_TABLE
global temp

from Bio import SeqIO
import Bio.bgzf

###################
### Definitions ###
###################

def _variant_caller_extern_get_bed(wildcards):

    in_file_pattern = svpoplib.rules.sample_table_entry(
        wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='bed', format_data=False
    )['DATA']

    in_file = in_file_pattern.format(**wildcards)

    if not os.path.isfile(in_file) and wildcards.vartype in {'sv', 'indel'} and wildcards.svtype in {'ins', 'del'}:
        in_file_svindel = svpoplib.util.format_cards(in_file_pattern, vartype='svindel').format(**wildcards)

        if os.path.isfile(in_file_svindel):
            in_file = in_file_svindel

    return in_file

def _variant_caller_extern_get_fa(wildcards):

    fa_input = svpoplib.rules.get_bed_fa_input(
        svpoplib.rules.sample_table_entry(
            wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='bed'
        ),
        wildcards,
        default=None
    )

    if fa_input is not None and not os.path.isfile(fa_input) and wildcards.vartype in {'sv', 'indel'} and wildcards.svtype in {'ins', 'del'}:
        fa_input_svindel = svpoplib.rules.get_bed_fa_input(
            svpoplib.rules.sample_table_entry(
                wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='bed',
            ),
            wildcards,
            default=None,
            vartype='svindel'
        )

        if fa_input_svindel is not None and os.path.isfile(fa_input_svindel):
            fa_input = fa_input_svindel

    if fa_input is None or not os.path.isfile(fa_input):
        return []

    return fa_input


#############
### Rules ###
#############

# variant_caller_extern_filter_bed
#
# Apply region filter to variant BED.
rule variant_caller_extern_get_bed:
    input:
        bed=_variant_caller_extern_get_bed,
        fa=_variant_caller_extern_get_fa
    output:
        bed=temp('temp/variant/caller/bed/{sourcename}/{sample}/all/all/bed/pre_filter/{vartype}_{svtype}.bed.gz'),
        fa=temp('temp/variant/caller/bed/{sourcename}/{sample}/all/all/bed/pre_filter/fa/{vartype}_{svtype}.fa.gz')
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|rgn|snv'
    run:

        # Get entry and FASTA file name
        sample_entry = svpoplib.rules.sample_table_entry(wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='bed')

        if sample_entry['TYPE'] != 'bed':
            raise RuntimeError('Cannot process non-bed type: ' + sample_entry['TYPE'])

        #fa_file_name = svpoplib.rules.get_bed_fa_input(sample_entry, wildcards)
        fa_file_name = input.fa

        if not fa_file_name:
            fa_file_name = None

        # Get parameters
        dedup = svpoplib.util.as_bool(sample_entry['PARAMS'].get('dedup', False), True)

        # Read
        df = pd.read_csv(input.bed, sep='\t')

        # Check for required fields
        required_cols = ['#CHROM', 'POS', 'END', 'SVTYPE']

        if wildcards.vartype == 'snv':
            required_cols += ['REF', 'ALT']
        else:
            required_cols += ['SVLEN']

        missing_cols = [col for col in required_cols if col not in df.columns]

        if missing_cols:
            raise RuntimeError('Missing {} column(s) from input BED "{}": {}'.format(
                len(missing_cols), input.bed, ', '.join(missing_cols))
            )

        # Set missing IDs
        if 'ID' not in df.columns or np.all(pd.isnull(df['ID'])):
            df['ID'] = svpoplib.variant.get_variant_id(df)

        na_cols = [col for col in required_cols + ['ID'] if np.any(pd.isnull(df[col]))]

        if na_cols:
            raise RuntimeError('Found {} required column(s) with NA values from input BED "{}": {}'.format(
                len(na_cols), input.bed, ', '.join(na_cols))
            )

        # Positive SVLEN (in case it slips by)
        if 'SVLEN' in df.columns:
            df['SVLEN'] = np.abs(df['SVLEN'])

        # Filter by type
        df = df.loc[df['SVTYPE'] == wildcards.svtype.upper()]

        if wildcards.vartype == 'sv' and wildcards.svtype in {'ins', 'del'}:
            df = df.loc[df['SVLEN'] >= 50]

        if wildcards.vartype == 'indel' and wildcards.svtype in {'ins', 'del'}:
            df = df.loc[df['SVLEN'] < 50]

        # Check for duplicate IDs
        if not dedup:
            dup_id = [col for col, count in collections.Counter(df['ID']).items() if count > 1]

            if dup_id:
                dup_id_list = ', '.join(sorted(dup_id)[:3])

                raise RuntimeError('Found {} duplicate ID values from input BED "{}": {}{}'.format(
                    len(dup_id), input.bed, ', '.join(sorted(dup_id)[:3]), '...' if len(dup_id) > 3 else '')
                )
        else:
            df.drop_duplicates('ID', inplace=True)

        # Write sequences
        if 'fa_incomplete' in sample_entry['PARAMS']:
            fa_incomplete = svpoplib.util.as_bool(sample_entry['PARAMS'].get('fa_incomplete'), none_val=True)
        else:
            fa_incomplete = False

        if fa_file_name is not None and 'SEQ' in df.columns:
            raise RuntimeError('Found variant sequences in FASTA "{}" and the SEQ column in the BED file (must select sequences from one source)')

        if fa_file_name is not None and os.stat(fa_file_name).st_size == 0:
            fa_file_name = None  # Emtpy FASTA is equivalent to missing

        if fa_file_name is not None:
            with Bio.bgzf.BgzfWriter(output.fa, 'wt') as out_file:
                SeqIO.write(
                    svpoplib.seq.fa_to_record_iter(
                        fa_file_name,
                        set(df['ID']),
                        require_all=not fa_incomplete
                    ),
                    out_file,
                    'fasta'
                )

        elif 'SEQ' in df.columns:
            with Bio.bgzf.BgzfWriter(output.fa, 'wt') as out_file:
                SeqIO.write(
                    svpoplib.seq.bed_to_seqrecord_iter(df),
                    out_file,
                    'fasta'
                )

            del(df['SEQ'])

        else:

            # Empty FASTA
            with open(output.fa, 'w'):
                pass

        # Write variant BED
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')