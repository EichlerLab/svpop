"""
Read from PAV BED files per haplotype.
"""

import os
import pandas as pd

import Bio.bgzf

import svpoplib

global SAMPLE_TABLE
global temp

def _variant_pavbedhap_get_version(sample_entry):

    if not pd.isnull(sample_entry['VERSION']):
        try:
            if svpoplib.util.cmp_ver(sample_entry['VERSION'], '3') >= 0:
                return 3
            else:
                return 1
        except ValueError as e:
            raise RuntimeError(f'Error comparing PAV version in sample table: {e}')

    # Assume latest version if not set
    return 3


def _variant_pavbedhap_get_var_bed(wildcards):

    # Check sample name
    sample_tok = wildcards.sample.rsplit('-', 1)

    if len(sample_tok) != 2:
        raise RuntimeError(f'Variant callset type "pavbedhap" requires sample names to be appended with the haplotype (e.g. sample "SAMPLE" with haplotype "h1" is "SAMPLE-h1"): Received sample with no dash separating the sample name and the haplotype name: "{wildcards.sample}"')

    sample = sample_tok[0]
    hap = sample_tok[1]

    # Set vartype and svtype
    if wildcards.svtype in {'ins', 'del'}:
        if wildcards.vartype not in {'sv', 'indel'}:
            raise RuntimeError(f'Unknown vartype "{wildcards.vartype}" for svtype "{wildcards.svtype}"')

        vartype = 'svindel'
    else:
        vartype = wildcards.vartype

    # Get entry from sample table
    sample_entry = svpoplib.rules.sample_table_entry(
        wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='pavbedhap'
    )

    # Set path
    if 'sample_pattern' in sample_entry['PARAMS']:
        sample_pattern = sample_entry['PARAMS']['sample_pattern']

        if '{sample}' not in sample_pattern:
            raise RuntimeError(f'Sample entry has a "sample_pattern" without a "{{sample}}" wildcard: {sample_pattern}')

        sample_fmt = sample_pattern.format(sample=sample)
    else:
        sample_fmt = sample

    # Get path elements
    pav_version = _variant_pavbedhap_get_version(sample_entry)

    if pav_version >= 3:
        return os.path.join(
            sample_entry['DATA'],
            'results',
            sample_fmt,
            'bed_hap', 'pass',
            hap,
            '{vartype}_{svtype}.bed.gz'.format(vartype=vartype, svtype=wildcards.svtype)
        )

    else:
        return os.path.join(
            sample_entry['DATA'],
            'results',
            sample_fmt,
            'bed', 'pre_merge',
            hap,
            '{vartype}_{svtype}.bed.gz'.format(vartype=vartype, svtype=wildcards.svtype)
        )


def _variant_pavbedhap_get_var_fa(wildcards):

    # Get entry from sample table
    sample_entry = svpoplib.rules.sample_table_entry(
        wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='pavbedhap'
    )

    pav_version = _variant_pavbedhap_get_version(sample_entry)

    if pav_version < 3 or wildcards.vartype == 'snv':
        return []

    # Check sample name
    sample_tok = wildcards.sample.rsplit('-', 1)

    if len(sample_tok) != 2:
        raise RuntimeError(f'Variant callset type "pavbedhap" requires sample names to be appended with the haplotype (e.g. sample "SAMPLE" with haplotype "h1" is "SAMPLE-h1"): Received sample with no dash separating the sample name and the haplotype name: "{wildcards.sample}"')

    sample = sample_tok[0]
    hap = sample_tok[1]

    # Set vartype and svtype
    if wildcards.svtype in {'ins', 'del'}:
        if wildcards.vartype not in {'sv', 'indel'}:
            raise RuntimeError(f'Unknown vartype "{wildcards.vartype}" for svtype "{wildcards.svtype}"')

        vartype = 'svindel'
    else:
        vartype = wildcards.vartype

    # Set path
    if 'sample_pattern' in sample_entry['PARAMS']:
        sample_pattern = sample_entry['PARAMS']['sample_pattern']

        if '{sample}' not in sample_pattern:
            raise RuntimeError(f'Sample entry has a "sample_pattern" without a "{{sample}}" wildcard: {sample_pattern}')

        sample_fmt = sample_pattern.format(sample=sample)
    else:
        sample_fmt = sample

    # Get path elements
    return os.path.join(
        sample_entry['DATA'],
        'results',
        sample_fmt,
        'bed_hap', 'pass',
        hap, 'fa',
        '{vartype}_{svtype}.fa.gz'.format(vartype=vartype, svtype=wildcards.svtype)
    )


# variant_pavbedhap_bed
#
# Get variants from PAV BED files.
rule variant_pavbedhap_bed:
    input:
        bed=_variant_pavbedhap_get_var_bed,
        fa=_variant_pavbedhap_get_var_fa
    output:
        bed=temp('temp/variant/caller/pavbedhap/{sourcename}/{sample}/all/all/bed/pre_filter/{vartype}_{svtype}.bed.gz'),
        fa=temp('temp/variant/caller/pavbedhap/{sourcename}/{sample}/all/all/bed/pre_filter/fa/{vartype}_{svtype}.fa.gz')
    run:

        # Get sample entry parameters
        sample_entry = svpoplib.rules.sample_table_entry(
            wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='pavbedhap'
        )

        pav_version = _variant_pavbedhap_get_version(sample_entry)

        # Get min/max SVLEN
        min_svlen = sample_entry['PARAMS'].get('min_svlen', None)
        max_svlen = sample_entry['PARAMS'].get('max_svlen', None)

        if min_svlen is not None:
            try:
                min_svlen = int(min_svlen)
            except ValueError:
                raise RuntimeError(f'Minimum SVLEN parameter ("min_svlen") is not a number: {min_svlen}')

        if max_svlen is not None:
            try:
                max_svlen = int(max_svlen)
            except ValueError:
                raise RuntimeError(f'Maximum SVLEN parameter ("max_svlen") is not a number: {max_svlen}')

        # Read
        df = pd.read_csv(input.bed, sep='\t', dtype={'#CHROM': str}, low_memory=False)

        # Input as svindel (separate)
        if wildcards.vartype == 'sv':
            df = df.loc[df['SVLEN'] >= 50]

        elif wildcards.vartype == 'indel':
            df = df.loc[df['SVLEN'] < 50]

        # Filter by size
        if wildcards.vartype in ('sv', 'indel'):
            if min_svlen is not None:
                df = df.loc[df['SVLEN'] >= min_svlen]

            if max_svlen is not None:
                df = df.loc[df['SVLEN'] <= min_svlen]

        # Write FASTA
        if wildcards.vartype != 'snv':

            if pav_version >= 3:
                # Subset FASTA file from PAV
                id_set = set(df['ID'])

                with Bio.bgzf.BgzfWriter(output.fa) as fa_file:
                    Bio.SeqIO.write(svpoplib.seq.fa_to_record_iter(input.fa, id_set), fa_file, 'fasta')

            else:
                # Pull from SEQ column in variant table
                with Bio.bgzf.BgzfWriter(output.fa) as fa_file:
                    Bio.SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(df), fa_file, 'fasta')

                del(df['SEQ'])

        else:
            # Empty file for SNVs
            with open(output.fa, 'wt') as out_file:
                pass

        # Write variant table
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')
