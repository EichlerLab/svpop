"""
Read from PAV BED files per haplotype.
"""

def _variant_pavbedhap_get_var_bed(wildcards):

    # Check sample name
    match = re.search(r'^(.*)-h(\d+)$', wildcards.sample)

    if not match:
        raise RuntimeError(f'Variant callset type "pavbedhap" requires sample names to be appended with the haplotype (e.g. "SAMPLE" haplotype 1 is "SAMPLE-h1"): Received sample: {wildcards.sample}')

    sample = match[1]
    hap = f'h{match[2]}'

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

    return os.path.join(
        sample_entry['DATA'],
        'results',
        sample_fmt,
        'bed', 'pre_merge',
        hap,
        '{vartype}_{svtype}.bed.gz'.format(vartype=vartype, svtype=wildcards.svtype)
    )

# variant_pavbedhap_bed
#
# Get variants from PAV BED files.
rule variant_pavbedhap_bed:
    input:
        bed=_variant_pavbedhap_get_var_bed
    output:
        bed=temp('temp/variant/caller/pavbedhap/{sourcename}/{sample}/all/all/bed/pre_filter/{vartype}_{svtype}.bed.gz'),
        fa=temp('temp/variant/caller/pavbedhap/{sourcename}/{sample}/all/all/bed/pre_filter/fa/{vartype}_{svtype}.fa.gz')
    run:

        # Get sample entry parameters
        sample_entry = svpoplib.rules.sample_table_entry(
            wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='pavbedhap'
        )

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
        df = pd.read_csv(input.bed, sep='\t')

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
            # Write sequence FASTA
            with Bio.bgzf.BgzfWriter(output.fa) as fa_file:
                Bio.SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(df), fa_file, 'fasta')

            del(df['SEQ'])

        else:
            # Empty file for SNVs
            with open(output.fa, 'wt') as out_file:
                pass

        # Write variant table
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')
