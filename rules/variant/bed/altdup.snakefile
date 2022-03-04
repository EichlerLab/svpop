# Make dup calls based on re-mapping INS sequences back to the reference. Useful for identifying the likely
# templated location and for matching INS calls with DUP calls.

def variant_bed_altdup_input_filename(wildcards, type):

    # Get BED (altmap input) or FA (original FASTA file)
    TYPE_PATTERN = {
        'altmap': 'results/variant/{sourcetype}/{sourcename}/{{sample}}/{{filter}}/all/anno/altmap/altmap-{aligner}_sv_ins.bed.gz',
        'bed': 'results/variant/{sourcetype}/{sourcename}/{{sample}}/{{filter}}/all/bed/sv_ins.bed.gz',
        'fa': 'results/variant/{sourcetype}/{sourcename}/{{sample}}/{{filter}}/all/bed/fa/sv_ins.fa.gz'
    }

    # Check type
    if type not in TYPE_PATTERN.keys():
        raise RuntimeError('Error getting input file for altdup input: Requested type "{type}" is unknown, should be "bed" or "fa"')

    # Get sample entry
    sample_entry = svpoplib.rules.sample_table_entry(
        wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, type='altdup'
    )

    # Get DATA tokens (sourcetype, sourcename, alt-map aligner)
    if (pd.isnull(sample_entry['DATA'])):
        raise RuntimeError('Alt-dup sample entry {} is missing the "DATA" field: Expected "sourcetype,sourcename,aligner"'.format(sample_entry['NAME']))

    data_tok = [val.strip() for val in sample_entry['DATA'].strip().split(',')]

    if len(data_tok) != 3 or np.any([val == '' for val in data_tok]):
        raise RuntimeError('Alt-dup sample entry {} has an unrecognized "DATA" field: Expected "sourcetype,sourcename,aligner": Found "{}"'.format(sample_entry['NAME'], sample_entry['DATA']))

    data_tok_dict = {key: val for key, val in zip(('sourcetype', 'sourcename', 'aligner'), data_tok)}

    # Return file
    return TYPE_PATTERN[type].format(**data_tok_dict)


# variant_bed_altdup_bed
#
# Get ALT-DUP BED file from INS.
rule variant_bed_altdup_bed:
    input:
        altmap=lambda wildcards: variant_bed_altdup_input_filename(wildcards, 'altmap'),
        bed=lambda wildcards: variant_bed_altdup_input_filename(wildcards, 'bed'),
        fa=lambda wildcards: variant_bed_altdup_input_filename(wildcards, 'altmap')
    output:
        bed=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/sv_dup.bed.gz'),
        fa=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/fa/sv_dup.fa.gz')
    run:

        # Read altmap and transform
        df = pd.read_csv(input.altmap, sep='\t')

        df['SVTYPE'] = 'DUP'
        df['SVLEN'] = df['MAP_LEN']

        del(df['MAP_LEN'])

        df = df.loc[df['SVLEN'] >= 50].copy()

        df = svpoplib.variant.order_variant_columns(df)
        df.sort_values(['#CHROM', 'POS', 'END', 'SVLEN'], inplace=True)

        df.index = list(df['ID'])  # Do not use set_index, it will change the index when ID is reassigned
        df.index.name = 'INDEX'

        df['ID'] = svpoplib.variant.get_variant_id(df)

        # Add BED columns
        df_bed = pd.read_csv(input.bed, sep='\t')
        df_bed.set_index('ID', inplace=True, drop=False)

        df_bed.columns = [f'INS_{col}' for col in df_bed.columns]

        df_bed = df_bed.reindex(df.index)

        df = pd.concat([df, df_bed], axis=1)

        del(df_bed)

        # Write BED
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

        # Write FASTA
        df = df[['#CHROM', 'POS', 'END', 'ID']].copy()

        with pysam.FastaFile(config['reference']) as in_file:
            df['SEQ'] = df.apply(lambda row:
                in_file.fetch(row['#CHROM'], row['POS'], row['END']),
                axis=1
            )

        with Bio.bgzf.BgzfWriter(output.fa) as out_file:
            SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(df), out_file, 'fasta')
