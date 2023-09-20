"""
Import external BED with variant calls for samples.
"""

#############
### Rules ###
#############

# variant_caller_extern_filter_bed
#
# Apply region filter to variant BED.
rule variant_caller_extern_get_bed:
    input:
        bed=lambda wildcards: svpoplib.rules.sample_table_entry(
            wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='bed'
        )['DATA'],
        fa=lambda wildcards: svpoplib.rules.get_bed_fa_input(
            svpoplib.rules.sample_table_entry(
                wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='bed'
            ),
            wildcards,
            default=[]
        )
    output:
        bed=temp('temp/variant/caller/bed/{sourcename}/{sample}/all/all/bed/pre_filter/{vartype}_{svtype}.bed.gz'),
        fa=temp('temp/variant/caller/bed/{sourcename}/{sample}/all/all/bed/pre_filter/fa/{vartype}_{svtype}.fa.gz')
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|rgn|snv',
    run:

        # Get entry and FASTA file name
        sample_entry = svpoplib.rules.sample_table_entry(wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='bed')

        if sample_entry['TYPE'] != 'bed':
            raise RuntimeError('Cannot process non-bed type: ' + sample_entry['TYPE'])

        fa_file_name = svpoplib.rules.get_bed_fa_input(sample_entry, wildcards)

        # Get parameters
        dedup = svpoplib.util.as_bool(sample_entry['PARAMS'].get('dedup', False), True)

        # Read
        df = pd.read_csv(input.bed, sep='\t')

        # Check for required fields
        required_cols = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN']

        if wildcards.vartype == 'snv':
            required_cols += ['REF', 'ALT']

        missing_cols = [col for col in required_cols if col not in df.columns]

        if missing_cols:
            raise RuntimeError('Missing {} column(s) from input BED "{}": {}'.format(
                len(missing_cols), input.bed, ', '.join(missing_cols))
            )

        na_cols = [col for col in required_cols if np.any(pd.isnull(df[col]))]

        if na_cols:
            raise RuntimeError('Found {} required column(s) with NA values from input BED "{}": {}'.format(
                len(na_cols), input.bed, ', '.join(na_cols))
            )

        # Positive SVLEN (in case it slips by)
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