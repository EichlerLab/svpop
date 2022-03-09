# Write callsets as a VCF

VARIANT_VCF_INFO_HEADER = {
    'ID': ('ID', '1', 'String', 'Variant ID'),
    'SVTYPE': ('SVTYPE', '1', 'String', 'Variant type'),
    'SVLEN': ('SVLEN', '.', 'String', 'Variant length'),
    'SEQ': ('SEQ', '.', 'String', 'Variant sequence')
}

VARIANT_VCF_ALT_HEADER = {
    'ins': ('INS', 'Sequence insertion'),
    'del': ('DEL', 'Sequence deletion'),
    'inv': ('INV', 'Inversion'),
    'dup': ('DUP', "Duplication")
}

VARIANT_VCF_SVTYPE_LIST = {
    'insdel': ('ins', 'del'),
    'insdelinv': ('ins', 'del', 'inv')
}


# vcf_write_vcf
#
# Make VCF headers.
rule vcf_write_vcf:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz',
        fa='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/fa/{vartype}_{svtype}.fa.gz',
        ref_tsv='data/ref/contig_info.tsv.gz'
    output:
        vcf='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/vcf/{vartype}_{svtype}_{altfmt}.vcf.gz'
    wildcard_constraints:
        altfmt='alt|sym|sym-noseq'
    run:

        # Check assembly name
        if wildcards.sample in {'#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'}:
            raise RuntimeError(f'Sample name conflicts with a VCF header column name: {wildcards.sample}')

        # Check alt format
        if wildcards.altfmt == 'alt':
            symbolic_alt = False
            alt_seq = True

        elif wildcards.altfmt == 'sym':
            symbolic_alt = True
            alt_seq = True

        elif wildcards.altfmt == 'sym-noseq':
            symbolic_alt = True
            alt_seq = True

        else:
            raise RuntimeError(f'Unknown alt format wildcard (altfmt): {wildcards.altfmt}')

        # Check for symbolic ALT compatibility
        if symbolic_alt:
            if wildcards.svtype not in {'ins', 'del', 'inv', 'insdel', 'insdelinv', 'dup'}:
                raise RuntimeError(f'Symbolic ALT output is not supported for input variant type: {wildcards.svtype}')
        else:
            if wildcards.svtype not in {'ins', 'del', 'insdel', 'snv'}:
                raise RuntimeError(f'ALT output (sequence in ALT column, not symbolic) is not supported for input variant type: {wildcards.svtype}')

        # Setup header lists
        info_header_id_list = list()

        # Read variants
        df = pd.read_csv(input.bed, sep='\t')

        # Read sequence from FASTA
        if alt_seq:

            # FASTA cannot be empty
            if os.stat(input.fa).st_size == 0:
                raise RuntimeError(f'Missing sequence data to construct VCF: Variant FASTA is empty: {input.fa}')

            # Get Sequence from FASTA and assign to SEQ column (not for SNVs)
            with gzip.open(input.fa, 'rt') as fa_in:
                df_seq_dict = {
                    record.name: str(record.seq) for record in SeqIO.parse(fa_in, 'fasta')
                }

            # Assign SEQ to records
            df.set_index('ID', inplace=True)

            df['SEQ'] = pd.Series(df_seq_dict)

            del(df_seq_dict)

            df.reset_index(inplace=True)

            # Check for missing sequences
            df_null_seq = df.loc[pd.isnull(df['SEQ'])]

            if df_null_seq.shape[0] > 0:
                id_list = ', '.join(df_null_seq.iloc[:3]['ID'])

                raise RuntimeError(
                    'Missing FASTA sequence for {} variants (vartype={}, svtype={}): {}{}'.format(
                        df_null_seq.shape[0],
                        wildcards.vartype,
                        wildcards.svtype,
                        ', '.join([str(val) for val in df_null_seq.iloc[:3]['ID']]),
                        '...' if df_null_seq.shape[0] > 3 else ''
                    )
                )

        # Add VARTYPE
        df['VARTYPE'] = wildcards.vartype.upper()

        # SVTYPE
        df['SVTYPE'] = df['SVTYPE'].apply(lambda val: val.upper())

        # Reformat fields for INFO
        df['SVLEN'] = df.apply(lambda row: np.abs(row['SVLEN']) * (-1 if row['SVTYPE'] == 'DEL' else 1), axis=1)

        # Add GT if missing
        if 'GT' not in df.columns:
            df['GT'] = './.'

        # INFO: Base
        df['INFO'] = df.apply(lambda row: 'ID={ID};SVTYPE={SVTYPE}'.format(**row), axis=1)

        info_header_id_list.append('ID')
        info_header_id_list.append('SVTYPE')

        # INFO: Add SV/INDEL annotations
        df['INFO'] = df.apply(lambda row: row['INFO'] + (';SVLEN={SVLEN}'.format(**row)) if row['SVTYPE'] != 'SNV' else '', axis=1)
        info_header_id_list.append('SVLEN')

        # REF
        df['REF'] = list(svpoplib.vcf.ref_base(df, config['reference']))

        # ALT
        if wildcards.svtype != 'snv':
            if symbolic_alt:
                df['ALT'] = df['SVTYPE'].apply(lambda val: f'<{val}>')

                if alt_seq:
                    df['INFO'] = df.apply(lambda row: row['INFO'] + ';SEQ={SEQ}'.format(**row), axis=1)
                    info_header_id_list.append('SEQ')

            else:

                df['REF'] = df.apply(lambda row:
                    ((row['REF'] + row['SEQ']) if row['POS'] > 0 else (row['SEQ'] + row['REF'])) if row['SVTYPE'] == 'DEL' else row['REF'],
                    axis=1
                )

                df['ALT'] = df.apply(lambda row:
                    ((row['REF'] + row['SEQ']) if row['POS'] > 0 else (row['SEQ'] + row['REF'])) if row['SVTYPE'] == 'INS' else row['REF'][0],
                    axis=1
                )

        else:
            # Fix position for SNVs (0-based BED to 1-based VCF)
            df['POS'] += 1

        # Remove SEQ
        if alt_seq:
            del df['SEQ']

        # No masked bases in REF/ALT
        df['REF'] = df['REF'].apply(lambda val: val.upper())
        df['ALT'] = df['ALT'].apply(lambda val: val.upper())

        # Save columns needed for VCF
        df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO', 'GT']].copy()

        # Sort
        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # INFO headers
        info_header_list = [VARIANT_VCF_INFO_HEADER[val] for val in info_header_id_list]

        # ALT headers
        if symbolic_alt:
            alt_header_list = [VARIANT_VCF_ALT_HEADER[svtype] for svtype in VARIANT_VCF_SVTYPE_LIST.get(wildcards.svtype, (wildcards.svtype, ))]
        else:
            alt_header_list = list()

        # QUAL, FILTER, FORMAT
        if 'QUAL' not in df.columns:
            df['QUAL'] = '.'

        if 'FILTER' in df.columns:
            raise RuntimeError('FILTER is defined in dataframe, but FILTER headers are not implemented')

        filter_header_list = list()
        df['FILTER'] = '.'

        df['FORMAT'] = 'GT'

        format_header_list = [
            ('GT', '1', 'String', 'Genotype')
        ]

        # VCF order
        df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'GT']]
        df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', wildcards.sample]

        # Write
        df_ref = pd.read_csv(input.ref_tsv, sep='\t')

        with Bio.bgzf.open(output.vcf, 'wt') as out_file:
            for line in svpoplib.vcf.header_list(
                df_ref,
                info_header_list,
                format_header_list,
                alt_header_list,
                filter_header_list,
                variant_source=f'SV-Pop {svpoplib.constants.VERSION}',
                ref_file_name=os.path.basename(config.get('reference'))
            ):
                out_file.write(line)

            out_file.write('\t'.join(df.columns))
            out_file.write('\n')

            for index, row in df.iterrows():
                out_file.write('\t'.join(row.astype(str)))
                out_file.write('\n')

        # Write tabix index if possible
        try:
            shell("""tabix {output.vcf} && touch -r {output.vcf} {output.vcf}.tbi""")
        except:
            pass
