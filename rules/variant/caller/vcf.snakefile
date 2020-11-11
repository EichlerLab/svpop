"""
Make VCF of variants
"""

def _variant_vcf_ref_base(subdf, ref_file_name):
    """
    Set REF to base before indel.

    :subdf: A subset of a dataframe to correct.
    """

    subdf = subdf.copy()

    # Open and update records
    with pysam.FastaFile(ref_file_name) as ref_file:
        for index, row in subdf.iterrows():
            subdf.loc[index, 'REF'] = ref_file.fetch(row['CHROM'], row['POS'] - 1, row['POS']).upper()

    # Return
    return subdf


# variant_vcf_create_vcf
#
# Create VCF for sample.
#
# wildcard ref: "hg" leaves reference contig names as their UCSC-style name (e.g. "chr1"), and "grc" translates reference
#   contigs to their GRC name (e.g. "1").
# wildcard fmt: "sv" writes SV annotations using the SVTYPE and SVLEN INFO fields, and "alt" encodes SVs in
#   "REF" and "ALT" fields.
rule variant_caller_vcf_create_vcf:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/bed/{sample}/{varset}/{filter}/byref/{vartype}_{svtype}.bed',
        fa=lambda wildcards: 'results/variant/{sourcetype}/{sourcename}/fasta/{sample}/{varset}/{filter}/{vartype}_{svtype}.fa.gz'.format(**wildcards) if wildcards.fmt == 'alt' else '',
        ref=config['reference']
    output:
        vcf=temp('temp/variant/{sourcetype}/{sourcename}/vcf/{sample}/{varset}/{filter}/{ref,grc|hg}_{fmt,sv|alt}/{vartype}_{svtype}.vcf')
    params:
        cores=4
    run:

        # Read variants
        df = pd.read_csv(input.bed, sep='\t', header=0)

        df.columns = [col if col != '#CHROM' else 'CHROM' for col in df.columns]  # Remove comment from column name

        df.sort_values(['CHROM', 'POS', 'SVTYPE'], inplace=True)

        # Check for acceptable variant types
        svtype_set = set(df['SVTYPE'])
        unsupported_set = svtype_set - {'INS', 'DEL', 'INV'}

        if unsupported_set:
            raise RuntimeError(
                'Unsupported SVTYPE values in variant BED: {} (supported: INS, DEL, INV)'.format(
                    ', '.join(sorted(unsupported_set))
                )
            )


        # Note: POS already to 1-based and the base before (same as BED with 0-based)

        # Determine if contig fields are present
        has_contig = np.all([val in df.columns for val in ('CONTIG', 'CONTIG_START', 'CONTIG_END')])
        has_contig_support = np.all([val in df.columns for val in ('CONTIG_SUPPORT', 'CONTIG_DEPTH')])

        # Set GT
        if 'GT' not in df:
            df['GT'] = '1/.'
        else:
            df['GT'].fillna('./.', inplace=True)

        # Set FILTER
        if 'FILTER' not in df:
            df['FILTER'] = '.'
        else:
            df['FILTER'].fillna('.', inplace=True)  # Note: Should set FILTER headers and check that filter values match

        # Set QUAL
        if 'QUAL' not in df:
            df['QUAL'] = '.'
        else:
            df['QUAL'] = [str(val.strip()) for val in df['QUAL'].fillna('.')]

        # Set SVLEN
        df['SVLEN'] = df.apply(lambda row: abs(row['SVLEN']) * (-1 if row['SVTYPE'] == 'DEL' else 1), axis=1)

        # Set REF to the base before the SV
        df = analib.pd.apply_parallel(df, _variant_vcf_ref_base, params.cores, kwds={'ref_file_name': input.ref})

        # Set index
        df.set_index('ID', inplace=True, drop=False)

        # Set ALT
        if wildcards.fmt == 'sv':
            df['ALT'] = df['SVTYPE'].apply(lambda val: '<{}>'.format(val))

        elif wildcards.fmt == 'alt':

            # Cannot run for inversions
            if 'INV' in svtype_set:
                raise RuntimeError('Cannot create a VCF in ALT format for variants with INV (need to use symbolic SV format, see "fmt" wildcard)')

            # Save reference base
            df['REF_BASE'] = df['REF']

            # Set SEQ
            df['SEQ'] = analib.seq.fa_to_series(input.fa)

            missing_list = set(df.loc[pd.isnull(df['SEQ']), 'ID'])

            if missing_list:
                raise RuntimeError('Missing sequence for {} ID(s): {}{}'.format(
                    len(missing_list),
                    ', '.join(sorted(missing_list)[0:3]),
                    '...' if len(missing_list) > 3 else ''
                ))

            # Set REF/ALT
            df['REF'] = df.apply(lambda row: (row['REF_BASE'] + (row['SEQ'] if row['SVTYPE'] == 'DEL' else '').upper()), axis=1)
            df['ALT'] = df.apply(lambda row: (row['REF_BASE'] + (row['SEQ'] if row['SVTYPE'] == 'INS' else '').upper()), axis=1)

        else:
            raise RuntimeError('Unknown format (fmt wildcard): {}'.format(wildcards.fmt))


        # Translate to GRC contig names
        if wildcards.ref == 'grc':
            df['CHROM'] = analib.ref.grc_to_hg_chrom(df['CHROM'], 'GRCh38', rev=True)


        # Get GT and sample fields
        if has_contig_support:
            df['FORMAT'] = 'GT:CS:CD'
            df['VCF_SAMPLE'] = df.apply(lambda row: '{GT}:{CONTIG_SUPPORT}:{CONTIG_DEPTH}'.format(**row), axis=1)
        else:
            df['FORMAT'] = 'GT'
            df['VCF_SAMPLE'] = df.apply(lambda row: '{GT}'.format(**row), axis=1)

        # Set INFO
        info_format = (
            """SVTYPE={SVTYPE};"""
            """SVLEN={SVLEN};"""
            """END={END};"""
        )

        if has_contig_support:
            info_format += (
                """CONTIG_SUPPORT={CONTIG_SUPPORT};"""
                """CONTIG_DEPTH={CONTIG_DEPTH};"""
            )

        if has_contig:
            info_format += (
                """CONTIG={CONTIG};"""
                """CONTIG_START={CONTIG_START};"""
                """CONTIG_END={CONTIG_END};"""
            )

        info_format += """BKPTID={ID}"""

        df['INFO'] = df.apply(lambda row: info_format.format(**row), axis=1)


        # Open output file
        with open(output.vcf, 'w') as out_file:
            # Inital metadata
            out_file.write("""##fileformat=VCFv4.2\n""")
            out_file.write("""##fileDate={}\n""".format(time.strftime('%Y%m%d')))
            out_file.write("""##source=svpop-1.0\n""")
            out_file.write("""##reference=file://{}\n""".format(config['reference']))

            # Contig metadata
            with open(config['reference'] + '.fai', 'r') as fai_file:
                for line in fai_file:
                    if not line:
                        continue

                    tok = line.split('\t')

                    out_file.write("""##contig=<ID={0},length={1}>\n""".format(*tok))

            # ALT tags

            if 'INS' in svtype_set:
                out_file.write("""##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">\n""")

            if 'DEL' in svtype_set:
                out_file.write("""##ALT=<ID=DEL,Description="Deletion relative to the reference">\n""")

            if 'INV' in svtype_set:
                out_file.write("""##ALT=<ID=INV,Description="Inversion of reference sequence">\n""")

            # Info field metadata
            out_file.write("""##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n""")
            out_file.write("""##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between the REF and ALT alleles">\n""")
            out_file.write("""##INFO=<ID=END,Number=1,Type=Integer,Description="Reference end position of the variant">\n""")

            if has_contig_support:
                out_file.write("""##INFO=<ID=CONTIG_SUPPORT,Number=1,Type=Integer,Description="Number of contigs supporting this variant">\n""")
                out_file.write("""##INFO=<ID=CONTIG_DEPTH,Number=1,Type=Integer,Description="Number of contigs at the variant locus">\n""")

            if has_contig:
                out_file.write("""##INFO=<ID=CONTIG,Number=1,Type=String,Description="Contig name variant was called from">\n""")
                out_file.write("""##INFO=<ID=CONTIG_START,Number=1,Type=Integer,Description="Start position in the contig for this variant">\n""")
                out_file.write("""##INFO=<ID=CONTIG_END,Number=1,Type=Integer,Description="End position in the contig for this variant">\n""")

            out_file.write("""##INFO=<ID=BKPTID,Number=1,Type=String,Description="ID of the assembled alternate allele in the assembly file">\n""")

            # Format metadata
            out_file.write("""##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n""")

            if has_contig_support:
                out_file.write("""##FORMAT=<ID=CS,Number=1,Type=Integer,Description="Number of local assembly contigs supporting this variant">\n""")
                out_file.write("""##FORMAT=<ID=CD,Number=1,Type=Integer,Description="Number of local assembly contigs at this location">\n""")

            # NOTE: Missing FILTER header

            # Column labels
            out_file.write("""#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n""".format(wildcards.sample))

            # Write variant records
            for index, row in df.iterrows():

                # Write record
                out_file.write(
                    '\t'.join([str(element) for element in row.loc[
                        ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'VCF_SAMPLE']]
                    ])
                )
                out_file.write('\n')
