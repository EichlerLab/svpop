# Make dup calls based on re-mapping INS sequences back to the reference. Useful for identifying the likely
# templated location and for matching INS calls with DUP calls.

def variant_bed_altdup_input_filename(wildcards, type):

    # Get BED (altmap input) or FA (original FASTA file)
    TYPE_PATTERN = {
        'altmap': 'results/variant/{sourcetype}/{sourcename}/{{sample}}/{filter}/all/anno/altmap/altmap-{aligner}_sv_ins.bed.gz',
        'bed': 'results/variant/{sourcetype}/{sourcename}/{{sample}}/{filter}/all/bed/sv_ins.bed.gz',
        'fa': 'results/variant/{sourcetype}/{sourcename}/{{sample}}/{filter}/all/bed/fa/sv_ins.fa.gz'
    }

    # Check type
    if type not in TYPE_PATTERN.keys():
        raise RuntimeError('Error getting input file for altdup input: Requested type "{type}" is unknown, should be "bed" or "fa"')

    # Get sample entry
    sample_entry = svpoplib.rules.sample_table_entry(
        wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='altdup'
    )

    # Get DATA tokens (sourcetype, sourcename, alt-map aligner)
    if (pd.isnull(sample_entry['DATA'])):
        raise RuntimeError('Alt-dup sample entry {} is missing the "DATA" field: Expected "sourcetype,sourcename,aligner"'.format(sample_entry['NAME']))

    data_tok = [val.strip() for val in sample_entry['DATA'].strip().split(',')]

    if len(data_tok) not in (3, 4) or np.any([val == '' for val in data_tok[:3]]):
        raise RuntimeError('Alt-dup sample entry {} has an unrecognized "DATA" field: Expected "sourcetype,sourcename,aligner (optional filter after aligner)": Found "{}"'.format(sample_entry['NAME'], sample_entry['DATA']))

    if len(data_tok) == 3:
        data_tok += ['all']
    elif data_tok[3] == '':
        data_tok[3] = 'all'

    data_tok_dict = {key: val for key, val in zip(('sourcetype', 'sourcename', 'aligner', 'filter'), data_tok)}

    # Return file
    return TYPE_PATTERN[type].format(**data_tok_dict)


def dup_seq_iter(fa_file_name, id_set, df_altmap):
    """
    Iterate over the insertion FASTA file (sequences mapped to produce DUP calls) and produce an iterator for the
    mapped portion of insertion sequences. If bases were clipped from either end of the insertion when aligning,
    those bases are removed by this function (i.e. only aligned portion of the original insertion is reported).
    Sequences are output in reference orientation.

    :param fa_file_name: Insertion FASTA file name.
    :param id_set: Set of insertion IDs to include.
    :param df_altmap: The altmap BED file after initial processing (ID is the duplication ID, INS_ID is the ID of
        the original insertion).

    :return: An iterator of Bio.SeqRecord.SeqRecord objects. May be fed directly to Bio.SeqIO.write.
    """

    df_align = df_altmap[['INS_ID', 'ID', 'IS_REV', 'QUERY_POS', 'QUERY_END']].set_index('INS_ID')

    for record in svpoplib.seq.fa_to_record_iter(fa_file_name, id_set):
        row_align = df_align.loc[record.id]

        seq = record.seq

        if row_align['IS_REV']:
            seq = seq.reverse_complement()

        seq = Bio.Seq.Seq(str(seq)[row_align['QUERY_POS']:row_align['QUERY_END']])

        yield Bio.SeqRecord.SeqRecord(seq, id=row_align['ID'], name='', description='')


#
# Rules
#

# variant_bed_altdup_inner
#
# Call inner-variants from duplicated regions using the duplicated copy aligned to the reference.
rule variant_bed_altdup_inner:
    input:
        altmap='results/variant/caller/{sourcename}/{sample}/all/all/altdup/align.bed.gz',
        fa=lambda wildcards: variant_bed_altdup_input_filename(wildcards, 'fa')
    output:
        bed_sv_ins=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/sv_ins.bed.gz'),
        fa_sv_ins=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/fa/sv_ins.fa.gz'),
        bed_sv_del=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/sv_del.bed.gz'),
        fa_sv_del=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/fa/sv_del.fa.gz'),
        bed_indel_ins=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/indel_ins.bed.gz'),
        fa_indel_ins=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/fa/indel_ins.fa.gz'),
        bed_indel_del=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/indel_del.bed.gz'),
        fa_indel_del=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/fa/indel_del.fa.gz'),
        bed_snv_snv=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/snv_snv.bed.gz'),
        fa_snv_snv=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/fa/snv_snv.fa.gz')
    run:

        # Read
        df_map = pd.read_csv(input.altmap, sep='\t')

        # Call inner-variants
        df_snv_list = list()
        df_svindel_list = list()

        with pysam.FastaFile(input.fa) as fa_qry, pysam.FastaFile(config['reference']) as fa_sub:

            for index, row in df_map.iterrows():
                pos_sub = row['POS']        # Subject (reference) position
                pos_qry = row['QUERY_POS']  # Query (SV insertion sequence) position

                chrom = row['#CHROM']
                ins_id = row['INS_ID']
                dup_id = row['ID']

                seq_qry = fa_qry.fetch(ins_id).upper()

                if row['IS_REV']:
                    seq_qry = str(Bio.Seq.Seq(seq_qry).reverse_complement())

                for cigar_rec in re.findall(r'\d+[^\d]', row['CIGAR']):
                    cigar_len = int(cigar_rec[:-1])
                    cigar_op = cigar_rec[-1]

                    if cigar_op == '=':  # Match
                        pos_sub += cigar_len
                        pos_qry += cigar_len

                    elif cigar_op == 'X':  # Mismatch
                        for i in range(cigar_len):

                            pos_snv_sub = pos_sub + i
                            pos_snv_qry = pos_qry + i

                            df_snv_list.append(
                                pd.Series(
                                    [
                                        chrom, pos_snv_sub, pos_snv_sub + 1,
                                        fa_sub.fetch(chrom, pos_snv_sub, pos_snv_sub + 1).upper(),
                                        seq_qry[pos_snv_qry:pos_snv_qry + 1],
                                        ins_id, dup_id
                                    ],
                                    index=['#CHROM', 'POS', 'END', 'REF', 'ALT', 'INS_ID', 'DUP_ID']
                                )
                            )

                        pos_sub += cigar_len
                        pos_qry += cigar_len

                    elif cigar_op == 'I':  # Insertion
                        df_svindel_list.append(
                            pd.Series(
                                [
                                    chrom, pos_sub, pos_sub + 1,
                                    'INS', cigar_len,
                                    seq_qry[pos_qry:pos_qry + cigar_len],
                                    ins_id, dup_id
                                ],
                                index=['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'SEQ', 'INS_ID', 'DUP_ID']
                            )
                        )

                        pos_qry += cigar_len

                    elif cigar_op == 'D':  # Deletion
                        df_svindel_list.append(
                            pd.Series(
                                [
                                    chrom, pos_sub, pos_sub + cigar_len,
                                    'DEL', cigar_len,
                                    fa_sub.fetch(chrom, pos_sub, pos_sub + cigar_len).upper(),
                                    ins_id, dup_id
                                ],
                                index=['#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'SEQ', 'INS_ID', 'DUP_ID']
                            )
                        )

                        pos_sub += cigar_len

                    elif cigar_op in {'S', 'H'}:  # Soft and hard clipping
                        continue  # Position already accounted for

                    elif cigar_op == 'M':  # Aligned bases (match & mismatch)
                        for i in range(cigar_len):

                            pos_snv_sub = pos_sub + i
                            pos_snv_qry = pos_qry + i

                            base_sub = fa_sub.fetch(chrom, pos_snv_sub, pos_snv_sub + 1).upper()
                            base_qry = seq_qry[pos_snv_qry:pos_snv_qry + 1]

                            if base_sub != base_qry:
                                df_snv_list.append(
                                    pd.Series(
                                        [
                                            chrom, pos_snv_sub, pos_snv_sub + 1,
                                            base_sub,
                                            base_qry,
                                            ins_id, dup_id
                                        ],
                                        index=['#CHROM', 'POS', 'END', 'REF', 'ALT', 'INS_ID', 'DUP_ID']
                                    )
                                )

                        pos_sub += cigar_len
                        pos_qry += cigar_len

                    elif cigar_op == 'P':  # Padding
                        continue

                    elif cigar_op == 'N':  # Skipped subject bases
                        pos_sub += cigar_len

                    else:
                        raise RuntimeError(f'Unrecognized CIGAR operation for record {dup_id}: {cigar_op}')

        # Merge and finish tables
        df_svindel = pd.concat(df_svindel_list, axis=1).T.sort_values(['#CHROM', 'POS', 'END', 'SVLEN'])
        del(df_svindel_list)

        df_snv = pd.concat(df_snv_list, axis=1).T.sort_values(['#CHROM', 'POS', 'END'])
        del(df_snv_list)

        df_snv['SVTYPE'] = 'SNV'
        df_snv['SVLEN'] = 1

        df_svindel['ID'] = svpoplib.variant.get_variant_id(df_svindel)
        df_snv['ID'] = svpoplib.variant.get_variant_id(df_snv)

        df_svindel = svpoplib.variant.order_variant_columns(df_svindel)
        df_snv = svpoplib.variant.order_variant_columns(df_snv)

        # Write SV & indel
        for vartype in ('sv', 'indel'):

            # Subset by vartype
            if vartype == 'sv':
                subdf_vartype = df_svindel.loc[df_svindel['SVLEN'] >= 50]
            elif vartype == 'indel':
                subdf_vartype = df_svindel.loc[df_svindel['SVLEN'] < 50]
            else:
                raise RuntimeError(f'Illegal variant type: {vartype}')

            for svtype in ('ins', 'del'):

                # Subset by svtype
                subdf = subdf_vartype.loc[subdf_vartype['SVTYPE'] == svtype.upper()].copy()

                # Write FA
                with Bio.bgzf.BgzfWriter(output[f'fa_{vartype}_{svtype}'], 'wt') as fa_file:
                    Bio.SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(subdf), fa_file, 'fasta')

                # Write table
                del(subdf['SEQ'])

                subdf.to_csv(output[f'bed_{vartype}_{svtype}'], sep='\t', index=False, compression='gzip')

        # Write SNV
        df_snv.to_csv(output.bed_snv_snv, sep='\t', index=False, compression='gzip')

        with Bio.bgzf.BgzfWriter(output[f'fa_snv_snv'], 'wt') as fa_file:
            Bio.SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(df_snv, seq_col='ALT'), fa_file, 'fasta')


# variant_bed_altdup_bed
#
# Get ALT-DUP BED file from INS.
rule variant_bed_altdup_bed:
    input:
        altmap=lambda wildcards: variant_bed_altdup_input_filename(wildcards, 'altmap'),
        bed=lambda wildcards: variant_bed_altdup_input_filename(wildcards, 'bed'),
        fa=lambda wildcards: variant_bed_altdup_input_filename(wildcards, 'fa')
    output:
        bed=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/sv_dup.bed.gz'),
        fa=temp('temp/variant/caller/altdup/{sourcename}/{sample}/all/all/bed/pre_filter/fa/sv_dup.fa.gz'),
        altmap='results/variant/caller/{sourcename}/{sample}/all/all/altdup/align.bed.gz'
    run:

        # Read altmap and transform
        df = pd.read_csv(input.altmap, sep='\t')

        df['SVTYPE'] = 'DUP'
        df['SVLEN'] = df['MAP_LEN']

        df['INS_ID'] = df['ID']

        df['ID'] = svpoplib.variant.get_variant_id(df)

        df = svpoplib.variant.order_variant_columns(df)

        del(df['MAP_LEN'])

        df.sort_values(['#CHROM', 'POS', 'END', 'SVLEN'], inplace=True)

        df_altmap = df.copy()

        # Clear unecessary altmap fields
        del(df['CIGAR'])
        del(df['MATCH_BP'])
        del(df['MISMATCH_BP'])
        del(df['INS_BP'])
        del(df['INS_N'])
        del(df['DEL_BP'])
        del(df['DEL_N'])
        del(df['CLIPPED_BP'])
        del(df['CLIPPED_N'])

        # Filter for SV-sized remapping
        df = df.loc[df['SVLEN'] >= 50].copy()

        # Add BED columns
        df.index = list(df['INS_ID'])  # Do not use set_index, it will change the index when ID is reassigned
        df.index.name = 'INDEX'
        del(df['INS_ID'])  # Gets merged back in from df_bed

        df_bed = pd.read_csv(input.bed, sep='\t', usecols=('#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN'))
        df_bed.set_index('ID', inplace=True, drop=False)

        df_bed.columns = [f'INS_{col}' for col in df_bed.columns]

        df_bed = df_bed.reindex(df.index)

        df = pd.concat([df, df_bed], axis=1)

        del(df_bed)

        # Write tables
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')
        df_altmap.to_csv(output.altmap, sep='\t', index=False, compression='gzip')

        # Write FASTA
        ins_id_set = set(df['INS_ID'])

        with Bio.bgzf.BgzfWriter(output.fa) as out_file:
            Bio.SeqIO.write(dup_seq_iter(input.fa, set(df['INS_ID']), df_altmap), out_file, 'fasta')
