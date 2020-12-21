"""
VCF writer.
"""

import datetime
import os
import pysam

INFO_FIELDS = [
    ('VARTYPE', 'A', 'String', 'Variant class (SV, INDEL, SNV)'),
    ('SVTYPE', 'A', 'String', 'Variant type (INS, DEL, INV, SNV)'),
    ('SVLEN', '.', 'Integer', 'Difference in length of ref and alt alleles'),
    ('END', 'A', 'Integer', 'Variant length (1 for SNVs, negative for DELs)'),
    ('LEAD_SAMPLE', 'A', 'String', 'Merged VCFs use the variant ID from the leading sample (first sample this variant was found it)'),
    ('ID', 'A', 'String', 'Representative variant for this record was taken from the lead sample (source of contig and sequence information)'),
    ('TIG_REGION', 'A', 'String', 'Assembly contig and coordinates (1-based, closed).'),
    ('TIG_STRAND', 'A', 'String', 'Assembly contig strand relative to the reference (+ or -).'),
    ('SEQ', 'A', 'String', 'Variant sequence in reference orientation.')
]

FORMAT_FIELDS = [
    ('GT', '1', 'String', 'Genotype')
]

ALT_FIELDS = [
    ('INS', 'Sequence insertion'),
    ('INS:ME', 'Mobile element insertion'),
    ('INS:ME:L1', 'L1 mobile element insertion'),
    ('INS:ME:ALU', 'Alu mobile element insertion'),
    ('DEL', 'Sequence deletion'),
    ('DEL:ME', 'Mobile element deletion'),
    ('DEL:ME:L1', 'L1 mobile element deletion'),
    ('DEL:ME:ALU', 'Alu mobile element deletion'),
    ('INV', 'Inversion'),
    ('DUP', 'Sequence duplication'),
    ('DUP:TANDEM', 'Tandem duplication'),
    ('CNV', 'Copy number variable region')
]

FILTER_FIELDS = [
    ('PASS', 'Passed all QV stages')
]

def header_list(
        df_ref,
        info_fields,
        format_fields,
        alt_fields = [],
        filter_fields = [],
        vcf_version='4.2',
        file_date=None,
        variant_source='SV-Pop',
        ref_file_name=None
    ):

    """
    Get a list of formatted header lines for a VCF. Each line ends with a newline. The column headings
    ("#CHROM", "POS", "REF"...) is not part of this list.

    For custom headers (tuples in `info_fields`, `format_fields`, `alt_fields`, `filter_fields`), do not include
    quotes in the string around the elements (quotes will be added and written, e.g. Description="...").

    :param df_ref: DataFrame of reference information (See analib.ref.get_ref_info or rule data_ref_contig_table).
    :param info_fields: List of INFO headers. May be 4-element tuples (ID, Number, Type, Description), built-in field
        names ("VARTYPE", "SVTYPE", "END", "ID", etc), or a mix of the two.
    :param format_fields: List of FORMAT headers. May be 4-element tuples (ID, Number, Type, Description), built-in field
        names ("GT"), or a mix of the two.
    :param alt_fields: List of ALT headers. May be 2-element tuples (ID, Description), built-in field
        names ("INS", "DEL", "INV", "INS:ME", etc), or a mix of the two.
    :param filter_fields: List of FILTER headers. May be 2-element tuples (ID, Description), built-in field
        names ("PASS", etc), or a mix of the two.
    :param vcf_version: VCF format version (should comply with VCF specs).
    :param file_date: Date string or `None` to generate from current date.
    :param variant_source: Variant caller or variant source.
    :param ref_file_name: Name of the reference file. If `None` or empty, then no "reference" header is written.

    :return: A list of header lines including newlines for each line.
    """

    header_list = list()

    # Get dict of known fields
    info_field_dict = {
        tup[0]: tup for tup in INFO_FIELDS
    }

    format_field_dict = {
        tup[0]: tup for tup in FORMAT_FIELDS
    }

    alt_field_dict = {
        tup[0]: tup for tup in ALT_FIELDS
    }

    filter_field_dict = {
        tup[0]: tup for tup in FILTER_FIELDS
    }

    # Base headers (format, date, source, reference file)
    if not vcf_version and not vcf_version.strip():
        raise RuntimeError('VCF version may not be emtpy or missing')

    header_list.append(f"""##fileformat=VCFv{vcf_version}\n""")

    if file_date != False:
        if file_date is None:
            timestamp = datetime.datetime.now()
            file_date = timestamp.strftime('%Y%m%d')

        header_list.append(f"""##fileDate={file_date}\n""")

    if variant_source and variant_source.strip():
        header_list.append(f"""##source={variant_source}\n""")

    if ref_file_name and ref_file_name.strip():
        ref_file_name = ref_file_name.strip()

        if '://' not in ref_file_name:
            header_list.append("""##reference=file://{}\n""".format(os.path.basename(ref_file_name.strip())))
        else:
            header_list.append("""##reference={}\n""".format(ref_file_name.strip()))

    # Reference contigs
    for index, row in df_ref.reset_index().iterrows():
        header_list.append("""##contig=<ID={CHROM},length={LEN},md5={MD5}>\n""".format(**row))

    # INFO
    info_element_set = set()

    for info_element in info_fields:
        if issubclass(info_element.__class__, tuple):
            if len(info_element) != 4:
                raise RuntimeError('INFO header element {} is not a 4-tuple'.format(info_element[0] if len(info_element) > 0 else '<EMPTY TUPLE>'))
        else:
            if info_element not in info_field_dict:
                raise RuntimeError('No built-in INFO header element {}'.format(info_element))

            info_element = info_field_dict[info_element]

        if info_element[0] in info_element_set:
            raise RuntimeError('Duplicate INFO header: {}'.format(info_element[0]))

        info_element_set.add(info_element[0])

        header_list.append("""##INFO=<ID={0},Number={1},Type={2},Description="{3}">\n""".format(*info_element))

    # FORMAT
    format_element_set = set()

    for format_element in format_fields:
        if issubclass(format_element.__class__, tuple):
            if len(format_element) != 4:
                raise RuntimeError('FORMAT header element {} is not a 4-tuple'.format(format_element[0] if len(format_element) > 0 else '<EMPTY TUPLE>'))
        else:
            if format_element not in format_field_dict:
                raise RuntimeError('No built-in FORMAT header element {}'.format(format_element))

            format_element = format_field_dict[format_element]

        if format_element[0] in format_element_set:
            raise RuntimeError('Duplicate FORMAT header: {}'.format(format_element[0]))

        format_element_set.add(format_element[0])

        header_list.append("""##FORMAT=<ID={0},Number={1},Type={2},Description="{3}">\n""".format(*format_element))

    # ALT
    alt_element_set = set()

    for alt_element in alt_fields:
        if issubclass(alt_element.__class__, tuple):
            if len(alt_element) != 2:
                raise RuntimeError('ALT header element {} is not a 2-tuple'.format(alt_element[0] if len(alt_element) > 0 else '<EMPTY TUPLE>'))
        else:
            if alt_element not in alt_field_dict:
                raise RuntimeError('No built-in ALT header element {}'.format(alt_element))

            alt_element = alt_field_dict[alt_element]

        if alt_element[0] in alt_element_set:
            raise RuntimeError('Duplicate FORMAT header: {}'.format(alt_element[0]))

        alt_element_set.add(alt_element[0])

        header_list.append("""##ALT=<ID={0},Description="{1}">\n""".format(*alt_element))

    # FILTER
    filter_element_set = set()

    for filter_element in filter_fields:
        if issubclass(filter_element.__class__, tuple):
            if len(filter_element) != 2:
                raise RuntimeError('ALT header element {} is not a 2-tuple'.format(filter_element[0] if len(filter_element) > 0 else '<EMPTY TUPLE>'))
        else:
            if filter_element not in filter_field_dict:
                raise RuntimeError('No built-in FILTER header element {}'.format(filter_element))

            filter_element = filter_field_dict[filter_element]

        if filter_element[0] in filter_element_set:
            raise RuntimeError('Duplicate FORMAT header: {}'.format(filter_element[0]))

        filter_element_set.add(filter_element[0])

        header_list.append("""##FILTER=<ID={0},Description="{1}">\n""".format(*filter_element))

    # Return headers
    return header_list


def ref_base(df, ref_fa):
    """
    Get reference base.

    :param df: Variant dataframe as BED.
    :param ref_file_name: Reference file.
    """

    # Open and update records
    with pysam.FastaFile(ref_fa) as ref_file:
        for index, row in df.iterrows():
            if row['VARTYPE'] in {'SV', 'INDEL'}:
                yield ref_file.fetch(row['#CHROM'], row['POS'] - 1, row['POS']).upper()

            elif row['VARTYPE'] == 'SNV':
                if 'REF' in row:
                    yield row['REF']
                else:
                    yield ref_file.fetch(row['#CHROM'], row['POS'], row['POS'] + 1).upper()

            else:
                raise RuntimeError('Unknown variant type: "{}" at index {}'.format(row['VARTYPE'], index))

def get_variant_seq(
    df, svpop_run_dir, fa_pattern, var_sv_type, fmt={}
):
    """
    Get fields from pre-merged variants for a set of samples.

    :param df: Merged data frame with variant ID ("ID" column) and discovery sample ("MERGE_SRC" column).
    :param svpop_run_dir: Directory SV-Pop was run from.
    :param bed_pattern: Pattern of BED files to search.
    :param var_sv_type: Variant/SV type keyword to `multivcflib.bed.VAR_SV_TYPE_LIST` to get a list of vartype, svtype
        tuples for all variant types (sv, indel, snv) and svtype (ins, del, snv, inv) that should be read and filled
        into "vartype" and "svtype" wildcards in `bed_pattern`.
    :param fmt: Additional format fields for `bed_pattern` ("sample", "vartype", and "svtype" are filled in
        for each sample and variant BED file).

    :return: A series keyed by IDs for sequences. SNVs and INVs are given np.nan sequence.
    """

    if var_sv_type not in multivcflib.bed.VAR_SV_TYPE_LIST.keys():
        raise RuntimeError('var_sv_type is not a key in VAR_SV_TYPE_LIST: {}'.format(var_sv_type))

    fmt = fmt.copy()

    # Process each sample
    df_list = list()

    # Collect dataframes for each vartype/svtype
    df_seq_list = list()

    for sample in set(df['MERGE_SRC']):

        fmt['sample'] = sample

        for vartype, svtype in multivcflib.bed.VAR_SV_TYPE_LIST[var_sv_type]:

            id_set = set(df.loc[
                (df['VARTYPE'] == vartype.upper()) &
                (df['SVTYPE'] == svtype.upper()) &
                (df['MERGE_SRC'] == sample),
                'ID'
            ])

            if svtype in {'ins', 'del', 'inv'}:

                fmt['vartype'] = vartype
                fmt['svtype'] = svtype

                with gzip.open(os.path.join(svpop_run_dir, fa_pattern.format(**fmt)), 'rt') as in_file:
                    df_seq = pd.Series({record.id: str(record.seq).upper() for record in SeqIO.parse(in_file, 'fasta') if record.id in id_set})

                if df_seq.shape[0] != len(id_set):
                    raise RuntimeError('Mismatch in sample {}: {} in df, {} in FASTA'.format(sample, len(id_set), df_seq.shape[0]))

                df_seq_list.append(df_seq)

            elif svtype in {'snv'}:
                df_seq_list.append(pd.Series(np.nan, index=list(id_set)))

            else:
                raise RuntimeError('Don\'t know how to handle SVTYPE {}'.format(svtype))

    return pd.concat(df_seq_list, axis=0)
