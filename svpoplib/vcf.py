"""
VCF writer.
"""

import Bio.SeqIO
import datetime
import numpy as np
import os
import pandas as pd
import pysam

import svpoplib

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

VARIANT_VCF_INFO_HEADER = {
    'ID': ('ID', '1', 'String', 'Variant ID'),
    'SVTYPE': ('SVTYPE', '1', 'String', 'Variant type'),
    'SVLEN': ('SVLEN', '.', 'Integer', 'Variant length'),
    'SEQ': ('SEQ', '.', 'String', 'Variant sequence')
}

VARIANT_VCF_ALT_HEADER = {
    'INS': ('INS', 'Sequence insertion'),
    'DEL': ('DEL', 'Sequence deletion'),
    'INV': ('INV', 'Inversion'),
    'DUP': ('DUP', "Duplication")
}

VARIANT_VCF_SVTYPE_LIST = {
    'insdel': ('INS', 'DEL'),
    'insdelinv': ('INS', 'DEL', 'INV')
}


def header_list(
        df_ref,
        info_fields,
        format_fields,
        alt_fields = tuple(),
        filter_fields = tuple(),
        vcf_version='4.2',
        file_date=None,
        variant_source='SV-Pop',
        ref_filename=None
    ):

    """
    Get a list of formatted header lines for a VCF. Each line ends with a newline. The column headings
    ("#CHROM", "POS", "REF"...) is not part of this list.

    For custom headers (tuples in `info_fields`, `format_fields`, `alt_fields`, `filter_fields`), do not include
    quotes in the string around the elements (quotes will be added and written, e.g. Description="...").

    :param df_ref: DataFrame of reference information (See svpoplib.ref.get_ref_info or rule data_ref_contig_table).
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
    :param ref_filename: Name of the reference file. If `None` or empty, then no "reference" header is written.

    :return: A list of header lines including newlines for each line.
    """

    header_list = list()

    if ref_filename is not None and not ref_filename.strip():
        ref_filename = None
    else:
        ref_filename = ref_filename.strip()

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

    if ref_filename:
        header_list.append("""##reference=file://{}\n""".format(os.path.basename(ref_filename)))

    # Reference contigs
    if df_ref is not None:
        for index, row in df_ref.reset_index().iterrows():
            header_list.append("""##contig=<ID={CHROM},length={LEN},md5={MD5}>\n""".format(**row))

    else:
        if ref_filename is None:
            raise RuntimeError('Either a reference table (CHROM, LEN, MD5 columns) or a reference FASTA file is required to write VCF headers.')

        ref_fai_name = ref_filename + '.fai'

        if not os.path.isfile(ref_fai_name):
            raise RuntimeError(f'Reference FAI is missing: {ref_fai_name}')

        df_fai = svpoplib.ref.get_df_fai(ref_fai_name)

        for index, row in df_fai.reset_index().iterrows():
            header_list.append("""##contig=<ID={CHROM},length={LEN}>\n""".format(**row))

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
    Get reference base preceding a variant (SV, indel) or at the point change (SNV).

    :param df: Variant dataframe as BED.
    :param ref_fa: Reference file.
    """

    # Open and update records
    with pysam.FastaFile(ref_fa) as ref_file:
        for index, row in df.iterrows():
            if row['SVTYPE'] in {'INS', 'DEL', 'INSDEL', 'DUP'}:
                yield ref_file.fetch(row['#CHROM'], row['POS'] + (-1 if row['POS'] > 0 else 0), row['POS']).upper()

            elif row['SVTYPE'] == 'SNV':
                if 'REF' in row:
                    yield row['REF']
                else:
                    yield ref_file.fetch(row['#CHROM'], row['POS'], row['POS'] + 1).upper()

            else:
                raise RuntimeError('Unknown variant type: "{}" at index {}'.format(row['VARTYPE'], index))


class VariantVcfTable:
    def __init__(self, df, sample, ref_filename, fa_filename=None, altfmt='alt'):

        self.sample = sample
        self.ref_filename = ref_filename
        self.fa_filename = fa_filename
        self.altfmt = altfmt

        # Check assembly name
        if sample in {'#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'}:
            raise RuntimeError(f'Sample name conflicts with a VCF header column name: {sample}')

        # Check alt format
        if altfmt == 'alt':
            symbolic_alt = False
            alt_seq = True

        elif altfmt == 'sym':
            symbolic_alt = True
            alt_seq = True

        elif altfmt == 'sym-noseq':
            symbolic_alt = True
            alt_seq = True

        else:
            raise RuntimeError(f'Unknown alt format wildcard (altfmt): {altfmt}')

        # Get SVTYPES
        svtype_set = set(df['SVTYPE'])

        # Check for symbolic ALT compatibility
        if symbolic_alt:
            unsupported = svtype_set - {'INS', 'DEL', 'INV', 'DUP'}

            if unsupported:
                n = len(unsupported)
                unsupported = ', '.join(sorted(unsupported)[:3]) + '..' if n > 3 else ''
                raise RuntimeError(f'Symbolic ALT output is not supported for input variant type(s): {unsupported}')

            is_snv = False

        else:
            unsupported = svtype_set - {'INS', 'DEL', 'SNV'}

            if unsupported:
                n = len(unsupported)
                unsupported = ', '.join(sorted(unsupported)[:3]) + '..' if n > 3 else ''
                raise RuntimeError(f'REF/ALT output is not supported for input variant type(s): {unsupported}')

            if 'SNV' in svtype_set:
                if svtype_set != {'SNV'}:
                    raise RuntimeError(f'Mixing SNVs and other variant types is not supported: {", ".join(sorted(svtype_set))}')

                is_snv = True
            else:
                is_snv = False

        # Setup header lists
        info_header_id_list = list()

        # Read sequence from FASTA
        read_seq = df.shape[0] > 0 and alt_seq and (
            'SEQ' not in df.columns or np.any(pd.isnull(df['SEQ']))
        )

        if read_seq:

            # FASTA must have been given
            if fa_filename is None:
                raise RuntimeError(f'Missing sequence data to construct VCF: No variant FASTA given and no SEQ column in variant table.')

            # FASTA cannot be empty
            if not os.path.isfile(fa_filename) or os.stat(fa_filename).st_size == 0:
                raise RuntimeError(f'Missing sequence data to construct VCF: Variant FASTA is missing or empty: {fa_filename}')

            # Get Sequence from FASTA and assign to SEQ column (not for SNVs)
            with svpoplib.seq.PlainOrGzReader(fa_filename, 'rt') as fa_in:
                df_seq_dict = {
                    record.name: str(record.seq) for record in Bio.SeqIO.parse(fa_in, 'fasta')
                }

            # Assign SEQ to records
            df.set_index('ID', inplace=True)

            df['SEQ'] = pd.Series(df_seq_dict)

            del(df_seq_dict)

            df.reset_index(inplace=True)

            # Check for missing sequences
            df_null_seq = df.loc[pd.isnull(df['SEQ'])]

            if df_null_seq.shape[0] > 0:
                raise RuntimeError(
                    'Missing FASTA sequence for {} variants: {}{}'.format(
                        df_null_seq.shape[0],
                        ', '.join([str(val) for val in df_null_seq.iloc[:3]['ID']]),
                        '...' if df_null_seq.shape[0] > 3 else ''
                    )
                )

        # # Add VARTYPE
        # df['VARTYPE'] = wildcards.vartype.upper()

        # SVTYPE
        if df.shape[0] > 0:
            df['SVTYPE'] = df['SVTYPE'].apply(lambda val: val.upper())

        # Reformat fields for INFO
        df['SVLEN'] = df.apply(lambda row: np.abs(row['SVLEN']) * (-1 if row['SVTYPE'] == 'DEL' else 1), axis=1)

        # Add GT if missing
        if 'GT' not in df.columns:
            df['GT'] = '1/.'

        # INFO: Base
        df['INFO'] = df.apply(lambda row: 'ID={ID};SVTYPE={SVTYPE}'.format(**row), axis=1)

        info_header_id_list.append('ID')
        info_header_id_list.append('SVTYPE')

        # INFO: Add SV/INDEL annotations
        if df.shape[0] > 0:
            df['INFO'] = df.apply(lambda row: row['INFO'] + (';SVLEN={SVLEN}'.format(**row)) if row['SVTYPE'] != 'SNV' else '', axis=1)

        info_header_id_list.append('SVLEN')

        # REF
        df['REF'] = list(svpoplib.vcf.ref_base(df, ref_filename))

        # ALT
        if not is_snv:
            if df.shape[0] > 0:
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
                df['ALT'] = np.nan

        else:
            # Fix position for SNVs (0-based BED to 1-based VCF)
            df['POS'] += 1

        # Remove SEQ
        if alt_seq and df.shape[0] > 0:
            del df['SEQ']

        # No masked bases in REF/ALT
        df['REF'] = df['REF'].apply(lambda val: val.upper())
        df['ALT'] = df['ALT'].apply(lambda val: val.upper())

        # Save columns needed for VCF
        df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO', 'GT']].copy()

        # Sort
        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # INFO headers
        info_header_list = [svpoplib.vcf.VARIANT_VCF_INFO_HEADER[val] for val in info_header_id_list]

        # ALT headers
        if symbolic_alt:
            alt_header_list = [svpoplib.vcf.VARIANT_VCF_ALT_HEADER[svtype] for svtype in sorted(svtype_set)]
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
        df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample]


        # Assign fields
        self.df = df
        self.info_header_list = info_header_list
        self.format_header_list = format_header_list
        self.alt_header_list = alt_header_list
        self.filter_header_list = filter_header_list

    def write(self, out_file, df_ref=None):
        for line in header_list(
            df_ref=df_ref,
            info_fields=self.info_header_list,
            format_fields=self.format_header_list,
            alt_fields=self.alt_header_list,
            filter_fields=self.filter_header_list,
            variant_source=f'SV-Pop {svpoplib.constants.VERSION}',
            ref_filename=self.ref_filename
        ):
            out_file.write(line)

        out_file.write('\t'.join(self.df.columns))
        out_file.write('\n')

        for index, row in self.df.iterrows():
            out_file.write('\t'.join(row.astype(str)))
            out_file.write('\n')