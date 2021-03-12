"""
Global rules for variants.
"""

###################
### Definitions ###
###################

def variant_global_sample_table_entry(sourcename, sample=None, wildcards=None):
    """
    Get an entry from SAMPLE_TABLE.

    :param sourcename: Caller source name.
    :param sample: Sample name. If not present, extract from "sample" entry in `wildcards`.
    :param wildcards: If present, parse the table entry with this value. Otherwise, parse with sample=sample.

    :return: A configuration entry (Pandas.Series).
    """

    if sample is None:
        if wildcards is None:
            raise RuntimeError('Cannot get global sampleset entry: sample and wildcards are both None')

        if 'sample' not in wildcards.keys():
            raise RuntimeError('Cannot get sample from wildcards: "sample" is not a wildcards key')

        sample = wildcards.sample

    if wildcards is not None and 'seq_set' in wildcards.keys():
        seq_set = wildcards.seq_set
    else:
        seq_set = 'DEFAULT'

    if (sourcename, sample, seq_set) in SAMPLE_TABLE.index:
        sample_entry = SAMPLE_TABLE.loc[(sourcename, sample, seq_set)]

    elif (sourcename, 'DEFAULT', seq_set) in SAMPLE_TABLE.index:
        sample_entry = SAMPLE_TABLE.loc[(sourcename, 'DEFAULT', seq_set)]

    else:
        raise RuntimeError('Cannot find sample table entry for sourcename "{}" and sample "{}" (or sample "DEFAULT") with seq-set {}'.format(
            sourcename, sample, ('"' + seq_set + '"') if not pd.isnull(seq_set) else '<Undefined>'
        ))

    if wildcards is not None:
        sample_entry['DATA'] = sample_entry['DATA'].format(**wildcards)
    else:
        sample_entry['DATA'] = sample_entry['DATA'].format(sample=sample)

    return sample_entry

def variant_global_sample_info_entry(sample):
    """
    Get an entry from SAMPLE_INFO_TABLE.

    :param sample: Sample name.

    :return: A configuration entry (Pandas.Series) or None if there is no sample information.
    """

    if sample not in SAMPLE_INFO_TABLE.index:
        return None

    return SAMPLE_INFO_TABLE.loc[sample]

def variant_global_sample_info_entry_element(sample, element, default=None, null_default=True):
    """
    Get an entry from SAMPLE_INFO_TABLE.

    :param sample: Sample name (sample info table row).
    :param element: Element name (sample info table column).
    :param default: Default value if sample or element is not defined.
    :param null_default: If the element is present but is null (`np.nan`), then return `default` instead of `np.nan`.

    :return: A configuration entry (Pandas.Series) or None if there is no sample information.
    """

    # Get entry
    sample_info = variant_global_sample_info_entry(sample)

    if sample_info is None:
        return default

    # Get element value
    if element not in sample_info:
        return default

    value = sample_info[element]

    # Translate null to default (if null_default)
    if null_default and pd.isnull(value):
        return default

    # Return value
    return value


#############
### Rules ###
#############

# variant_global_filter_fa
rule variant_global_filter_fa:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        fa='temp/variant/caller/{sourcename}/{sample}/all/all/bed/pre_filter/fa/{vartype}_{svtype}.fa.gz'
    output:
        fa='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz'
    run:

        # Read variant IDs
        id_set = set(
            pd.read_csv(input.bed, sep='\t', usecols=('ID', ), squeeze=True)
        )

        # Filter
        with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
            SeqIO.write(svpoplib.seq.fa_to_record_iter(input.fa, id_set), out_file, 'fasta')


# variant_global_filter_region
#
# Apply a BED filter to SVs.
rule variant_global_filter_region:
    input:
        bed='temp/variant/caller/{sourcename}/{sample}/all/all/bed/pre_filter/{vartype}_{svtype}.bed.gz',
        filter=lambda wildcards: svpoplib.variant.get_filter_bed(wildcards.filter, UCSC_REF_NAME, config, SVPOP_DIR)
    output:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        bed_filt='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/filter_dropped/{vartype}_{svtype}_dropped.bed.gz'
    wildcard_constraints:
        filter='((?!all)[^\/]*)|all.+',
        svtype='ins|del|inv|dup|rgn|sub'
    run:

        if wildcards.filter != 'all':
            shell(
                """bedtools intersect -wa -v -sorted -a {input.bed} -b {input.filter} -header | gzip > {output.bed}; """
                """bedtools intersect -wa -u -sorted -a {input.bed} -b {input.filter} -header | gzip > {output.bed_filt}; """
            )
        else:
            shell(
                """cp {input.bed} {output.bed}; """
                """touch {output.bed_filt}; """
            )

# variant_global_byref_to_bylen
#
# Variant byref to bylen.
rule variant_global_byref_to_bylen:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/bed/{sample}/{svset}/{filter}/byref/{vartype}_{svtype}.bed'
    output:
        bed='results/variant/{sourcetype}/{sourcename}/bed/{sample}/{svset}/{filter}/bylen/{vartype}_{svtype}.bed'
    run:

        if wildcards.svtype == 'ins':
            # Read
            df = pd.read_csv(input.bed, sep='\t', header=0)

            df['END'] = df['POS'] + df['SVLEN']
            df.to_csv(output.bed, sep='\t', index=False)

        else:
            shell(
                """ln -sf ../byref/{wildcards.vartype}_{wildcards.svtype}.bed {output.bed}; """
            )

# variant_global_vcf_gz
#
# Compress VCF files
rule variant_global_vcf_gz:
    input:
        vcf='temp/{sourcetype}/{sourcename}/vcf/{sample}/{varset}/{filter}/{ref}_{fmt}/{vartype}_{svtype}.vcf'
    output:
        vcf='results/{sourcetype}/{sourcename}/vcf/{sample}/{varset}/{filter}/{ref,grc|hg}_{fmt,sv|alt}/{vartype}_{svtype}.vcf.gz',
        tbi='results/{sourcetype}/{sourcename}/vcf/{sample}/{varset}/{filter}/{ref,grc|hg}_{fmt,sv|alt}/{vartype}_{svtype}.vcf.gz.tbi'
    shell:
        """bgzip -c {input.vcf} > {output.vcf}; """
        """sleep 5; """
        """tabix {output.vcf}"""

# variant_global_uncompress_fa
#
# Uncompress for tools that cannot read a gzipped FASTA.
rule variant_global_uncompress_fa:
    input:
        fa='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz'
    output:
        fa=temp('temp/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa')
    run:

        if os.stat(input.fa).st_size > 0:
            shell("""zcat {input.fa} > {output.fa}""")
        else:
            shell("""> {output.fa}""")

# variant_global_variant_fai
#
# Create FAI for variant FASTA files.
rule variant_global_variant_fai:
    input:
        fa='results/variant/{sourcetype}/{sourcename}/fasta/{sample}/{varset}/{filter}/{vartype}_{svtype}.fa.gz'
    output:
        fai='results/variant/{sourcetype}/{sourcename}/fasta/{sample}/{varset}/{filter}/{vartype}_{svtype}.fa.gz.fai'
    shell:
         """samtools faidx {input.fa}"""
