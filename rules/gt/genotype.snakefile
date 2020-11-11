"""
Rules for processing and summarizing genotypes.
"""


###################
### Definitions ###
###################

# Convert a genotype call to a keyword
GENOTYPE_TO_AC = {
    './.': -1,
    '0/0': 0,
    '1/0': 1,
    '0/1': 1,
    '1/1': 2,
    '.|.': -1,
    '0|0': 0,
    '1|0': 1,
    '0|1': 1,
    '1|1': 2
}

def _gt_get_input_vcf(wildcards):
    """
    Get input genotype VCF.
    """

    # Check
    if wildcards.gtset not in config['genotype']['sets']:
        raise RuntimeError('Genotype set "{gtset}" is not configured (config["genotype"]["sets"]["{gtset}"] (_gt_get_input_vcf)'.format(**wildcards))

    if 'gt_vcf' not in config['genotype']['sets'][wildcards.gtset]:
        raise RuntimeError('Section "gt_vcf" missing from the configuration for set "{gtset}" (config["genotype"]["sets"]["{gtset}"]["gt_vcf"] (_gt_get_input_vcf)'.format(**wildcards))

    # Return VCF
    return config['genotype']['sets'][wildcards.gtset]['gt_vcf']

def _gt_get_sv_bed(wildcards):
    """
    Get SV BED file for a genotype set.

    'wildcards' must contain all or none of the elements "svset", "filter", and "svtype". If these fields are in
    `wildcards`, then they are used to locate the SV BED file. If they are not, then the genotype set configuration
    is used to locate the BED file. These three elements are used to locate BED files for filtering genotype calls, and
    they are left out (defined by the configuration) for reading the SV calls to filter the genotype calls. In both
    cases, the sample name is found in the genotype configuration using the "sample" element.
    """

    # Setup param dict
    param_dict = dict()

    # Get GT Set config
    if wildcards.gtset not in config['genotype']['sets']:
        raise RuntimeError('Genotype set "{gtset}" is not configured (config["genotype"]["sets"]["{gtset}"] (gt_get_sv_bed)'.format(**wildcards))

    config_gtset = config['genotype']['sets'][wildcards.gtset]

    if 'sv_bed' not in config_gtset:
        raise RuntimeError('Param "sv_bed" is not in the configuration for gtset {}'.format(wildcards.gtset))

    return config_gtset['sv_bed']

def _gt_get_sample_info_table(wildcards):
    """
    Get a table of sample information (population, superpopulation, sex).
    """

    # Check
    if wildcards.gtset not in config['genotype']['sets']:
        raise RuntimeError('Genotype set "{gtset}" is not configured (config["genotype"]["sets"]["{gtset}"] (_gt_get_sample_info_table)'.format(**wildcards))

    config_gtset = config['genotype']['sets'][wildcards.gtset]

    if 'sample_table' not in config_gtset:
        raise RuntimeError('Missing "sample_table" in genotype config for "{gtset}" (config["genotype"]["sets"]["{gtset}"]["sample_table"] (_gt_get_sample_info_table)'.format(**wildcards))

    return config_gtset['sample_table']




#############
### Rules ###
#############


# gt_info_pop_vs_pop
#
# Get F statistics for a population group vs others.
rule gt_info_pop_vs_pop:
    input:
        call='results/genotype/{gtset}/tables/{svset}/{filter}/call_{callpct}/{svtype}/call_table.tab.gz',
        sample_info=_gt_get_sample_info_table
    output:
        tab='results/genotype/{gtset}/tables/{svset}/{filter}/call_{callpct}/{svtype}/pop_diff/{region_a}_vs_{region_b}/pop_diff_table.tab',
        summary='results/genotype/{gtset}/tables/{svset}/{filter}/call_{callpct}/{svtype}/pop_diff/{region_a}_vs_{region_b}/pop_diff_summary.tab'
    run:

        raise RuntimeError('Genotyping changed from symbolic (NO_CALL, HET, etc) to numeric (-1, 1, etc). This rule (gt_info_pop_vs_pop) needs to be updated and tested ')

        # Get region A
        tok = wildcards.region_a.strip().split('_', 1)

        if len(tok) == 1 or tok[0].upper() not in {'POP', 'SUPERPOP'}:
            raise RuntimeError('Population A must begin with "pop_" or "superpop_": {}'.format(wildcards.region_a))

        region_a_level = tok[0].upper()
        region_a = tok[1].upper()

        # Get region B
        tok = wildcards.region_b.strip().split('_', 1)

        if len(tok) == 1:
            if tok[0].upper() != 'ALL':
                raise RuntimeError('Population B must be "all" or begin with "pop_" or "superpop_": Found single word that is not "all": {}'.format(wildcards.region_b))

            region_b_level = 'ALL'
            region_b = 'ALL'

        elif tok[0].upper() in {'POP', 'SUPERPOP'}:
            region_b_level = tok[0].upper()
            region_b = tok[1].upper()

        else:
            raise RuntimeError('Population B must be "all" or begin with "pop_" or "superpop_": {}'.format(wildcards.region_b))


        # Read
        sample_info = pd.read_csv(
            input.sample_info, sep='\t', header=0, comment='#', index_col='SAMPLE', usecols=('SAMPLE', 'POP', 'SUPERPOP')
        )

        sample_info = sample_info.fillna('UNKNOWN')

        if len(set(sample_info.index)) != sample_info.shape[0]:
            raise RuntimeError('Sample info contains duplicate sample records: {}'.format(input.sample_info))

        df_call = pd.read_csv(input.call, sep='\t', index_col='ID', header=0)

        if len(set(df_call.columns)) != df_call.shape[1]:
            raise RuntimeError('Genotype call table contains duplicate sample records: {}'.format(input.call))

        n_sv = df_call.shape[0]

        # Report missing samples
        sample_table_set = set(sample_info.index)

        missing_samples = [gt_sample for gt_sample in df_call.columns if gt_sample not in sample_table_set]

        if len(missing_samples) > 0:
            missing_sample_out = ', '.join(missing_samples[0:3])

            if len('missing_samples') > 3:
                missing_sample_out += '...'

            raise RuntimeError('Genotype table contains samples missing in the sample table: {}'.format(missing_sample_out))

        # Subset sample info table
        sample_info = sample_info.loc[df_call.columns]

        # Group a must be a subset of b, or a and b must be mutually exclusive
        pop_a = list(sample_info.loc[sample_info.loc[: , region_a_level] == region_a].index)

        if region_b_level == 'ALL':
            pop_b = list(sample_info.index)
        else:
            pop_b = list(sample_info.loc[sample_info.loc[: , region_b_level] == region_b].index)

        pop_a_set = set(pop_a)
        pop_b_set = set(pop_b)

        if pop_a_set == pop_b_set:
            # Cannot be equal
            raise RuntimeError('Population regions yield the same genotype samples: a = {region_a}, b = {region_b}'.format(**wildcards))

        if pop_a_set < pop_b_set:  # Subset, but not same set
            # A is a subset of B and A is not B
            pop_b_set = pop_b_set.difference(pop_a_set)

            pop_b = [sample for sample in pop_b if sample in pop_b_set]

        elif len(pop_a_set & pop_b_set) != 0:
            raise RuntimeError('Population regions yield overlapping genotype samples and A is not strictly a subset of B: a = {region_a}, b = {region_b}'.format(**wildcards))

        del(pop_a_set)
        del(pop_b_set)

        # Check group a and b (neither may be empty)
        if len(pop_a) == 0:
            raise RuntimeError('Found 0 genotype samples in group A ({})'.format(wildcards.region_a))

        if len(pop_b) == 0:
            raise RuntimeError('Found 0 genotype samples in group B ({})'.format(wildcards.region_b))


        # Subset genotype table
        df_a = df_call.loc[:, pop_a]
        df_b = df_call.loc[:, pop_b]

        del(df_call)


        # Get number of samples and filter no-call
        n_a = df_a.shape[1]
        n_b = df_b.shape[1]

        # Get summary stats
        df_stat = pd.concat(
            [
                n_a - df_a.apply(lambda col: np.sum(col == 'NO_CALL'), axis=1),
                df_a.apply(lambda col: np.sum(col == 'NO_CALL'), axis=1),
                df_a.apply(lambda col: np.sum(col == 'HOM_REF'), axis=1),
                df_a.apply(lambda col: np.sum(col == 'HET'), axis=1),
                df_a.apply(lambda col: np.sum(col == 'HOM_ALT'), axis=1),
                n_b - df_b.apply(lambda col: np.sum(col == 'NO_CALL'), axis=1),
                df_b.apply(lambda col: np.sum(col == 'NO_CALL'), axis=1),
                df_b.apply(lambda col: np.sum(col == 'HOM_REF'), axis=1),
                df_b.apply(lambda col: np.sum(col == 'HET'), axis=1),
                df_b.apply(lambda col: np.sum(col == 'HOM_ALT'), axis=1)
            ],
            axis=1
        )

        df_stat.columns = pd.MultiIndex.from_tuples(
            list(itertools.product(('POP_A', 'POP_B'), ('N', 'NO_CALL', 'HOM_REF', 'HET', 'HOM_ALT'))),
            names=('POP', 'STAT')
        )

        # Set allele frequency (AF)
        df_stat.insert(
            5,
            ('POP_A', 'AF'),
            df_stat.loc[:, 'POP_A'].apply(
                lambda row: (row['HET'] + 2 * row['HOM_ALT']) / (2 * (row['N'])) if row['N'] > 0 else np.nan
                , axis=1
            )
        )

        df_stat.insert(
            11,
            ('POP_B', 'AF'),
            df_stat.loc[:, 'POP_B'].apply(
                lambda row: (row['HET'] + 2 * row['HOM_ALT']) / (2 * (row['N'])) if row['N'] > 0 else np.nan
                , axis=1
            )
        )

        # Calculate FST
        df_stat[('STAT', 'FST')] = df_stat.apply(analib.gt.fst_wc, axis=1)
        df_stat[('STAT', 'AF_DIFF')] = df_stat[('POP_A', 'AF')] - df_stat[('POP_B', 'AF')]

        # Make summary table
        fst_list = df_stat.loc[:, ('STAT', 'FST')]

        fst_list = fst_list.loc[~ pd.isnull(fst_list)]

        df_summary = pd.Series(
            [
                n_a,
                n_b,
                n_sv,
                fst_list.shape[0],
                np.mean(fst_list),
                np.std(fst_list),
                np.max(fst_list),
                np.sum(fst_list >= 0.2),
                np.sum(fst_list >= 0.5),
                np.sum(fst_list >= 0.8)
            ],
            index=['N_A', 'N_B', 'N_SV', 'FST_N', 'FST_MEAN', 'FST_STD', 'FST_MAX', 'N_FST_20', 'N_FST_50', 'N_FST_80']
        )

        df_summary.name = 'VAL'
        df_summary.index.name = 'STAT'

        # Write
        df_stat.index.name = ('SV', 'ID')
        df_stat.reset_index().to_csv(output.tab, sep='\t', index=False)

        df_summary.to_csv(output.summary, sep='\t', index=True, header=True)


# gt_info_filter
#
# Filter each variant by the proportion of genotype samples a call could be made in and any other filter that could be
# applied to the SV bed (e.g. svset and svtype).
rule gt_info_filter:
    input:
        bed=_gt_get_sv_bed,
        info='results/genotype/{gtset}/gt_data/info_table.bed.gz',
        call='results/genotype/{gtset}/gt_data/call_table.tab.gz'
    output:
        info='results/genotype/{gtset}/tables/{svset}/{filter}/call_{callpct}/{svtype}/info_table.bed.gz',
        call='results/genotype/{gtset}/tables/{svset}/{filter}/call_{callpct}/{svtype}/call_table.tab.gz',
        info_dropped='results/genotype/{gtset}/tables/{svset}/{filter}/call_{callpct}/{svtype}/info_table_dropped.bed.gz'
    run:

        # Get call proportion
        call_prop = np.float64(np.int16(wildcards.callpct) / 100)

        if not (0.0 <= call_prop <= 1.0):
            raise RuntimeError('Call proportion must be between 0 and 1 (inclusive): {}'.format(call_prop))

        # Read SV IDs
        sv_id = list(pd.read_csv(input.bed, sep='\t', header=0, usecols=('ID', ), squeeze=True))

        # Read info
        df_info = pd.read_csv(input.info, sep='\t', header=0)
        df_info.set_index('ID', inplace=True, drop=False)

        if not all([id in df_info['ID'] for id in sv_id]):
            raise RuntimeError('SV BED contains variant IDs not in genotype info table: {}'.format(input.info))

        df_info = df_info.loc[sv_id]


        # Read calls
        df_call = pd.read_csv(input.call, sep='\t', header=0)
        df_call.set_index('ID', inplace=True, drop=False)

        if not all([id in df_call['ID'] for id in sv_id]):
            raise RuntimeError('SV BED contains variant IDs not in genotype call table: {}'.format(input.call))

        df_call = df_call.loc[sv_id]


        # Subset SV IDs by genotype proportion, then by filtered SVs
        callable_id = list(df_info.loc[df_info['PROP_CALL'] >= call_prop, 'ID'])

        callable_id = [id for id in callable_id if id in sv_id]


        # Drop genotype records
        df_info.loc[
            [id not in callable_id for id in df_info.index]
        ].to_csv(
            output.info_dropped, sep='\t', index=False, compression='gzip'
        )

        df_info = df_info.loc[callable_id]
        df_call = df_call.loc[callable_id]


        # Write
        df_info.to_csv(output.info, sep='\t', index=False, compression='gzip')
        df_call.to_csv(output.call, sep='\t', index=False, compression='gzip')


# gt_info_table
#
# Summarize genotype information for each variant.
rule gt_info_table:
    input:
        bed=_gt_get_sv_bed,
        call='results/genotype/{gtset}/gt_data/call_table.tab.gz'
    output:
        info='results/genotype/{gtset}/gt_data/info_table.bed.gz',
        nc_samples='results/genotype/{gtset}/gt_data/no_call_samples.txt'
    run:

        # Read variant information
        df_sv = pd.read_csv(input.bed, sep='\t', header=0)

        # Subset columns
        head_cols = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN']
        tail_cols = list()

        for col in (
            'CALLERSET_N', 'CALLERSET_SOURCE', 'CALLERSET_LIST',
            'MERGE_SRC', 'MERGE_AC', 'MERGE_AF', 'MERGE_SAMPLES', 'DISC_CLASS'
        ):
            if col in df_sv.columns:
                tail_cols.append(col)

        df_sv = df_sv.loc[:, head_cols + tail_cols]

        # Set index
        df_sv.set_index('ID', inplace=True, drop=False)

        # No duplicate SV IDs in SV BED
        dup_list = [id for id, count in collections.Counter(df_sv['ID']).items() if count > 1]

        if len(dup_list) > 0:
            raise RuntimeError('Found {} SV IDs in variant calls with duplicate entries: "{}"...'.format(len(dup_list), dup_list[0]))

        # Read genotype info
        df_gt = pd.read_csv(input.call, sep='\t', header=0, index_col='ID')

        # No duplicated SV IDs in genotypes
        dup_list = [id for id, count in collections.Counter(df_gt.index).items() if count > 1]

        if len(dup_list) > 0:
            raise RuntimeError('Found {} SV IDs in genotype calls with duplicate entries: "{}"...'.format(len(dup_list), dup_list[0]))

        # All genotype SV IDs must be in the SV table
        missing_list = sorted(set(df_gt.index) - set(df_sv['ID']))

        if missing_list:
            raise RuntimeError('Found {} SV IDs in genotype table that are not present in the SV IDs: {}'.format(
                len(missing_list),
                ', '.join(missing_list[:3]) + (', ...' if len(missing_list) > 3 else '')
            ))

        # Arrange genotypes
        df_sv = df_sv.loc[df_sv['ID'].apply(lambda val: val in df_gt.index)]
        df_gt = df_gt.loc[df_sv['ID']]

        # Write list of no-call samples
        no_call_list = df_gt.columns[df_gt.apply(lambda col: np.all(col == -1), axis=0)].tolist()

        with open(output.nc_samples, 'w') as out_file:
            out_file.write('\n'.join(no_call_list))

            if len(no_call_list) > 0:
                out_file.write('\n')

        # Get counts
        n_nc = df_gt.apply(lambda row: sum(row == -1), axis=1)
        n_nc.name = 'N_NO_CALL'

        n_snc = df_gt.apply(lambda row: sum(row != -1), axis=1)
        n_snc.name = 'N_HAS_CALL'

        n_hom_ref = df_gt.apply(lambda row: sum(row == 0), axis=1)
        n_hom_ref.name = 'N_HOM_REF'

        n_het = df_gt.apply(lambda row: sum(row == 1), axis=1)
        n_het.name = 'N_HET'

        n_hom_alt = df_gt.apply(lambda row: sum(row == 2), axis=1)
        n_hom_alt.name = 'N_HOM_ALT'

        n_hom = n_hom_ref + n_hom_alt
        n_hom.name = 'N_HOM'

        # Get genotyping stats
        stat_ac = (n_hom_alt * 2 + n_het)
        stat_ac.name = 'AC'

        stat_af = stat_ac / (n_snc * 2)
        stat_af.name = 'AF'

        stat_ns = n_nc + n_snc
        stat_ns.name = 'NS'

        stat_an = n_snc * 2
        stat_an.name = 'AN'

        # Get proportions
        prop_nc = n_nc / (n_nc + n_snc)
        prop_nc.name = 'PROP_NO_CALL'

        prop_snc = n_snc / (n_nc + n_snc)
        prop_snc.name = 'PROP_CALL'

        prop_het = n_het / n_snc
        prop_het.name = 'PROP_HET'

        exp_prop_het = (2 * stat_af * (1 - stat_af))
        exp_prop_het.name = 'EXP_PROP_HET'

        # Write
        pd.concat(
            [df_sv, n_nc, n_snc, n_het, n_hom_alt, n_hom_ref, n_hom, stat_ac, stat_ns, stat_an, stat_af, prop_nc, prop_snc, prop_het, exp_prop_het],
            axis=1
        ).to_csv(output.info, sep='\t', index=False, na_rep='NA', compression='gzip')


# gt_call_table
#
# Make call table.
rule gt_call_table:
    input:
        vcf='results/genotype/{gtset}/gt_data/gt_vcf_records.tab.gz'
    output:
        calls='results/genotype/{gtset}/gt_data/call_table.tab.gz'
    run:

        # Read
        df = pd.read_csv(input.vcf, sep='\t', header=0, index_col='ID')
        df = df.loc[:, [col for col in df.columns if col not in ('#CHROM', 'REF', 'POS', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')]]

        # Get genotype string
        for col_name in df.columns:
            df[col_name] = df[col_name].apply(lambda val: GENOTYPE_TO_AC[val.split(':')[0]]).astype(np.int32)

        # Write
        df.to_csv(output.calls, sep='\t', compression='gzip')

# gt_uncompress_vcf
#
# Uncompress VCF.
rule gt_uncompress_vcf:
    input:
        vcf=_gt_get_input_vcf
    output:
        tab='results/genotype/{gtset}/gt_data/gt_vcf_records.tab.gz',
        headers='results/genotype/{gtset}/gt_data/gt_vcf_headers.txt'
    run:

        with gzip.open(input.vcf, 'rt') as in_file:

            line = next(in_file)

            # Write headers
            with open(output.headers, 'w') as out_file:
                while line and line.startswith('##'):
                    out_file.write(line)

                    try:
                        line = next(in_file)
                    except StopIteration:
                        line = ''
                        break

            # Write VCF records
            with gzip.open(output.tab, 'w') as out_file:
                while line:
                    out_file.write(line.encode())

                    try:
                        line = next(in_file)
                    except StopIteration:
                        line = ''
                        break
