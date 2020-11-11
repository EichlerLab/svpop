"""
Make duplication calls based on sample read depth.
"""


###################
### Definitions ###
###################

def _variant_smrtsv_dup_get_config(map_profile, config):

    # Check config section
    if 'smrtsv_map' not in config:
        raise RuntimeError('Missing config section: "smrtsv_map"')

    if 'dup' not in config['smrtsv_map']:
        raise RuntimeError('Missing config section: "dup" in "smrtsv_map"')

    if map_profile not in config['smrtsv_map']['dup']:
        raise RuntimeError('Missing config section: "{}" in "smrtsv_map/dup"'.format(map_profile))

    dup_entry = config['smrtsv_map']['dup'][map_profile]

    # Check entry
    missing_elements = [element for element in ['win_size', 'mapq', 'distance', 'flank', 'autosome_z', 'chrx_m_z', 'fallout_max', 'min_len'] if element not in dup_entry]

    if missing_elements:
        raise RuntimeError('Missing element(s) from config entry (params/smrtsv/dup): {}'.format(', '.join(missing_elements)))

    # Set types
    dup_entry['win_size'] = np.int32(dup_entry['win_size'])
    dup_entry['mapq'] = np.int32(dup_entry['mapq'])
    dup_entry['distance'] = np.int32(dup_entry['distance'])
    dup_entry['flank'] = np.int32(dup_entry['flank'])
    dup_entry['autosome_z'] = np.float32(dup_entry['autosome_z'])
    dup_entry['chrx_m_z'] = np.float32(dup_entry['chrx_m_z'])
    dup_entry['fallout_max'] = np.int32(dup_entry['fallout_max'])
    dup_entry['min_len'] = np.int32(dup_entry['min_len'])

    # Return
    return dup_entry


def _bedcnv_merge_consecutive_regions(df, fallout_max):
    """
    Merge copy-number regions by chromosome.

    :param df: Dataframe of one chromosome with column 'SIG' to identify regions of significance that should be
        merged. Frame index must be a list of consecutive integers.
    :param sig_list: List of True/False values the size of df where True indicates that
    """

    # Set type list
    df_merged_dtype = {
        'POS': np.int32,
        'END': np.int32,
        'Z_MEAN': np.float32,
        'DEPTH_MEAN': np.float32,
        'DEPTH_MIN': np.float32,
        'DEPTH_MAX': np.float32
    }

    # Reset index
    df.reset_index(drop=True, inplace=True)

    # Get a list of indicies for significant loci
    sig_index_list = list(df.index[df['SIG']])

    if len(sig_index_list) == 0:
        return pd.DataFrame(columns=['POS', 'END', 'DEPTH_MEAN', 'DEPTH_MIN', 'DEPTH_MAX', 'Z_MEAN']).astype(df_merged_dtype)

    # Initialize for range traversal (merges consecutive indices into a continuous range)
    range_list = list()  # List of ranges merged from indices

    sig_start = sig_index_list[0]  # Start of the current range
    sig_last = sig_start           # Last index in the current range

    # Traverse
    for sig_index in sig_index_list[1:]:

        # Next range if index is too far away from last index
        if sig_index - sig_last > fallout_max:
            range_list.append((sig_start, sig_last))
            sig_start = sig_index

        sig_last = sig_index

    range_list.append((sig_start, sig_last))

    # Get records for each range
    df_merged_list = list()

    for sig_start, sig_end in range_list:
        z_list = df.loc[list(range(sig_start, sig_end + 1)), 'Z']
        depth_list = df.loc[list(range(sig_start, sig_end + 1)), 'DEPTH']

        df_merged_list.append(pd.Series(
            [
                df.loc[sig_start, 'POS'],
                df.loc[sig_end, 'END'],
                np.mean(depth_list),
                np.min(depth_list),
                np.max(depth_list),
                np.mean(z_list)
            ],
            index=['POS', 'END', 'DEPTH_MEAN', 'DEPTH_MIN', 'DEPTH_MAX', 'Z_MEAN']
        ))

    # Merge
    return pd.concat(df_merged_list, axis=1).T.astype(df_merged_dtype)


def _variant_caller_smrtsv_dup_call_input_bed(wildcards):
    """
    Get input file for DUP caller. This is a BED file of windows with read depth and z-scores calculated for each. The
    BED file is already filtered by a filter region as well as DUP and GAP UCSC tracks.
    """

    # Check
    if 'map_profile' not in wildcards.keys():
        raise RuntimeError('Cannot get input file for DUP calling: Missing "map_profile" in wildcards (_variant_caller_smrtsv_dup_call_input_bed())')

    # Input file pattern
    file_pattern = 'results/variant/caller/smrtsv/dup/counts/{sample}/{filter}_{distance}_{flank}/counts_{win_size}_{mapq}.bed.gz'

    # Get config entry (has parameters for DUP calls, such as Z-score and window size)
    dup_entry = _variant_smrtsv_dup_get_config(wildcards.map_profile, config)

    # Add wildcards
    dup_entry['sample'] = wildcards.sample
    dup_entry['filter'] = wildcards.filter

    # Parse
    return file_pattern.format(**dup_entry)


#############
### Rules ###
#############


# variant_caller_smrtsv_dup_copy_default
#
# Copy final call to output.
#
# Note: Rule "variant_smrtsv_bed_fa" writes the final BED and FA files after this rule.
rule variant_caller_smrtsv_dup_copy_default:
    input:
        bed='results/variant/caller/smrtsv/dup/call/{sample}/{filter}/profile_default/sv_dup.bed'
    output:
        bed=temp('temp/variant/caller/smrtsv/bed/{sample}/all/{filter}/byref/sv_dup.bed')
    shell:
        """cp -l {input.bed} {output.bed}"""

# variant_caller_smrtsv_dup_call
#
# Call regions of elevated copy number
rule variant_caller_smrtsv_dup_call:
    input:
        ref=config['reference'],
        bed=_variant_caller_smrtsv_dup_call_input_bed
    output:
        bed='results/variant/caller/smrtsv/dup/call/{sample}/{filter}/profile_{map_profile}/sv_dup.bed'
    run:

        # Get config entry (has parameters for DUP calls, such as Z-score and window size)
        dup_entry = _variant_smrtsv_dup_get_config(wildcards.map_profile, config)

        autosome_z = dup_entry['autosome_z']
        chrx_z = dup_entry['chrx_m_z']
        chry_z = dup_entry['chrx_m_z']

        # Adjust chrX z-score based on sex
        if variant_global_sample_info_entry_element(wildcards.sample, 'SEX', 'U').upper() != 'M':
            chrx_z = autosome_z

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Warn if chrX was not found. Could be against the wrong reference (e.g. "X" instead of "chrX")
        if not np.any(df['#CHROM'] == 'chrX') and df.shape[0] > 0:
            raise RuntimeError('Expecting GRCh38 chromosome names in BED file, but no "chrX" was found: {}'.format(input.bed))

        # Flag regions of significance
        df['SIG'] = df['Z'] >= df['#CHROM'].apply(
            lambda chrom: (chrx_z if chrom == 'chrX' else (chry_z if chrom == 'chrY' else autosome_z))
        )

        # Merge consecutive significant regions into DUP calls
        df = df.groupby('#CHROM').apply(_bedcnv_merge_consecutive_regions, fallout_max=dup_entry['fallout_max'])

        # Reset index
        df.index = df.index.droplevel(1)  # Drop the original index (numeric) and keep #CHROM
        df.reset_index(inplace=True)

        # Set SV columns
        df['SVLEN'] = df['END'] - df['POS']
        df['SVTYPE'] = 'DUP'
        df['ID'] = analib.variant.get_variant_id(df)

        # Filter by length
        df = df.loc[df['SVLEN'] >= dup_entry['min_len']]

        # Set type and ID
        df['SVTYPE'] = 'DUP'
        df['ID'] = analib.variant.get_variant_id(df)

        # Set sequence
        with pysam.FastaFile(input.ref) as fasta_file:
            df['SEQ'] = df.apply(lambda row: fasta_file.fetch(row['#CHROM'], row['POS'], row['END']), axis=1)

        # Sort columns
        df = analib.variant.order_variant_columns(df)

        # Write
        df.to_csv(output.bed, sep='\t', float_format='%0.4f', index=False)


# variant_caller_smrtsv_dup_score_depth
#
# Z-Score depth per window.
rule variant_caller_smrtsv_dup_score_depth:
    input:
        bed='temp/variant/caller/smrtsv/dup/data/{sample}/{filter}_{distance}_{flank}/depth_{win_size}_{mapq}.bed.gz'
    output:
        bed='results/variant/caller/smrtsv/dup/counts/{sample}/{filter}_{distance}_{flank}/counts_{win_size}_{mapq}.bed.gz'
    run:

        # Define data types
        DATA_TYPES = {
            '#CHROM': np.object,
            'POS': np.int32,
            'END': np.int32,
            'DEPTH': np.float32
        }

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0, dtype=DATA_TYPES)

        # Get Z
        mean = np.mean(df['DEPTH'])
        std = np.std(df['DEPTH'])

        df['Z'] = ((df['DEPTH'] - mean) / std).astype(np.float32)

        # Write all
        df.to_csv(output.bed, sep='\t', index=False, float_format='%0.4f', compression='gzip')

# variant_caller_smrtsv_dup_apply_filter
#
# Filter alignment depth BED.
rule variant_caller_smrtsv_dup_apply_filter:
    input:
        bed='results/depth/sample/{sample}/depth/depth_{win_size}_{mapq}.bed.gz',
        bed_filt='results/variant/caller/smrtsv/dup/data/{filter}_{distance}_{flank}/filter.bed'
    output:
        bed=temp('temp/variant/caller/smrtsv/dup/data/{sample}/{filter}_{distance}_{flank}/depth_{win_size}_{mapq}.bed.gz')
    shell:
        """bedtools intersect -a {input.bed} -b {input.bed_filt} -sorted -header -v | """
        """gzip """
        """> {output.bed}"""

# variant_caller_smrtsv_dup_make_filter
#
# Make DUP filter including centromeres, gaps, and explicitly filtered loci.
rule variant_caller_smrtsv_dup_make_filter:
    input:
        fai=config['reference'] + '.fai',
        bed_filt=os.path.join(SVPOP_DIR, 'files/filter/{}/{{filter}}.bed'.format(UCSC_REF_NAME)),
        bed_gap='data/anno/gap/gap_regions_{distance}_0.bed',
        bed_cen='data/anno/cen/cen_regions_{distance}_0.bed'
    output:
        bed='results/variant/caller/smrtsv/dup/data/{filter}_{distance}_{flank}/filter.bed'
    shell:
        """echo -e "#CHROM\tPOS\tEND" > {output.bed}; """
        """awk -vOFS="\t" '($1 !~ /^#/) {{print $1, $2, $3}}' {input.bed_filt} {input.bed_gap} {input.bed_cen} | """
        """sort -k1,1 -k2,2n | """
        """bedtools merge -d {wildcards.distance} | """
        """bedtools slop -g {input.fai} -b {wildcards.flank} """
        """>> {output.bed}"""
