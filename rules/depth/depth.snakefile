"""
Get alignment depth in windows over the genome.
"""

#
# Definitions
#

def _depth_data_sample_get_depth_merge_input(wildcards):

    file_pattern = 'temp/depth/sample/{sample}/depth/depth_{win_size}_{mapq}/{index}.bed.gz'

    # Get a list of samples
    index_list = list(pd.read_csv(
        'results/depth/sample/{sample}/tables/bam_table.tab'.format(**wildcards),
        sep='\t',
        usecols=('INDEX', ),
        squeeze=True
    ))

    # Parse into file pattern
    format_dict = {
        'sample': wildcards.sample,
        'win_size': wildcards.win_size,
        'mapq': wildcards.mapq
    }

    file_list = list()

    for index in index_list:
        format_dict['index'] = str(index)
        file_list.append(file_pattern.format(**format_dict))

    # Return file list
    return file_list


#
# Depth wiggle
#

# depth_data_sample_bigwig
#
# Wiggle to BigWig. The BigWig file output can be loaded into UCSC to show a read depth track.
rule depth_data_sample_bigwig:
    input:
        wig='temp/depth/sample/{sample}/depth/wig/depth_{win_size}_{mapq}.wig',
        fai=config['reference'] + '.fai'
    output:
        bw='results/depth/sample/{sample}/depth/wig/depth_{win_size}_{mapq}.bw'
    shell:
        """wigToBigWig {input.wig} {input.fai} {output.bw}"""

# depth_data_sample_wiggle
#
# Make wiggle file from alignment depths. This file can be turned into a UCSC track by converting it to BigWig format.
rule depth_data_sample_wiggle:
    input:
        bed='results/depth/sample/{sample}/depth/depth_{win_size}_{mapq}.bed.gz'
    output:
        wig=temp('temp/depth/sample/{sample}/depth/wig/depth_{win_size}_{mapq}.wig')
    run:

        win_size = int(wildcards.win_size)

        # Read
        df = pd.read_csv(input.bed, sep='\t')

        # Write wiggle
        last_chrom = ''

        with open(output.wig, 'w') as out_file:

            for index, row in df.iterrows():

                # Write CHROM header
                if row['#CHROM'] != last_chrom:
                    out_file.write(
                        'fixedStep chrom={0} start={1} step={2}\n'.format(
                            row['#CHROM'], row['POS'] + 1, wildcards.win_size
                        )
                    )

                    last_chrom = row['#CHROM']

                # Write record
                if (row['END'] - row['POS']) == win_size:
                    out_file.write('{0:.2f}\n'.format(row['DEPTH']))


#
# Depth
#

# depth_data_sample_depth_merge
#
# Merge depth.
rule depth_data_sample_depth_merge:
    input:
        tab='results/depth/sample/{sample}/tables/bam_table.tab',
        bed=_depth_data_sample_get_depth_merge_input
    output:
        bed='results/depth/sample/{sample}/depth/depth_{win_size}_{mapq}.bed.gz'
    run:

        # Definitions
        column_types = {
            '#CHROM': np.object,
            'POS': np.int32,
            'END': np.int32,
            'DEPTH': np.float32
        }

        # Init
        df = None

        # Read each BED - merge counts
        for bed_file in input.bed:
            df_next = pd.read_csv(
                bed_file, sep='\t', header=0,
                index_col=('#CHROM', 'POS', 'END'),
                usecols=('#CHROM', 'POS', 'END', 'DEPTH'),
                dtype=column_types,
                squeeze=True
            )

            if df is not None:
                df += df_next

            else:
                df = df_next

        # As DataFrame
        df = pd.DataFrame(df).reset_index()

        df['POS'] = df['POS'].astype(np.int32)
        df['END'] = df['END'].astype(np.int32)
        df['DEPTH'] = df['DEPTH'].astype(np.float32)

        # Sort
        df.sort_values(['#CHROM', 'POS', 'END'], inplace=True)

        # Write
        df.to_csv(output.bed, sep='\t', index=False, header=True, float_format='%0.4f', compression='gzip')


# depth_data_sample_depth_per_index
#
# Estimate alignment depth for each alignment file in each window.
rule depth_data_sample_depth_per_index:
    input:
        bed_align='results/depth/sample/{sample}/reads/aligned_reads_{index}.bed.gz',
        bed_win='results/depth/data/windows_{win_size}.bed.gz'
    output:
        bed=temp('temp/depth/sample/{sample}/depth/depth_{win_size}_{mapq}/{index}.bed.gz')
    run:

        # Get minimum MAPQ
        try:
            mapq = int(wildcards.mapq)
        except ValueError:
            raise RuntimeError('Non-integer MAPQ: "{}" (wildcards.mapq)'.format(wildcards.mapq))

        if mapq < 0:
            raise RuntimeError('MAPQ may not be less than 0: {} (wildcards.mapq)'.format(wildcards.mapq))

        # Definitions
        last_chrom = None
        chrom_tree = None

        # Get a table of reads
        df = pd.read_csv(input.bed_align, sep='\t', header=0)
        df = df.loc[df['MAPQ'] >= mapq]

        # Traverse windows
        with gzip.open(input.bed_win, 'r') as in_file:
            with gzip.open(output.bed, 'w') as out_file:
                out_file.write('#CHROM\tPOS\tEND\tDEPTH\n'.encode())

                for line in in_file:

                    # Decode line
                    line = line.decode('utf-8').strip()

                    if not line or line.startswith('#'):
                        continue

                    chrom, pos, end = line.split('\t')
                    pos = int(pos)
                    end = int(end)

                    # Setup tree
                    if chrom != last_chrom:
                        chrom_tree = intervaltree.IntervalTree()

                        for index, row in df.loc[df['#CHROM'] == chrom].iterrows():
                            chrom_tree[row['POS']:row['END']] = 1

                        last_chrom = chrom

                    # Count bp
                    bp_count = 0

                    for interval in chrom_tree[pos:end]:
                        bp_count += min(interval.end, end) - max(interval.begin, pos)

                    # Write
                    out_file.write(
                        '{}\t{:d}\t{:d}\t{:0.4f}\n'.format(
                            chrom, pos, end,
                            bp_count / (end - pos)
                        ).encode()
                    )


#
# Alignment record BED
#

# depth_data_sample_read_bed
#
# Get alignment BED for one part (one aligned cell or split BAM) in one sample.
rule depth_data_sample_read_bed:
    input:
        tab='results/depth/sample/{sample}/tables/bam_table.tab'
    output:
        bed='results/depth/sample/{sample}/reads/aligned_reads_{index}.bed.gz'
    run:

        # Get alignment file name
        align_file_name = pd.read_csv(input.tab, sep='\t', index_col='INDEX', squeeze=True)[int(wildcards.index)]

        # Get records
        clip_l = 0
        clip_r = 0

        record_list = list()

        if os.stat(align_file_name).st_size > 100:
            with pysam.AlignmentFile(align_file_name, 'rb') as in_file:
                for record in in_file:

                    # Skipped unmapped reads
                    if record.is_unmapped:
                        continue

                    # Read tags
                    tags = dict(record.get_tags())

                    # Get clipping
                    cigar_tuples = record.cigartuples

                    l_index = 0 if cigar_tuples[0][0] != 5 else 1
                    r_index = -1 if cigar_tuples[-1][0] != 5 else -2

                    clip_l = cigar_tuples[l_index][1] if cigar_tuples[l_index][0] == 4 else 0
                    clip_r = cigar_tuples[r_index][1] if cigar_tuples[r_index][0] == 4 else 0

                    # Save record
                    record_list.append(
                    pd.Series(
                        [
                            record.reference_name,
                            record.reference_start,
                            record.reference_end,
                            record.query_name,

                            tags['RG'] if 'RG' in tags else 'NA',

                            record.mapping_quality,
                            clip_l,
                            clip_r,

                            str(record.is_reverse),
                            '0x{:04x}'.format(record.flag)
                        ],
                        index=[
                            '#CHROM', 'POS', 'END', 'ID',
                            'RG',
                            'MAPQ', 'CLIP_L', 'CLIP_R',
                            'REV', 'FLAGS'
                        ]
                    ))

            # Merge records
            df = pd.concat(record_list, axis=1).T

        else:
            df = pd.DataFrame(
                [],
                columns=[
                    '#CHROM', 'POS', 'END', 'ID',
                    'RG',
                    'MAPQ', 'CLIP_L', 'CLIP_R',
                    'REV', 'FLAGS'
                ]
            )

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


#
# Get a table of input alignment files per sample
#

# depth_data_sample_bam_table
#
# Make a table with each aligned file (FILE) and assign a unique numeric ID (INDEX) to each.
rule depth_data_sample_bam_table:
    output:
        tab='results/depth/sample/{sample}/tables/bam_table.tab'
    run:

        # Get sample table entry.
        sample_params = variant_global_sample_table_entry('smrtsv', wildcards.sample)

        if not isinstance(sample_params, pd.Series):
            raise RuntimeError('Configuration for SMRT-SV sample {} is not a Series'.format(wildcards.sample))

        # Check version
        if sample_params['VERSION'] != "2":
            raise RuntimeError('Cannot get depth for SMRT-SV sample with version {}: Currently only supports version 2'.format(sample_params['VERSION']))

        # Get base path
        align_fofn = os.path.join(sample_params['DATA'], 'align/alignments.fofn')

        # Get a list of files
        with open(os.path.join(sample_params['DATA'], 'align/alignments.fofn'), 'r') as in_file:
            file_list = [line.strip() for line in in_file if line.strip()]

        # Get a dataframe with two columns (INDEX and FILE) where the index is sequential starting from 0
        df = pd.concat(
            [pd.Series([index, file], index=['INDEX', 'FILE']) for index, file in zip(range(len(file_list)), file_list)],
            axis=1
        ).T

        # Write
        df.to_csv(output.tab, sep='\t', index=False)


#
# Make windows
#

# depth_data_make_windows
#
# Make windows for alignment depth.
rule depth_data_make_windows:
    input:
        fai=config['reference'] + '.fai'
    output:
        bed='results/depth/data/windows_{win_size}.bed.gz'
    shell:
        """{{\n"""
        """    echo -e "#CHROM\tPOS\tEND"; \n"""
        """    bedtools makewindows -g {input.fai} -w {wildcards.win_size} \n"""
        """ }} | """
        """gzip """
        """>{output.bed}"""
