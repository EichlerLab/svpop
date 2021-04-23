"""
Find duplication donor sites for insertions.
"""

# variant_caller_anno_altmap_sv_distance
#
# Filter insertion alignment sites by distance.
rule variant_caller_anno_altmap_sv_distance:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/altmap/altmap-{mapper}_sv_ins.bed.gz'
    output:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/altmap/altmap-{mapper}-distance-{dist_prop}-{ro}-{map_prop}_sv_ins.bed.gz'
    wildcard_constraints:
        dist_prop='\\d+(\\.\\d+)?',
        #ro='\\d+|(any)',
        #map_prop='\\d+|(any)',
        mapper='minimap2'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Filter by SV call distance to mapped location
        if wildcards.dist_prop.lower() != 'any':

            # Filter for same chromosome first
            df = df.loc[~ pd.isnull(df['DISTANCE'])]

            # Filter by distance
            if wildcards.dist_prop.lower() != 'samechr':
                df = df.loc[df['DISTANCE_PROP'] <= float(wildcards.dist_prop)]

        # Filter by length reciprocal-overlap
        if wildcards.ro.lower() != 'any':
            df = df.loc[df['LEN_RO'] <= float(wildcards.ro) / 100]

        # Filter by proportion of the INS matching the alignment
        if wildcards.map_prop.lower() != 'any':
            df = df.loc[df['MATCH_BP'] / df['SVLEN'] <= float(wildcards.map_prop) / 100]

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# variant_caller_anno_altmap_bed
#
# Make BED of INS alignments.
#
# Note: Requires "=" and "X" CIGAR operations.
rule variant_caller_anno_altmap_bed:
    input:
        sv_bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        bam='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/altmap/altmap-{mapper}_{vartype}_{svtype}.bam'
    output:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/altmap/altmap-{mapper}_{vartype}_{svtype}.bed.gz'
    wildcard_constraints:
        mapper='minimap2',
        vartype='\w+',
        svtype='\w+'
    run:

        bed_record_list = list()
        record_count = 0

        # Define headers for writing empty files for 0 alignments.
        header_list = [
            '#CHROM', 'POS', 'END', 'ID', 'MAPQ',
            'FLAGS', 'IS_REV',
            'MATCH_BP', 'MISMATCH_BP', 'INS_BP', 'INS_N', 'DEL_BP', 'DEL_N',
            'CLIPPED_BP', 'CLIPPED_N',
            'MAP_LEN', 'SV_CHROM',
            'SV_POS', 'SV_END', 'SVLEN', 'LEN_RO',
            'DISTANCE', 'DISTANCE_PROP'
        ]

        # Write no records in input BAM is empty
        if os.stat(input.bam).st_size == 0:

            with open(output.bed, 'w') as out_file:
                out_file.write('\t'.join(header_list) + '\n')

            return

        # Convert BAM records
        with pysam.AlignmentFile(input.bam, 'r') as in_file:
            for record in in_file:
                record_count += 1

                cigar_stat_bp, cigar_stat_n = record.get_cigar_stats()

                if cigar_stat_bp[0] > 0:
                    raise RuntimeError(
                        'Detected match records (M) in CIGAR string for record {}: Expected only "=" and "X" CIGAR operations'.format(record_count)
                    )

                if not record.is_unmapped:
                    bed_record_list.append(pd.Series(
                        [
                            record.reference_name,
                            record.reference_start,
                            record.reference_start + record.reference_length,
                            record.query_name,
                            record.mapq,
                            record.flag,
                            record.is_reverse,
                            cigar_stat_bp[CIGAR_INT['=']],
                            cigar_stat_bp[CIGAR_INT['X']],
                            cigar_stat_bp[CIGAR_INT['I']],
                            cigar_stat_n[CIGAR_INT['I']],
                            cigar_stat_bp[CIGAR_INT['D']],
                            cigar_stat_n[CIGAR_INT['D']],
                            cigar_stat_bp[CIGAR_INT['S']] + cigar_stat_bp[CIGAR_INT['H']],
                            cigar_stat_n[CIGAR_INT['S']] + cigar_stat_n[CIGAR_INT['H']],
                            record.reference_length
                        ],
                        index=[
                            '#CHROM', 'POS', 'END', 'ID', 'MAPQ',
                            'FLAGS', 'IS_REV',
                            'MATCH_BP', 'MISMATCH_BP',
                            'INS_BP', 'INS_N', 'DEL_BP', 'DEL_N',
                            'CLIPPED_BP', 'CLIPPED_N', 'MAP_LEN'
                        ]
                    ))

        # Write an emtpy file if there are no records
        if len(bed_record_list) == 0:
            with open(output.bed, 'w') as out_file:
                out_file.write('\t'.join(header_list) + '\n')

            return

        # Merge
        df = pd.concat(bed_record_list, axis=1).T

        df.sort_values(['#CHROM', 'POS'], inplace=True)

        df.set_index('ID', inplace=True, drop=False)

        # Annotate with SVs
        if any([val > 1 for val in collections.Counter(df['ID']).values()]):
            raise RuntimeError('Alignment BED contains multiple records for a single SV')

        df_sv = pd.read_csv(input.sv_bed, sep='\t', header=0, usecols=('#CHROM', 'POS', 'END', 'ID', 'SVLEN'), index_col='ID')
        df_sv = df_sv.loc[df['ID']]

        df['SV_CHROM'] = df_sv['#CHROM']
        df['SV_POS'] = df_sv['POS']
        df['SV_END'] = df_sv['END']
        df['SVLEN'] = df_sv['SVLEN']

        # Annotate distance
        df['LEN_RO'] = df.apply(lambda row: np.min([row['MAP_LEN'] / row['SVLEN'], row['SVLEN'] / row['MAP_LEN']]), axis=1)

        def sv_align_distance(row):
            if row['#CHROM'] != row['SV_CHROM']:
                return np.nan

            if row['POS'] < row['SV_POS']:
                return row['SV_POS'] - row['POS']
            elif row['END'] > row['SV_END']:
                return row['END'] - row['SV_END']
            else:
                return 0

        df['DISTANCE'] = df.apply(sv_align_distance, axis=1).astype(np.float32)

        df['DISTANCE_PROP'] = df['DISTANCE'] / df_sv['SVLEN']

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

#
# Minimap2
#

# variant_caller_anno_altmap_minimap2_sort
#
# Sort BAM.
rule variant_caller_anno_altmap_minimap2_sort:
    input:
        bam='temp/variant/caller/{sourcename}/{sample}/{filter}/all/anno/altmap/altmap-minimap2_{vartype}_{svtype}/align.bam'
    output:
        bam='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/altmap/altmap-minimap2_{vartype}_{svtype}.bam'
    shell:
        """if [ -s {input.bam} ]; then """
            """SORT_TEMP=$(dirname {input.bam})/donor_align_sort; """
            """samtools view -h {input.bam} | """
            """samtools sort -T ${{SORT_TEMP}} -o {output.bam}; """
            """rm -f ${{SORT_TEMP}}*; """
            """sleep 10; """
            """samtools index {output.bam}; """
        """else """
            """ touch {output.bam}; """
        """fi"""

# variant_caller_anno_altmap_minimap2_map
#
# Map insertions to the reference genome.
rule variant_caller_anno_altmap_minimap2_map:
    input:
        fa='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz',
        ref=config['reference']
    output:
        bam=temp('temp/variant/caller/{sourcename}/{sample}/{filter}/all/anno/altmap/altmap-minimap2_{vartype}_{svtype}/align.bam')
    shell:
        """if [ -s {input.fa} ]; then """
        """minimap2 """
            """-x asm20 -H --secondary=no -r 2k -Y -a --eqx -L -t 4 """
            """{input.ref} {input.fa} | """
        """samtools view """
            """-h -F 0x800 -O BAM -o {output.bam}; """
        """else """
            """ touch {output.bam}; """
        """fi"""
