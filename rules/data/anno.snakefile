"""
Prepare annotation input files.
"""


#############
### Rules ###
#############

#
# Region Merge
#

# data_merge_bed_regions
#
# Merge all BED records within a set distance to produce regions covered by the BED. First, regions are merged
# by "distance", then the merged regions are expanded by "flank". If any regions overlap after expanding, then they
# are merged.
rule data_merge_bed_regions:
    input:
        bed='data/anno/{annotype}/{annoname}.bed.gz',
        fai=config['reference'] + '.fai'
    output:
        bed='data/anno/{annotype}/{annoname}_regions_{distance}_{flank}.bed.gz'
    wildcard_constraints:
        distance=r'\d+',
        flank=r'\d+'
    shell:
        """awk -vOFS="\\t" '{{print $1, $2, $3}}' <(zcat {input.bed}) | """
        """python3 {SVPOP_DIR}/scripts/filter_chrom.py -g {input.fai} | """
        """bedtools merge -d {wildcards.distance} -header | """
        """bedtools slop -b {wildcards.flank} -g {input.fai} -header | """
        """awk '$2 >= 0' | """
        """bedtools merge -header | """
        """gzip > {output.bed}"""

#
# ORegAnno
#

# data_oreganno
#
# Prepare ORegAnno BED.
rule data_oreganno:
    input:
        txt='temp/data/anno/ucsc/database/oreganno.txt.gz'
    output:
        bed='data/anno/oreganno/oreganno.bed.gz'
    shell:
        """{{\n"""
        """    echo -e "#CHROM\tPOS\tEND\tID\tSTRAND\tNAME";\n"""
        """    awk -vOFS="\t" '{{print $2, $3, $4, $5, $6, $7}}' <(zcat {input.txt});\n"""
        """}} | gzip > {output.bed}"""

# data_oreganno_table_gz
#
# Compress ORegAnno annotation table.
rule data_oreganno_table_gz:
    input:
        tsv='temp/data/anno/oreganno/oreganno-table.tsv'
    output:
        tsv='data/anno/oreganno/oreganno-table.tsv.gz'
    shell:
        """gzip -c {input.tsv} > {output.tsv}"""

# data_oreganno_table
#
# Get a table of oreganno annotations.
rule data_oreganno_table:
    output:
        tsv=temp('temp/data/anno/oreganno/oreganno-table.tsv')
    params:
        url=config.get('data', dict()).get('oreganno', dict()).get('path', None)
    run:

        if params.url is None:
            raise RuntimeError('Missing ORegAnno table path in config["data"]["oreganno"]["path"]')

        shell("""wget {params.url} -O {output.tsv}""")

#
# ENCODE 2020 (DHS: Vierstra 2020, CCRE: Consortium 2020)
#

# data_encode_dhs2020_minscore
#
# Get table with minimum score
rule data_encode_dhs2020_minscore:
    input:
        bed='data/anno/dhs2020/dhs_all.bed.gz'
    output:
        bed='data/anno/dhs2020/dhs_{score}.bed.gz'
    wildcard_constraints:
        score='[0-9]+'
    run:

        min_score = int(wildcards.score)

        df = pd.read_csv(input.bed, sep='\t')

        df = df.loc[df['DHS_SIG_MEAN'] >= min_score]

        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# data_encode_dhs2020_header
#
# Add header to DHS 2020
#
# Legend:
#   https://resources.altius.org/~jvierstra/projects/footprinting.2020/consensus.index/consensus_footprints_and_motifs_legend.txt
rule data_encode_ccre2020_header:
    input:
        bed='temp/data/anno/ccre2020/ccre2020.bed'
    output:
        bed='data/anno/ccre2020/ccre_all.bed.gz'
    run:

        df = pd.read_csv(input.bed, sep='\t', header=None)

        df.columns = [
            '#CHROM', 'POS', 'END',
            'EH38D', 'EH38E',
            'CCRE2020'
        ]

        df.sort_values(['#CHROM', 'POS'], inplace=True)

        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# data_encode_dhs2020_dl
#
# Download DHS data.
rule data_encode_ccre2020_dl:
    output:
        bed=temp('temp/data/anno/ccre2020/ccre2020.bed')
    shell:
        """wget """
            """http://gcp.wenglab.org/GRCh38-ccREs.bed """
            """-O {output.bed}"""

# data_encode_dhs2020_header
#
# Add header to DHS 2020
#
# Legend:
#   https://resources.altius.org/~jvierstra/projects/footprinting.2020/consensus.index/consensus_footprints_and_motifs_legend.txt
rule data_encode_dhs2020_header:
    input:
        bed='temp/data/anno/dhs2020/dhs2020.bed.gz'
    output:
        bed='data/anno/dhs2020/dhs_all.bed.gz'
    run:

        df = pd.read_csv(input.bed, sep='\t', header=None)

        df.columns = [
            '#CHROM', 'POS', 'END',
            'DHS_ID',
            'DHS_SIG_MEAN', 'SAMPLE_N', 'FPS_N',
            'WIDTH', 'SUMMIT',
            'CORE_START', 'CORE_END',
            'MOTIF_CLUSTERS'
        ]

        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# data_encode_dhs2020_dl
#
# Download DHS data.
rule data_encode_dhs2020_dl:
    output:
        bed='temp/data/anno/dhs2020/dhs2020.bed.gz'
    shell:
        """wget """
            """https://resources.altius.org/~jvierstra/projects/footprinting.2020/consensus.index/consensus_footprints_and_motifs_hg38.bed.gz """
            """-O {output.bed}"""

#
# Encode DHS (DNAse hypersensitivity) - UCSC browser
#

# data_encode_dhs_cluster_score
#
# Subset DHS by minimum score.
rule data_encode_dhs_cluster_score:
    input:
        bed='data/anno/dhs_cluster/dhs_cluster_all.bed.gz'
    output:
        bed='data/anno/dhs_cluster/dhs_cluster_{score,[^(all)].*}.bed.gz'
    run:

        min_score = int(wildcards.score)

        df = pd.read_csv(input.bed, sep='\t', header=0)

        df = df.loc[df['DHS_SCORE'] >= min_score]

        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# data_encode_dhs_cluster
#
# Get DHS cluster BED.
rule data_encode_dhs_cluster:
    input:
        txt='temp/data/anno/ucsc/database/wgEncodeRegDnaseClustered.txt.gz'
    output:
        bed='data/anno/dhs_cluster/dhs_cluster_all.bed.gz'
    run:

        df = pd.read_csv(
            input.txt, sep='\t', header=None,
            names=('bin', '#CHROM', 'POS', 'END', 'DHS_NAME', 'DHS_SCORE', 'DHS_SOURCE_COUNT', 'DHS_SOURCE_IDS', 'DHS_SOURCE_SCORES'),
            usecols=('#CHROM', 'POS', 'END', 'DHS_NAME', 'DHS_SCORE', 'DHS_SOURCE_COUNT', 'DHS_SOURCE_IDS', 'DHS_SOURCE_SCORES')
        )

        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


#
# Encode histone modification (H3K4Me1, H3K4Me3, H3K27Ac)
#

# data_encode_merge_cells
#
# Merge all cells.
rule data_encode_merge_cells:
    input:
        bed=expand(
            'data/anno/encode/threshold/{{mark}}_{cell_line}_{{threshold}}.bed.gz',
            cell_line=('GM12878', 'H1-hESC', 'HSMM', 'HUVEC', 'K562', 'NHEK', 'NHLF')
        )
    output:
        bed='data/anno/encode/threshold/{mark}_all_{threshold}.bed.gz'
    shell:
        """{{\n"""
        """    head -n 1 <(zcat {input.bed[0]});\n"""
        """    cat <(zcat {input.bed}) | egrep -v "^#" | sort -k 1,1 -k2,2n;\n"""
        """}} | """
        """bedtools merge -i stdin -header | """
        """gzip > {output.bed}"""

# data_encode_threshold_merge
#
# Filter by threshold and merge windows within 100 bp.
rule data_encode_threshold_merge:
    input:
        bed='data/anno/encode/full/{mark}_{cell_line}_all.bed.gz'
    output:
        bed='data/anno/encode/threshold/{mark}_{cell_line,[^(all)].*}_{threshold,[^(all)].*}.bed.gz'
    shell:
        """{{\n"""
        """    zcat {input.bed} | \n"""
        """    awk -vOFS="\\t" '(NR == 1 || $4 >= {wildcards.threshold}) {{print $1, $2, $3}}' | \n"""
        """    bedtools merge -i stdin -d 100 -header;"""
        """}} | gzip > {output.bed}"""

# data_encode_bw_to_bed
#
# BigWig to BedGraph for encode.
rule data_encode_bw_to_bed:
    input:
        bw='temp/data/anno/encode/bw/encode_{mark}_{cell_line}.bigWig'
    output:
        bed='data/anno/encode/full/{mark}_{cell_line}_all.bed.gz',
        bed_tmp=temp('temp/data/anno/encode/full/{mark}_{cell_line}_all.bed.gz')
    shell:
        """bigWigToBedGraph {input.bw} {output.bed_tmp}; """
        """{{\n"""
        """    echo -e "#CHROM\tPOS\tEND\tENCODE_SIGNAL";\n"""
        """    cat {output.bed_tmp}\n"""
        """}} | gzip > {output.bed}"""

# data_encode_dl
#
# Download Encode data
rule data_encode_dl:
    output:
        bw=temp('temp/data/anno/encode/bw/encode_{mark}_{cell_line}.bigWig')
    run:

        # Get "PascalCase" version of the cell line name and encode mark
        cell_line_pascal = wildcards.cell_line.lower().replace('-', '')
        cell_line_pascal = cell_line_pascal[0].upper() + cell_line_pascal[1:]

        mark_pascal = wildcards.mark.lower().replace('-', '')
        mark_pascal = mark_pascal[0].upper() + mark_pascal[1:]

        shell(
            """rsync rsync://hgdownload.cse.ucsc.edu:/gbdb/{UCSC_REF_NAME}/bbi/wgEncodeReg/wgEncodeRegMark{mark_pascal}/wgEncodeBroadHistone{cell_line_pascal}{mark_pascal}StdSig.bigWig {output.bw}"""
        )

#
# CpG Islands
#

# data_cpg_unmasked
#
# CpG unmasked (includes records in RepeatMasked regions).
rule data_cpg_unmasked:
    input:
        txt='temp/data/anno/ucsc/database/cpgIslandExtUnmasked.txt.gz'
    output:
        bed='data/anno/cpgisland/cpgisland_unmasked.bed.gz'
    shell:
        """{{\n"""
        """    echo -e "#CHROM\tPOS\tEND\tNAME\tLENGTH\tCPG_NUM\tGC_NUM\tCPG_PCT\tGT_PCT\tOBS_EXP";\n"""
        """    awk -vOFS="\t" '{{print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}}' <(zcat {input.txt});\n"""
        """}} | gzip > {output.bed}"""


#
# RepeatMasker
#

# An alias of pattern keywords for RMSK filters. Each value is a tuple [0] pattern, and [1] a boolean indicating if the
# pattern is a plain string or a regular expression.
_ANNO_RMSK_FAM_PATTERN = {
    'any': None,
    'simple': {
        'class': r'Simple_repeat'
    },
    'low-complexity': {
        'class': r'Low_complexity'
    },
    'satellite': {
        'class': r'Satellite'
    },
    'alu': {
        'fam': r'ALU'
    },
    'l1': {
        'fam': r'L1'
    },
    'l1hs': {
        'name': r'L1HS'
    },
    'sva': {
        'fam': r'SVA'
    },
    'aluy': {
        'name': r'AluY.*'
    }
}

# data_rmsk_subset_bed
#
# Subset the RepeatMasker BED file.
rule data_rmsk_subset_bed:
    input:
        bed='data/anno/rmsk/full_table/rmsk_all.bed.gz'
    output:
        bed='data/anno/rmsk/rmsk-{filter_spec}-{ident}.bed.gz'
    wildcard_constraints:
        filter_spec='[a-z0-9]+',
        ident='((lt|gt)[\d]+)|([\d]+to[\d]+)|all'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t')

        # Get identity (%) filter
        if re.match('^gt[\d]+$', wildcards.ident):
            min_ident = int(wildcards.ident[2:])
            max_ident = None

        elif re.match('^lt[\d]+$', wildcards.ident):
            min_ident = None
            max_ident = int(wildcards.ident[2:])

        elif re.match('^[\d]+to[\d]+$', wildcards.ident):

            match_obj = re.match('^([\d]+)to([\d]+)$', wildcards.ident)

            min_ident = match_obj[1]
            max_ident = match_obj[2]

        elif wildcards.ident == 'all':
            min_ident = None
            max_ident = None

        else:
            raise RuntimeError('Unknown identity filter specification: ' + wildcards.ident)

        # Do filter
        if min_ident is not None:
            if min_ident > 100:
                raise RuntimeError('Minimum identity must not be greater than 100: ' + wildcards.ident)

            if max_ident is not None and min_ident > max_ident:
                raise RuntimeError('Minimum identity must not be greater than maximum identity: ' + wildcards.ident)

            df = df.loc[df['MILLI_DIV'] >= min_ident * 10]

        if max_ident is not None:
            if max_ident > 100:
                raise RuntimeError('Maximum identity must not be greater than 100: ' + wildcards.ident)

            df = df.loc[df['MILLI_DIV'] <= min_ident * 10]

        if df.shape[0] == 0:
            raise RuntimeError('No RMSK records after filtering by identity: {}'.format(wildcards.ident))

        # Get filter specification
        if wildcards.filter_spec not in _ANNO_RMSK_FAM_PATTERN:
            raise RuntimeError('Unknown repeat filter_spec: {}: (May be missing from _ANNO_RMSK_FAM_PATTERN)'.format(wildcards.filter_spec))

        filter_spec = _ANNO_RMSK_FAM_PATTERN[wildcards.filter_spec]

        # Do filter
        if filter_spec is not None:  # None is set if filter_spec == "any"

            # Check filter_spec
            if len(set(filter_spec.keys())) == 0:
                raise RuntimeError('No filter spec keywords for filter_spec {}'.format(wildcards.filter_spec))

            unknown_spec_keyword = sorted(set(filter_spec.keys()) - {'class', 'fam', 'name'})

            if unknown_spec_keyword:
                raise RuntimeError('Unknown filter specification keywords for filter_spec {}: {}'.format(wildcards.filter_spec, ', '.join(unknown_spec_keyword)))

            # Process each keyword in filter_spec
            for spec_kw in ('class', 'fam', 'name'):

                # Skip unspecified filter specifications
                if spec_kw not in filter_spec.keys():
                    continue

                # Do filter
                re_pattern = re.compile(r'^' + filter_spec[spec_kw] + r'$')

                df = df.loc[df['REP_{}'.format(spec_kw.upper())].apply(lambda val: bool(re.match(re_pattern, val)))]

        if df.shape[0] == 0:
            raise RuntimeError('No RMSK records after filtering by filter_spec: {}'.format(wildcards.filter_spec))

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')



# data_rmsk_to_bed
#
# Make RepeatMasker annotation BED file.
rule data_rmsk_to_bed:
    input:
        txt='temp/data/anno/ucsc/database/rmsk.txt.gz'
    output:
        bed='data/anno/rmsk/full_table/rmsk_all.bed.gz'
    run:

        # Read
        df = pd.read_csv(
            input.txt,
            sep='\t',
            header=None,
            names=(
                'bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns', 'genoName', 'genoStart', 'genoEnd', 'genoLeft',
                'strand', 'repName', 'repClass', 'repFamily', 'repStart', 'repEnd', 'repLeft', 'id'
            ),
            usecols=('genoName', 'genoStart', 'genoEnd', 'swScore', 'strand', 'repClass', 'repFamily', 'repName', 'milliDiv', 'milliDel', 'milliIns')
        )

        # Assign unique ID
        df['ID'] = df.apply(lambda row: '{}-{}-{}-{}'.format(row['genoName'], row['genoStart'] + 1, row['repClass'], row['genoEnd'] - row['genoStart']), axis=1)

        # Rearrange columns
        df = df.loc[:, ('genoName', 'genoStart', 'genoEnd', 'ID', 'swScore', 'strand', 'repClass', 'repFamily', 'repName', 'milliDiv', 'milliDel', 'milliIns')]
        df.columns = ('#CHROM', 'POS', 'END', 'ID', 'SCORE', 'STRAND', 'REP_CLASS', 'REP_FAM', 'REP_NAME', 'MILLI_DIV', 'MILLI_DEL', 'MILLI_INS')

        df.sort_values(['#CHROM', 'POS', 'END', 'ID'], inplace=True)

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


#
# Segmental Duplications
#

# data_sd_max_merge
#
# Concatenate consecutive records if they have the same identity.
rule data_sd_max_merge:
    input:
        bed='temp/data/anno/sd/sd-max-{match_type}.bed.gz'
    output:
        bed='data/anno/sd/sd-max-{match_type}.bed.gz'
    wildcard_constraints:
        match_type='frac|fracindel'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', names=['#CHROM', 'POS', 'END', 'MATCH'])

        # Merge records
        record_list = list()

        chrom = None
        pos = 0
        end = 0
        match = 0

        for index, row in df.iterrows():

            # Switch chromosomes
            if row['#CHROM'] != chrom:

                # Write last chromosome record
                if chrom is not None:

                    record_list.append(pd.Series(
                        [chrom, pos, end, match],
                        index=['#CHROM', 'POS', 'END', 'MATCH']
                    ))

                # Initialize for this chromosome with the first record
                chrom = row['#CHROM']

                pos = row['POS']
                end = row['END']
                match = row['MATCH']

                # Skip
                continue

            # Extend or new record
            if row['POS'] == end - 1 and row['MATCH'] == match:

                # Extend last record
                end = row['END']

            else:

                # Write last record
                record_list.append(pd.Series(
                    [chrom, pos, end, match],
                    index=['#CHROM', 'POS', 'END', 'MATCH']
                ))

                # Start new record
                pos = row['POS']
                end = row['END']
                match = row['MATCH']

        # Write final record
        if chrom is not None:

            record_list.append(pd.Series(
                [chrom, pos, end, match],
                index=['#CHROM', 'POS', 'END', 'MATCH']
            ))

        # Merge and write
        df_merge = pd.concat(record_list, axis=1).T

        df_merge.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# data_sd_max
#
# Get maximum identity per window.
rule data_sd_max:
    input:
        bed='data/anno/sd/sd.bed.gz',
        frag_bed='temp/data/anno/sd/sd-max-fragmented.bed.gz'
    output:
        bed=temp('temp/data/anno/sd/sd-max-{match_type}.bed.gz')
    params:
        frac_field=lambda wildcards: 'FRAC_MATCH' if wildcards.match_type == 'frac' else ('FRAC_MATCH_INDEL' if wildcards.match_type == 'fracindel' else 'ERR-UNKNOWN-FIELD')
    wildcard_constraints:
        match_type='(frac|fracindel)'
    shell:
        """zcat {input.bed} | awk -vOFS="\\t" '"""
            """NR==1 {{ for (i=1; i<=NF; i++) {{ f[$i] = i }} }} """
            """{{ print $(f["#CHROM"]), $(f["POS"]), $(f["END"]), $(f["{params.frac_field}"]) }} """
        """' | """
        """bedtools map -a {input.frag_bed} -b stdin -c 4 -o max | """
        """gzip """
        """> {output.bed}"""

# data_sd_max_fragment
#
# Fragment SD records to non-overlapping segments.
rule data_sd_max_fragment:
    input:
        bed='data/anno/sd/sd.bed.gz'
    output:
        bed=temp('temp/data/anno/sd/sd-max-fragmented.bed.gz')
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', low_memory=False)

        # Definitions
        def split_list_by_region(target_list, pos, end):

            in_region = list()
            out_region = list()

            for next_pos in target_list:
                if pos <= next_pos < end:
                    in_region.append(next_pos)
                else:
                    out_region.append(next_pos)

            return in_region, out_region


        # Open output file
        with gzip.open(output.bed, 'wt') as out_file:

            # Process each chromosome
            for chrom in sorted(set(df['#CHROM'])):

                df_chrom = df.loc[df['#CHROM'] == chrom]

                pos_list = sorted(set(df_chrom['POS']) | set(df_chrom['END'] - 1))

                # Find max regions
                max_region_list = list()

                start_pos = 0
                max_end = 0

                for index, row in df_chrom.iterrows():

                    if row['POS'] > max_end:

                        # List last record
                        if start_pos < max_end:
                            max_region_list.append((start_pos, max_end))

                        start_pos = row['POS']
                        max_end = row['END']

                    else:
                        max_end = np.max([max_end, row['END']])

                # List last max record
                if start_pos < max_end:
                    max_region_list.append((start_pos, max_end))

                # Write breaks by max records
                for max_pos, max_end in max_region_list:
                    this_pos_list, pos_list = split_list_by_region(pos_list, max_pos, max_end)

                    pos_start = this_pos_list[0]

                    for pos_end in this_pos_list[1:]:
                        out_file.write('{}\t{}\t{}\n'.format(chrom, pos_start, pos_end + 1))

                        pos_start = pos_end

# data_sd_to_bed
#
# SD track to BED.
rule data_sd_to_bed:
    input:
        txt='temp/data/anno/ucsc/database/genomicSuperDups.txt.gz'
    output:
        bed='data/anno/sd/sd.bed.gz'
    shell:
        """{{\n"""
        """    echo -e "#CHROM\\tPOS\\tEND\\tNAME\\tSTRAND\\tOTHER_CHROM\\tOTHER_POS\\tOTHER_END\\tOTHER_SIZE\\tFRAC_MATCH\\tFRAC_MATCH_INDEL"; \n"""
        """    awk -vOFS="\\t" '{{print $2, $3, $4, $5, $7, $8, $9, $10, $11, $27, $28}}' <(zcat {input.txt}) \n"""
        """}} | gzip > {output.bed}"""


#
# TRF
#

# data_get_trf_txt_to_bed
#
# TRF to BED.
rule data_get_trf_txt_to_bed:
    input:
        txt='temp/data/anno/ucsc/database/simpleRepeat.txt.gz'
    output:
        bed='data/anno/trf/trf.bed.gz'
    shell:
        """{{\n"""
        """    echo -e "#CHROM\tPOS\tEND";\n"""
        """    awk 'BEGIN {{OFS="\\t"}} {{print $2, $3, $4}}' <(zcat {input.txt});\n"""
        """}} | gzip > {output.bed}"""


#
# RefSeq
#

# data_anno_refseq_updown_bed
#
# Get a BED file of regions upstream and downstream of refseq annotations.
rule data_anno_refseq_updown_bed:
    input:
        bed='data/anno/refseq/refseq.bed.gz'
    output:
        bed_up='data/anno/refseq/region/flank_{flank}/refseq_up.bed.gz',
        bed_dn='data/anno/refseq/region/flank_{flank}/refseq_dn.bed.gz'
    run:

        flank = int(wildcards.flank)

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0, usecols=('#chrom', 'txStart', 'txEnd', 'name', 'strand', 'comName'))

        # Make upstream table
        df_up = df.copy()
        df_up['txStart'] = df_up.apply(lambda row: row['txStart'] - flank if row['strand'] == '+' else row['txEnd'], axis=1)
        df_up['txEnd'] = df_up.apply(lambda row: row['txStart'] + flank if row['strand'] == '+' else row['txStart'] + flank, axis=1)

        df_up.loc[df_up['txStart'] < 0, 'txStart'] = 0

        df_up = df_up.loc[df_up.apply(lambda row: (row['txEnd'] - row['txStart']) > 0, axis=1)]  # Remove start=0 && end=0 records

        df_up.sort_values(['#chrom', 'txStart'], inplace=True)

        # Make downstream table
        df_dn = df.copy()
        df_dn['txStart'] = df_dn.apply(lambda row: row['txEnd'] if row['strand'] == '+' else row['txStart'] - flank, axis=1)
        df_dn['txEnd'] = df_dn.apply(lambda row: row['txStart'] + flank if row['strand'] == '+' else row['txStart'] + flank, axis=1)

        df_dn.loc[df_dn['txStart'] < 0, 'txStart'] = 0

        df_dn = df_dn.loc[df_dn.apply(lambda row: (row['txEnd'] - row['txStart']) > 0, axis=1)]  # Remove start=0 && end=0 records

        df_dn.sort_values(['#chrom', 'txStart'], inplace=True)

        # Write
        df_up.to_csv(output.bed_up, sep='\t', index=False, compression='gzip')
        df_dn.to_csv(output.bed_dn, sep='\t', index=False, compression='gzip')

# data_anno_refseq_bed
#
# RefSeq BED.
rule data_anno_refseq_bed:
    input:
        refseq='temp/data/anno/ucsc/database/refGene.txt.gz'
    output:
        bed='data/anno/refseq/refseq.bed.gz'
    shell:
        """{{\n"""
            """echo -e "#chrom\ttxStart\ttxEnd\tname\tstrand\tcomName\tcdsStart\tcdsEnd\texonStarts\texonEnds"; \n"""
            """awk -vOFS="\\t" '{{print $3, $5, $6, $2, $4, $13, $7, $8, $10, $11}}' <(zcat {input.refseq}) | \n"""
            """sort -k1,1 -k2,2n \n"""
        """}} | gzip > {output.bed}"""


#
# Chromosome bands
#

# data_get_chrom_bands
#
# Get chromosome bands and add a header line
rule data_get_chrom_bands:
    input:
        txt='temp/data/anno/ucsc/database/cytoBand.txt.gz'
    output:
        bed='data/anno/bands/bands.bed.gz'
    shell:
        """{{"""
        """    echo -e "#chrom\tstart\tend\tname\tgieStain";\n"""
        """    zcat {input.txt} | sort -k 1,1 -k 2,2n;\n"""
        """}} | gzip > {output.bed}"""


#
# AGP (Golden path)
#

# data_get_agp
#
# Get AGP file for the assembly. A header line is added to the original file. This shows the tiling path ("golden
# path") through scaffolds to create the primary assembly.
rule data_get_agp:
    input:
        txt=lambda wildcards: 'temp/data/anno/ucsc/bigZips/{UCSC_REF_NAME}.agp.gz'.format(UCSC_REF_NAME=UCSC_REF_NAME)
    output:
        bed='data/anno/agp/agp.bed.gz'
    shell:
        """{{\n"""
        """   echo -e "object\tobject_beg\tobject_end\tpart_number\tcomponent_type\tcomponent_id\tcomponent_beg\tcomponent_end\torientation"\n"""
        """   zcat {input.txt};\n"""
        """}} | gzip > {output.bed}; """

# variant_anno_agp_switch_bed
#
# Get a BED file of intervals around AGP switchpoints (where contigs in the reference are joined).
rule variant_anno_agp_switch_bed:
    input:
        bed='data/anno/agp/agp.bed.gz',
        fai=config['reference'] + '.fai'
    output:
        bed='data/anno/agp/agp_switch_{flank}.bed.gz'
    run:

        # NOTE: Using the BED record CHROM and START as AGP switchpoints. END will be ignored.
        # Must filter records at the start of chromosomes or that border a gap.

        # Get flank argument
        flank = int(wildcards.flank)

        if flank < 1:
            raise RuntimeError('Flank may not be less than 1: {}'.format(wildcards.flank))

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0, usecols=('object', 'object_beg', 'object_end', 'component_type', 'component_id'))
        df.columns = ('#CHROM', 'POS', 'END', 'TYPE', 'AGP_CONTIG_DOWN')

        # Annotate the downstream contig name
        df['AGP_CONTIG_UP'] = ['XXX'] + list(df['AGP_CONTIG_DOWN'][0:-1])

        # Annotate the type that came before (for filtering)
        df['TYPE_PREV'] = ['N'] + list(df['TYPE'][0:-1])

        # Filter (by gap/chr-start, then by unplaced/unlocalized)
        df = df.loc[[all(x) for x in zip(df['TYPE'] != 'N', df['TYPE_PREV'] != 'N', df['POS'] > 1)]]
        df = df.loc[~ df['#CHROM'].str.contains('_')]

        # Drop TYPE fields
        del(df['TYPE'])
        del(df['TYPE_PREV'])

        # Make copy (df is currently a slice of the original)
        df = df.copy()

        # Set END, then POS
        df['END'] = df['POS'] + flank
        df['POS'] = df['POS'] - flank

        # Read FAI (Series, length keyed by chromosome)
        fai = pd.read_csv(
            input.fai,
            sep='\t',
            names=('CHROM', 'LEN', 'START', 'LINE_BASE', 'LINE_BYTES'),
            usecols=('CHROM', 'LEN'),
            index_col='CHROM',
            squeeze=True
        )

        # Filter records within that reach the chromosome ends (should not filter anything in the reference within 500 bp flank)
        df = df.loc[df.apply(lambda row: row['POS'] > 0 and row['END'] < fai[row['#CHROM']], axis=1)]

        # Write
        df = df.loc[: , ('#CHROM', 'POS', 'END', 'AGP_CONTIG_UP', 'AGP_CONTIG_DOWN')]
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

#
# Gaps
#

# data_gap_txt_to_bed
#
# Gap to BED.
rule data_gap_txt_to_bed:
    input:
        txt='temp/data/anno/ucsc/database/gap.txt.gz'
    output:
        bed='data/anno/gap/gap.bed.gz'
    shell:
        """{{\n"""
        """    echo -e "#CHROM\tSTART\tEND\tNAME\tSIZE\tTYPE\tBRIDGE";\n"""
        """    awk -vOFS="\\t" '{{print $2, $3, $4, $2 "-" ($3 + 1) "-" $7 "-" $8 "-" $9, $7, $8, $9}}' <(zcat {input.txt});\n"""
        """}} | gzip > {output.bed}"""


#
# Centromeres
#
rule data_anno_cen:
    input:
        txt='temp/data/anno/ucsc/database/centromeres.txt.gz'
    output:
        bed='data/anno/cen/cen.bed.gz'
    run:

        df = pd.read_csv(input.txt, sep='\t', header=None, names=('BIN', '#CHROM', 'POS', 'END', 'ID'), usecols=('#CHROM', 'POS', 'END', 'ID'))

        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

#
# Get UCSC Track
#

# data_get_ucsc_track_database
#
# Download annotations from UCSC.
rule data_anno_dl_ucsc:
    output:
        txt=temp('temp/data/anno/ucsc/{subdir}/{basename}.gz')
    wildcard_constraints:
        subdir='database|bigZips'
    shell:
        """wget -P $(dirname {output.txt}) http://hgdownload.soe.ucsc.edu/goldenPath/{UCSC_REF_NAME}/{wildcards.subdir}/{wildcards.basename}.gz"""
