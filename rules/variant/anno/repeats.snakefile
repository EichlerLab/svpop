"""
Run RepeatMasker and TRF annotations on sequences.
"""

###################
### Definitions ###
###################


#############
### Rules ###
#############


#
# RepeatMasker
#

# variant_anno_repeat_rmsk_table
#
# Make RepeatMask table.
rule variant_anno_repeat_rmsk_table:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        out='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/rmsk/data/rmsk-table_{vartype}_{svtype}/rmsk.out.gz'
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/rmsk/rmsk-table_{vartype}_{svtype,ins|del|inv|dup}.tsv.gz'
    run:

        # Create empty tab file if rmsk output is empty
        if os.stat(input.out).st_size == 0:
            df = pd.DataFrame(
                [],
                columns=[
                    'ID', 'REPEAT_FAMILY', 'REPEAT_TYPE',
                    'REP_PROP', 'SV_PROP', 'ALIGN_SCORE',
                    'DIFF_PROP', 'DEL_PROP', 'INS_PROP',
                    'SV_MATCH_START', 'SV_MATCH_END', 'SV_MATCH_LEFT',
                    'STRAND', 'REP_START', 'REP_END', 'REP_LEFT'
                ]
            )

            df.to_csv(output.tab, sep='\t', index=False)

            return

        # Read records
        record_list = list()

        with gzip.open(input.out, 'rt') as in_file:
            for i in range(3):
                line = next(in_file)

            line_count = 3

            for line in in_file:
                line_count += 1

                # Split
                tok = re.split('\s+', line.strip())

                # Check number of tokens
                if len(tok) < 14:
                    raise RuntimeError('Parse error on line {}: Expected at least 14 whitespace-delimited fields: Found {}'.format(line_count, len(tok)))

                tok = tok[:14]

                # Remove parenthesis
                tok[7] = re.sub('\((.*)\)', '\\1', tok[7])
                tok[11] = re.sub('\((.*)\)', '\\1', tok[11])
                tok[13] = re.sub('\((.*)\)', '\\1', tok[13])

                # Percent to proportion (div, del, ins fields)
                tok[1] = float(tok[1]) / 100.0
                tok[2] = float(tok[2]) / 100.0
                tok[3] = float(tok[3]) / 100.0

                # Change complement to "-"
                if tok[8] == 'C':
                    tok[8] = '-'

                # Append to token list
                record_list.append(pd.Series(
                    tok,
                    index=[
                        'ALIGN_SCORE', 'DIFF_PROP', 'DEL_PROP', 'INS_PROP', 'ID', 'SV_MATCH_START', 'SV_MATCH_END',
                        'SV_MATCH_LEFT', 'STRAND', 'REPEAT_TYPE', 'REPEAT_FAMILY', 'REP_START', 'REP_END', 'REP_LEFT'
                    ]
                ))

        # Merge
        df = pd.concat(record_list, axis=1).T
        df.set_index('ID', inplace=True, drop=False)

        # Fix fields
        #df['BETTER_OVERLAP_HIT'] = df['BETTER_OVERLAP_HIT'].apply(lambda val: val == '*')

        df['ALIGN_SCORE'] = df['ALIGN_SCORE'].astype(np.int32)
        df['DIFF_PROP'] = df['DIFF_PROP'].astype(np.float32)
        df['DEL_PROP'] = df['DEL_PROP'].astype(np.float32)
        df['INS_PROP'] = df['INS_PROP'].astype(np.float32)
        df['SV_MATCH_START'] = df['SV_MATCH_START'].astype(np.int32)
        df['SV_MATCH_END'] = df['SV_MATCH_END'].astype(np.int32)
        df['SV_MATCH_LEFT'] = df['SV_MATCH_LEFT'].astype(np.int32)
        df['REP_START'] = df['REP_START'].astype(np.int32)
        df['REP_END'] = df['REP_END'].astype(np.int32)
        df['REP_LEFT'] = df['REP_LEFT'].astype(np.int32)

        # Add proportions
        df['SVLEN'] = pd.read_csv(input.bed, sep='\t', header=0, usecols=('ID', 'SVLEN'), index_col='ID', squeeze=True)

        df['REP_PROP'] = (df['REP_START'] - 1 + df['REP_LEFT']) / (df['REP_END'] + df['REP_LEFT'])
        df['SV_PROP'] = (df['SV_MATCH_START'] - 1 + df['SV_MATCH_LEFT']) / df['SVLEN']

        del(df['SVLEN'])

        # Rearrange
        head_cols = ['ID', 'REPEAT_FAMILY', 'REPEAT_TYPE', 'REP_PROP', 'SV_PROP']

        df = df.loc[:, head_cols + [col for col in df.columns if col not in head_cols]]

        # Write
        df.to_csv(output.tsv, sep='\t', index=False)

# variant_anno_repeat_rmsk_run
#
# Run RepeatMasker.
rule variant_anno_repeat_rmsk_run:
    input:
        fasta='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz'
    output:
        out='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/rmsk/data/rmsk-table_{vartype}_{svtype}/rmsk.out.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup|rgn|sub'
    shadow:
        'shallow'
    params:
        threads=4
    run:
        if os.stat(input.fasta).st_size > 0:
            shell(
                """RMSK_TMP=$(dirname {output.out}); """
                """zcat {input.fasta} > sv.fa; """
                """RepeatMasker """
                    """-species "Homo sapiens" """
                    """-dir ./ """
                    """-xsmall """
                    """-no_is """
#                    """-e wublast """
                    """-s """
                    """-pa {params.threads} """
                    """sv.fa; """
                """ls -al; """
                """echo "####"; """
                """gzip -c sv.fa.out > {output.out}"""
            )

        else:
            shell("""touch {output.out}""")


#
# TRF
#

# variant_anno_repeat_trf_table
#
# Make table of repeat annotations.
rule variant_anno_repeat_trf_table:
    input:
        out='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/trf/data/trf-table_{vartype}_{svtype}/trf.out'
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/trf/trf-table_{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup'
    run:

        # Create empty tab file if trf output is empty
        if os.stat(input.out).st_size == 0:
            df = pd.DataFrame(
                [],
                columns=[
                    'ID', 'SV_POS', 'SV_END', 'PERIOD', 'COPIES', 'CONSENSUS_SIZE',
                    'MATCH_PCT', 'INDEL_PCT', 'ALIGN_SCORE',
                    'BASE_A_PCT', 'BASE_C_PCT', 'BASE_G_PCT', 'BASE_T_PCT',
                    'ENTROPY', 'CONSENSUS',
                    'LONGEST_PERFECT_N', 'LONGEST_PERFECT_BP'
                ]
            )

            df.to_csv(output.tab, sep='\t', index=False)

            return

        # Summarize records
        record_list = list()

        sv_id = None

        line_count = 0

        with open(input.out, 'r') as in_file:
            for line in in_file:
                line_count += 1

                line = line.strip()

                if not line:
                    continue

                # Get SV ID
                if line.startswith('@'):
                    sv_id = line[1:]
                    continue

                elif sv_id is None:
                    raise RuntimeError('First line in TRF output file must begin with "@"')

                # Parse
                tok = re.split('\s+', line)

                if len(tok) != 17:
                    raise RuntimeError('Expected 17 fields on line {}: Received {}'.format(line_count, len(tok)))

                # BED coordinates
                tok[0] = int(tok[0]) - 1

                # Add record
                record_list.append(pd.Series(
                    [sv_id] + tok,
                    index=[
                        'ID',
                        'SV_POS', 'SV_END',
                        'PERIOD', 'COPIES', 'CONSENSUS_SIZE',
                        'MATCH_PCT', 'INDEL_PCT', 'ALIGN_SCORE',
                        'BASE_A_PCT', 'BASE_C_PCT', 'BASE_G_PCT', 'BASE_T_PCT', 'ENTROPY',
                        'CONSENSUS', 'SEQ',
                        'FLANK_L', 'FLANK_R'
                    ]
                ))

        # Merge dataframe
        df = pd.concat(record_list, axis=1).T

        df['PERIOD'] = df['PERIOD'].astype(np.int64)

        # Define: Function to find longest perfect match. Start with lowest position index of the consensus repeat in
        # the string, then find all other sequential positions.
        def longest_match(row):
            match_set = {m.start() for m in re.finditer(row['CONSENSUS'], row['SEQ'])}

            max_count = 0

            while match_set:
                match_count = 1
                next_index = min(match_set)

                match_set = match_set - {next_index}

                next_index += row['PERIOD']

                while next_index in match_set:
                    match_count += 1
                    match_set = match_set - {next_index}
                    next_index += row['PERIOD']

                if match_count > max_count:
                    max_count = match_count

            return max_count

        # Find longest perfect match
        df['LONGEST_PERFECT_N'] = df.apply(longest_match, axis=1)
        df['LONGEST_PERFECT_BP'] = df['LONGEST_PERFECT_N'] * df['PERIOD']

        # Trim sequence and flanking sequence
        del(df['SEQ'])
        del(df['FLANK_L'])
        del(df['FLANK_R'])

        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

# variant_anno_repeat_trf_run
#
# Run TRF annotations.
#
# TRF Parameters:
# * Match = 2
# * Mismatch = 7
# * Delta (indel penalty) = 7
# * PM (prob. match) = 80
# * PI (prop. indel) = 10
# * Minscore = 20
# * MaxPeriod = 2000
rule variant_anno_repeat_trf_run:
    input:
        fa='temp/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa'
    output:
        out='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/trf/data/trf-table_{vartype}_{svtype}/trf.out'
    wildcard_constraints:
        svtype='ins|del|inv|dup|rgn|sub'
    run:

        if os.stat(input.fa).st_size > 0:
            shell(
                """INPUT_FA=$(readlink -f {input.fa}); """
                """cd $(dirname {output.out}); """
                """{SVPOP_DIR}/scripts/anno/trf ${{INPUT_FA}} 2 7 7 80 10 20 2000 -ngs -h > $(basename {output.out})"""
            )
        else:
            shell("""touch {output.out}""")
