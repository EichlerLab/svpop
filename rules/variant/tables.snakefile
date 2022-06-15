"""
Data summary tables.
"""

# variant_tables_caller_xlsx
#
# TSV to Excel.
rule variant_tables_caller_xlsx:
    input:
        tab='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/tables/variant_summary/{vartype}_{svtype}.tsv.gz'
    output:
        xlsx='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/tables/variant_summary/{vartype}_{svtype}.xlsx'
    run:

        pd.read_csv(
            input.tab, sep='\t'
        ).to_excel(
            output.xlsx, index=False
        )

# variant_tables_caller_summary_indel
#
# Get a table line summarizing variant calls.
rule variant_tables_caller_summary_insdel:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz'
    output:
        tsv='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/tables/variant_summary/{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
        vartype='sv|indel|sub|rgn'
    run:

        # Read summary data
        df = pd.read_csv(input.bed, sep='\t', header=0)

        if df.shape[0] > 0:

            # Summarize
            df_summary = df.groupby('SVTYPE')['SVLEN'].apply(
                lambda vals:
                    pd.DataFrame(
                        pd.Series(
                            [wildcards.sourcetype, wildcards.sourcename, wildcards.sample, len(vals), np.mean(vals), np.sum(vals)],
                            index=['SOURCETYPE', 'SOURCENAME', 'SAMPLE', 'N', 'MEAN', 'BP']
                        )
                    ).T.astype(
                        {'SOURCETYPE': str, 'SOURCENAME': str, 'SAMPLE': str, 'N': np.int, 'MEAN': np.float, 'BP': np.int}
                    )
            )

            if df.shape[0] > 1:
                df_summary.index = df_summary.index.droplevel(1)  # Removes index level artifact from Series
                df_summary.reset_index(inplace=True)
            else:
                df_summary.insert(0, 'SVTYPE', df.loc[0, 'SVTYPE'])

            df_summary = df_summary.loc[:, ['SOURCETYPE', 'SOURCENAME', 'SAMPLE', 'SVTYPE'] + list(df_summary.columns)[4:]]

            # Sort
            df_summary['SVTYPE'] = pd.Categorical(df_summary['SVTYPE'], ['INS', 'DEL', 'INV', 'DUP'])
            df_summary.sort_values('SVTYPE', inplace=True)

            # Write
            df_summary.to_csv(output.tsv, sep='\t', index=False)

        else:
            pd.DataFrame(columns=['SOURCETYPE', 'SOURCENAME', 'SAMPLE', 'SVTYPE', 'N', 'MEAN', 'BP']).to_csv(output.tsv, sep='\t', index=False)


# variant_tables_caller_summary_snv
#
# Get a table line summarizing SNV variant calls.
rule variant_tables_caller_summary_snv:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/snv_snv.bed.gz'
    output:
        tsv='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/tables/variant_summary/snv_snv.tsv.gz'
    run:

        # Read summary data
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Normalize REF/ALT
        df['REF'].fillna('N', inplace=True)
        df['ALT'].fillna('N', inplace=True)

        df['REF'] = df['REF'].apply(lambda val: val.upper())
        df['ALT'] = df['ALT'].apply(lambda val: val.upper())

        ref_non_n = df['REF'].apply(lambda val: val in {'A', 'C', 'G', 'T'})
        alt_non_n = df['ALT'].apply(lambda val: val in {'A', 'C', 'G', 'T'})

        n_snv_with_non_acgt = df.shape[0]

        # Filter non-ACGT in REF or ALT
        df = df.loc[ref_non_n & alt_non_n]

        n_snv = df.shape[0]

        # Compute TI/TV
        df_refalt_count = df.groupby(['REF', 'ALT'])['ID'].count()

        ti_sum = np.sum(df_refalt_count.loc[[
            ('A', 'G'), ('G', 'A'),
            ('C', 'T'), ('T', 'C')
        ]])

        tv_sum = n_snv - ti_sum

        if df.shape[0] > 0:
            df_stats = pd.Series(
                [
                    wildcards.sourcetype,
                    wildcards.sourcename,
                    wildcards.sample,
                    'SNV',

                    n_snv,
                    (ti_sum / tv_sum) if tv_sum > 0 else np.nan,

                    df_refalt_count[('A', 'C')] / n_snv,
                    df_refalt_count[('A', 'G')] / n_snv,
                    df_refalt_count[('A', 'T')] / n_snv,

                    df_refalt_count[('C', 'A')] / n_snv,
                    df_refalt_count[('C', 'G')] / n_snv,
                    df_refalt_count[('C', 'T')] / n_snv,

                    df_refalt_count[('G', 'A')] / n_snv,
                    df_refalt_count[('G', 'C')] / n_snv,
                    df_refalt_count[('G', 'T')] / n_snv,

                    df_refalt_count[('T', 'A')] / n_snv,
                    df_refalt_count[('T', 'C')] / n_snv,
                    df_refalt_count[('T', 'G')] / n_snv,

                    n_snv_with_non_acgt - n_snv
                ],
                index=[
                    'SOURCETYPE', 'SOURCENAME', 'SAMPLE', 'SVTYPE',
                    'N', 'TITV_RATIO',
                    'A_C', 'A_G', 'A_T',
                    'C_A', 'C_G', 'C_T',
                    'G_A', 'G_C', 'G_T',
                    'T_A', 'T_C', 'T_G',
                    'NON_ACGT'
                ]
            )

            df_stats = pd.DataFrame(df_stats).T

            df_stats['N'] = df_stats['N'].astype(np.int32)
            df_stats['NON_ACGT'] = df_stats['NON_ACGT'].astype(np.int32)

            df_stats.to_csv(output.tsv, sep='\t', index=False)


        else:
            pd.DataFrame(
                columns=[
                    'N', 'TITV_RATIO',
                    'A_C', 'A_G', 'A_T',
                    'C_A', 'C_G', 'C_T',
                    'G_A', 'G_C', 'G_T',
                    'T_A', 'T_C', 'T_G',
                    'NON_ACGT'
                ]
            ).to_csv(output.tsv, sep='\t', index=False)
