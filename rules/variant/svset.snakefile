"""
Filter variants based on variant characteristics.
"""

# variant_svsetfilter_bed
#
# Filter by SV set.
rule variant_svsetfilter_bed:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        filter_files=lambda wildcards: svpoplib.svset.get_filter_input_files(wildcards.svset, wildcards, config)
    output:
        bed='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz'
    wildcard_constraints:
        svset='((?!all).*|all.+)',  # Do not allow "all"
        svtype='ins|del|inv|dup|sub|rgn|snv|insdel|insdelinv'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t')

        # Filter
        df = svpoplib.svset.apply_svset_filter(
            pd.read_csv(input.bed, sep='\t'),
            wildcards.svset,
            wildcards,
            config
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# variant_svsetfilter_anno
#
# Filter variant annotations by svset.
rule variant_svsetfilter_anno:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz',
        anno='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/anno/{annodir}/{annotype}_{vartype}_{svtype}.{ext}.gz'
    output:
        anno='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/{annodir}/{annotype}_{vartype}_{svtype}.{ext}.gz'
    wildcard_constraints:
        svset='((?!all).*|all.+)',  # Do not allow "all"
        svtype='ins|del|inv|dup|sub|rgn|snv|insdel|insdelinv',
        ext='tsv|bed'
    run:

        id_set = set(pd.read_csv(input.bed, sep='\t', usecols=('ID',))['ID'])

        df_iter = pd.read_csv(input.anno, sep='\t', chunksize=5000, iterator=True)

        write_header = True

        with gzip.open(output.anno, 'wt') as out_file:
            for df in df_iter:

                df = df.loc[df['ID'].isin(id_set)]

                if df.shape[0] > 0 or write_header:
                    df.to_csv(out_file, sep='\t', index=False, header=write_header)
                    write_header = False


# variant_svsetfilter_fa
#
# Filter variant sequence FASTA by svset.
rule variant_svsetfilter_fa:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz',
        fa='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz'
    output:
        fa='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/fa/{vartype}_{svtype}.fa.gz'
    wildcard_constraints:
        svset='((?!all).*|all.+)',  # Do not allow "all"
        svtype='ins|del|inv|dup|sub|rgn|snv|insdel|insdelinv'
    run:

        id_set = set(pd.read_csv(input.bed, sep='\t', usecols=('ID',))['ID'])

        # Filter
        with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
            SeqIO.write(svpoplib.seq.fa_to_record_iter(input.fa, id_set), out_file, 'fasta')
