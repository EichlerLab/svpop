"""
Annotations based on sequence content.
"""

###################
### Definitions ###
###################


#############
### Rules ###
#############


#
# GC content
#

# variant_anno_caller_seqcontent_gc
#
# Annotate GC content.
rule variant_anno_seqcontent_gc:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        fa='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz'
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/gc/gc_content_{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup|rgn|sub'
    run:

        # Read variants
        df = pd.read_csv(input.bed, sep='\t', header=0, usecols=('ID', 'SVTYPE', 'SVLEN'), index_col='ID')

        # Get GC
        if os.stat(input.fa).st_size > 0:

            gc_list = list()

            # Read table
            with gzip.open(input.fa, 'rt') as in_file:
                for record in Bio.SeqIO.parse(in_file, 'fasta'):
                    seq = str(record.seq).upper()
                    seq = re.sub('[^ACGT]', '', seq)

                    if len(seq) > 0:
                        gc = len(re.sub('A|T', '', seq.upper())) / len(seq)
                    else:
                        gc = np.nan

                    gc_list.append(pd.Series(
                        [record.id, gc],
                        index=['ID', 'GC_CONTENT']
                    ))

            # Assign GC to df
            df['GC_CONTENT'] = pd.concat(gc_list, axis=1).T.set_index('ID').squeeze()

        else:
            # Empty FASTA, assign NA to all
            df['GC_CONTENT'] = np.nan

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, header=True, na_rep='NA', compression='gzip')


# variant_anno_caller_seqcontent_seq
#
# Get sequence from FASTA and pull into a table.
rule variant_anno_seqcontent_seq:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        fa='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz'
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/seq/seq_{vartype}_{svtype}.tsv.gz'
    run:

        # Read variants
        id_set = set(pd.read_csv(input.bed, sep='\t', header=0, usecols=('ID', ))['ID'])

        # Get SEQ
        if os.stat(input.fa).st_size > 0:

            seq_list = list()

            # Read table
            with gzip.open(input.fa, 'rt') as in_file:
                for record in Bio.SeqIO.parse(in_file, 'fasta'):
                    seq = str(record.seq).upper()

                    if len(seq) == 0:
                        seq = np.nan

                    seq_list.append(pd.Series(
                        [record.id, seq],
                        index=['ID', 'SEQ']
                    ))

            # Assign GC to df
            df_seq = pd.concat(seq_list, axis=1).T

            df_seq = df_seq.loc[df_seq['ID'].apply(lambda val: val in id_set)]

        else:
            # Empty FASTA, assign NA to all
            df_seq = pd.DataFrame([], columns=['ID', 'SEQ'])

        # Write
        df_seq.to_csv(output.tsv, sep='\t', index=False, header=True, compression='gzip')

# variant_anno_seqcontent_hom_ref
#
# Breakpoint homology vs reference for sequence-resolved SVs
rule variant_anno_seqcontent_hom_ref:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        fa='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz'
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/seq/break_hom_ref_{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
        svtype='ins|del'
    run:

        chrom_name = None
        chrom_seq = None

        # Read variants
        df_sv = pd.read_csv(input.bed, sep='\t', header=0, usecols=('#CHROM', 'POS', 'END', 'ID', 'SVTYPE'), index_col='ID')
        df_sv = df_sv[['#CHROM', 'POS', 'END', 'SVTYPE']]

        # Get SEQ
        if os.stat(input.fa).st_size > 0:

            seq_list = list()

            # Read table
            with gzip.open(input.fa, 'rt') as in_file, pysam.FastaFile(config['reference']) as ref_fa:
                for record in Bio.SeqIO.parse(in_file, 'fasta'):
                    seq = str(record.seq).upper()

                    if len(seq) == 0:
                        continue

                    chrom, pos, end, svtype = df_sv.loc[record.id]

                    # Get reference sequence
                    if chrom != chrom_name:
                        chrom_seq = ref_fa.fetch(chrom).upper()
                        chrom_name = chrom

                    # Get homology
                    seq_list.append(pd.Series(
                        [
                            record.id,
                            svpoplib.variant.left_homology(pos - 1, chrom_seq, seq),
                            svpoplib.variant.right_homology(end - (1 if svtype == 'INS' else 0), chrom_seq, seq)
                        ],
                        index=['ID', 'HOM_REF_L', 'HOM_REF_R']
                    ))

            # Merge
            if len(seq_list) > 0:
                df_hom = pd.concat(seq_list, axis=1).T
            else:
                df_hom = None  # Signal empty dataframe

        else:
            # Empty FASTA, assign NA to all
            df_hom = None  # Signal empty dataframe

        # Empty dataframe
        if df_hom is None:
            df_hom = pd.DataFrame([], columns=['ID', 'HOM_REF_L', 'HOM_REF_R'])

        # Write
        df_hom.to_csv(output.tsv, sep='\t', index=False, header=True, compression='gzip')