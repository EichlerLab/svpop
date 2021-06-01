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
rule variant_anno_caller_seqcontent_gc:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        fa='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz'
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/gc/gc_content_{vartype}_{svtype}.tsv.gz'
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
rule variant_anno_caller_seqcontent_seq:
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
