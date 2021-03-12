"""
Reference related tasks
"""

# data_ref_contig_table
#
# Contig table.
rule data_ref_contig_table:
    output:
        tsv='data/ref/contig_info.tsv.gz'
    run:

        svpoplib.ref.get_ref_info(
            config['reference']
        ).to_csv(
            output.tsv, sep='\t', index=True, compression='gzip'
        )
