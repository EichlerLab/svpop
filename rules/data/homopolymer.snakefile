"""
Make homopolymer reference annotations. Makes a BED file with homopolymer locations (2+ bp) including
#CHROM, POS, END, BASE, and COUNT (homopolymer length).
"""

# hpref_merge_bed
#
# Merge chrom BEDs.
rule hpref_merge_bed:
    input:
        bed=expand(
            'temp/data/anno/homopolymer/{{rep_len}}_regions_{chrom}.bed.gz',
            chrom=sorted(list(svpoplib.ref.get_df_fai(config['reference_fai']).index))
        )
    output:
        bed='data/anno/homopolymer/{rep_len}_regions_all.bed'
    wildcard_constraints:
        rep_len='homopolymer|dinucleotide'
    run:

        with open(output.bed, 'wt') as out_file:
        #with gzip.open(output.bed, 'wt') as out_file:
            out_file.write('#CHROM\tPOS\tEND\tBASE\tCOUNT\n')

            for bed_file_name in input.bed:
                with gzip.open(bed_file_name, 'rt') as in_file:
                    for line in in_file:
                        out_file.write(line)

# hpref_bed_chrom
#
# Make homopolymer BED for one chrom
rule hpref_bed_chrom:
    output:
        bed=temp('temp/data/anno/homopolymer/homopolymer_regions_{chrom}.bed.gz')
    run:

        chrom = wildcards.chrom

        with gzip.open(output.bed, 'wt') as out_file:
            with pysam.FastaFile(config['reference']) as in_file:

                # Get sequence
                seq = in_file.fetch(chrom).upper()
                seq_len = len(seq)

                # Init
                last_base = 'N'
                pos = 0
                count = 0

                # Crawl sequences
                for base in seq:

                    if base == last_base:
                        count += 1

                    else:

                        if count > 1 and last_base in {'A', 'C', 'G', 'T'}:
                            out_file.write('{:}\t{:d}\t{:d}\t{}\t{:d}\n'.format(chrom, pos - count, pos, last_base, count))

                        count = 1
                        last_base = base

                    # Next position
                    pos += 1

                # Write last position (if homopolymer)
                if count > 1 and last_base in {'A', 'C', 'G', 'T'}:
                    out_file.write('{:}\t{:d}\t{:d}\t{}\t{:d}\n'.format(chrom, pos - count, pos, last_base, count))

# hpref_bed_chrom
#
# Make dinucleotide BED for one chrom
rule hpref_bed_chrom_dinucl:
    output:
        bed=temp('temp/data/anno/homopolymer/dinucleotide_regions_{chrom}.bed.gz')
    run:

        chrom = wildcards.chrom

        with gzip.open(output.bed, 'wt') as out_file:
            with pysam.FastaFile(config['reference']) as in_file:

                # Get sequence
                seq = in_file.fetch(chrom).upper()
                seq_len = len(seq)

                # Init
                i = 0

                while i < len(seq) - 2:

                    # Skip if dinucleotide repeat does not start at pos i
                    if seq[i + 2] != seq[i] or seq[i + 1] == seq[i]:
                        i += 1
                        continue

                    # Get bases
                    a = seq[i]
                    b = seq[i + 1]

                    j = i + 2

                    while j < len(seq) - 1 and seq[j] == a and seq[j + 1] == b:
                        j += 2

                    rep_len = (j - i) // 2

                    if rep_len > 1:
                        out_file.write('{:}\t{:d}\t{:d}\t{}\t{:d}\n'.format(chrom, i, j, a + b, rep_len))

                    i = j - 1

