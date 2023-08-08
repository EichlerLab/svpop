"""
Intersect SVs by region for multiple annotations.
"""

###################
### Definitions ###
###################



#############
### Rules ###
#############

#
# Regions - Multiple annotations
#

# variant_anno_caller_region_intersect
#
# Intersect annotations by region.
rule variant_anno_caller_region_intersect:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        anno_bed='data/anno/{annotype}/{annoname}_regions_{distance}_{flank}.bed.gz'
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/{annotype}/{annoname}_regions_{distance}_{flank}_{overlap}_{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
        vartype='sv|indel|snv|rgn|sub',
        svtype='ins|del|inv|dup|snv|rgn|sub'
    run:

        # Get overlap parameters
        if wildcards.overlap != 'any':
            overlap_proportion = int(wildcards.overlap) / 100

            if overlap_proportion <= 0 or overlap_proportion > 1:
                raise ValueError(
                    'Overlap length must be between 0 (exclusive) and 100 (inclusive): {}'.format(wildcards.overlap)
                )

            overlap_params = '-f {}'.format(overlap_proportion)

        else:
            overlap_params = ''

        # Do intersect
        shell(
            """{{\n"""
            """    echo "ID";\n"""
            """    bedtools intersect -a <(zcat {input.bed} | cut -f1-4) -b {input.anno_bed} -wa {overlap_params} -u | cut -f 4;\n"""
            """}} |"""
            """gzip > {output.tsv}"""
        )

# variant_anno_caller_region_intersect
#
# Intersect annotations by region.
rule variant_anno_caller_region_intersect_basecount:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        anno_bed='data/anno/{annotype}/{annoname}_regions_{distance}_{flank}.bed.gz'
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/{annotype}/{annoname}_basecount_{distance}_{flank}_{overlap}_{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    run:

        # Get overlap parameters
        if wildcards.overlap != 'any':
            overlap_proportion = int(wildcards.overlap) / 100

            if overlap_proportion <= 0 or overlap_proportion > 1:
                raise ValueError(
                    'Overlap length must be between 0 (exclusive) and 100 (inclusive): {}'.format(wildcards.overlap)
                )

            overlap_params = '-f {}'.format(overlap_proportion)

        else:
            overlap_params = ''

        # Do intersect
        shell(
            """{{\n"""
            """    echo -e "ID\tBP";\n"""
            """    bedtools intersect -a <(zcat {input.bed} | cut -f1-4) -b <(zcat {input.anno_bed} | cut -f1-3) -wao {overlap_params} | \n"""
            """    cut -f4,8 | \n"""
            """    sort | \n"""
            """    awk '{{arr[$1]+=$2}} END {{for (key in arr) printf("%s\\t%s\\n", key, arr[key])}}'; \n"""
            """}} |"""
            """gzip > {output.tsv}"""
        )

#
# Chromosome Band
#

# variant_anno_caller_region_band
#
# Get band for each variant. Excludes SVs on unplaced and unlocalized chromosomes.
rule variant_anno_caller_region_band:
    input:
        bed='temp/variant/caller/{sourcename}/{sample}/{filter}/all/anno/bands/bands_{vartype}_{svtype}.bed.gz'
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/bands/bands_{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)

        # Write empty table
        if df.shape[0] == 0:
            pd.DataFrame([], columns=['ID', 'BAND']).to_csv(output.tab, sep='\t', index=False)
            return

        # Get band prefix (chr name without "chr"), np.nan for unplaced/unlocalized contigs
        df['BAND_PREFIX'] = df['#CHROM'].apply(lambda val: re.match('chr([0-9]*[0-9XY])$', val))
        df['BAND_PREFIX'] = df['BAND_PREFIX'].apply(lambda item: item.group(1) if item else np.nan)
        df.dropna(inplace=True)

        df['BAND'] = df.apply(lambda row: '{BAND_PREFIX}{name}'.format(**row), axis=1)

        # Group by ID (incase variant affects more than one band), join with commas, and write
        df.groupby('ID')['BAND'].apply(lambda vals: ','.join(vals)).to_csv(
            output.tsv, sep='\t', index=True, header=True, compression='gzip'
        )

# variant_anno_caller_region_band_intersect
#
# Intersect SVs with chromosome bands.
rule variant_anno_caller_region_band_intersect:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        bands='data/anno/bands/bands.bed.gz'
    output:
        bed=temp('temp/variant/caller/{sourcename}/{sample}/{filter}/all/anno/bands/bands_{vartype}_{svtype}.bed.gz')
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    shell:
        """{{\n"""
        """    echo -e $(cut -f 1-4 <(zcat {input.bed}) | head -n 1) "\\t" $(head -n 1 <(zcat {input.bands})) | sed -re 's/\s+/\t/g';\n"""
        """    cut -f 1-4 <(zcat {input.bed}) | bedtools intersect -a stdin -b {input.bands} -sorted -loj;\n"""
        """}} | gzip > {output.bed}"""

#
# CpG sites
#

# variant_anno_caller_region_cpg_sites
#
# Find CpG sites in the reference.
rule variant_anno_caller_region_cpg_sites:
    input:
        ref=config['reference']
    output:
        bed='data/anno/cpgsites/cpgsites.bed.gz'
    run:

        # Crawl bases for CpG sites
        record_list = list()

        with open(input.ref, 'r') as in_file:
            for record in SeqIO.parse(in_file, 'fasta'):

                # Init sequence
                seq = str(record.seq).upper()
                chrom = record.id

                print(chrom)

                # Find positions of CpG's
                fwd_pos = [m.start() for m in re.finditer('CG', seq)]  # Thanks https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring
                rev_pos = [m.start() for m in re.finditer('GC', seq)]

                for pos in fwd_pos:
                    record_list.append(pd.Series(
                        [
                            chrom,
                            pos,
                            pos + 2,
                            '+'
                        ],
                        index=['#CHROM', 'POS', 'END', 'STRAND']
                    ))

                for pos in rev_pos:
                    record_list.append(pd.Series(
                        [
                            chrom,
                            pos,
                            pos + 2,
                            '-'
                        ],
                        index=['#CHROM', 'POS', 'END', 'STRAND']
                    ))

        # Merge and sort
        df = pd.concat(record_list, axis=0).sort_values(['#CHROM', 'POS'])

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

#
# SD
#

# variant_anno_caller_sd_max
#
# Intersect with SD annotations and determine the maximum identity
rule variant_anno_caller_sd_max:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        sd_bed='data/anno/sd/sd-max-{match_type}.bed.gz'
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/sd/sd-max-{match_type}_{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    shell:
        """{{\n"""
        """    echo -e "ID\\tMATCH_MIN\\tMATCH_MAX";\n"""
        """    cut -f1-4 <(zcat {input.bed}) |\n"""
        """    bedtools map -a stdin -b {input.sd_bed} -c 4 -o min,max |\n"""
        """    cut -f4-6\n"""
        """}} | gzip > {output.tsv}"""


#
# RMSK (RepeatMasker)
#

rule variant_anno_rmsk:
    input:
        bed='results/variant/caller/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz',
        tsv='temp/variant/caller/{sourcename}/{sample}/{filter}/all/anno/rmsk/rmsk-{filter_spec}-{ident}_intersect_{vartype}_{svtype}.tsv.gz'
    output:
        tsv='results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/rmsk/rmsk-{filter_spec}-{ident}_intersect_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        df = pd.read_csv(input.tsv, sep='\t')

        # Remove BED columns from RMSK
        del(df['#CHROM.1'])
        del(df['POS.1'])
        del(df['END.1'])

        df.columns = [col if col != 'ID.1' else 'ID_RMSK' for col in df.columns]

        # Get coordinates relative to the original range
        df_coord = pd.read_csv(input.bed, sep='\t', index_col='ID', usecols=['ID', 'POS', 'END'])

        df['VAR_POS'] = df.apply(lambda row: row['POS'] - df_coord.loc[row['ID'], 'POS'], axis=1)
        df['VAR_END'] = df['VAR_POS'] + (df['END'] - df['POS'])

        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

# variant_anno_rmsk
#
# Intersect (LOJ) with RMSK.
rule variant_anno_rmsk_intersect:
    input:
        bed='temp/variant/caller/{sourcename}/{sample}/{filter}/all/bed4/{vartype}_{svtype}.bed.gz',
        rmsk_bed='data/anno/rmsk/rmsk-{filter_spec}-{ident}.bed.gz'
    output:
        tsv=temp('temp/variant/caller/{sourcename}/{sample}/{filter}/all/anno/rmsk/rmsk-{filter_spec}-{ident}_intersect_{vartype}_{svtype}.tsv.gz')
    shell:
        """{{\n"""
        """    {{\n"""
        """        head -n 1 <(zcat {input.bed});\n"""
        """        head -n 1 <(zcat {input.rmsk_bed});\n"""
        """    }} | \n"""
        """    awk -vORS="\\t" {{print}} | \n"""
        """    sed -re 's/\\s+/\\t/g' | \n"""
        """    sed -re 's/\\s*$/\\n/';\n"""
        """    bedtools intersect -a {input.bed} -b {input.rmsk_bed} -wb;\n"""
        """}} | gzip > {output.tsv}"""
