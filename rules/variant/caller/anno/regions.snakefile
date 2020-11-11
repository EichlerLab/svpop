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
        bed='results/variant/{sourcetype}/{caller}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        anno_bed='data/anno/{annotype}/{annoname}_regions_{distance}_{flank}.bed'
    output:
        tab='results/variant/{sourcetype}/{caller}/anno/{sample}/all/{filter}/{annotype}/{annoname}_regions_{distance}_{flank}_{overlap}_{vartype}_{svtype}.tab'
    wildcard_constraints:
        sourcetype='caller|varset',
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
            """echo "ID" > {output.tab}; """
            """bedtools intersect -a {input.bed} -b {input.anno_bed} -wa {overlap_params} -u | """
            """cut -f4 """
            """>> {output.tab}"""
        )


#
# Chromosome Band
#

# variant_anno_caller_region_band
#
# Get band for each variant. Excludes SVs on unplaced and unlocalized chromosomes.
rule variant_anno_caller_region_band:
    input:
        bed='temp/variant/{sourcetype}/{caller}/anno/{sample}/all/{filter}/bands/bands_{vartype}_{svtype}.bed'
    output:
        tab='results/variant/{sourcetype}/{caller}/anno/{sample}/all/{filter}/bands/bands_{vartype}_{svtype}.tab'
    wildcard_constraints:
        sourcetype='caller|varset',
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
            output.tab, sep='\t', index=True, header=True
        )

# variant_anno_caller_region_band_intersect
#
# Intersect SVs with chromosome bands.
rule variant_anno_caller_region_band_intersect:
    input:
        bed='results/variant/{sourcetype}/{caller}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        bands='data/anno/bands/bands.bed'
    output:
        bed=temp('temp/variant/{sourcetype}/{caller}/anno/{sample}/all/{filter}/bands/bands_{vartype}_{svtype}.bed')
    wildcard_constraints:
        sourcetype='caller|varset',
        svtype='ins|del|inv|dup|snv|rgn|sub'
    shell:
        """echo -e $(cut -f 1-4 {input.bed} | head -n 1) "\\t" $(head -n 1 {input.bands}) | """
        """sed -re 's/\s+/\t/g' """
        """>{output.bed};"""
        """cut -f 1-4 {input.bed} | """
        """bedtools intersect -a stdin -b {input.bands} -sorted -loj """
        """>>{output.bed}"""

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
        bed='data/anno/cpgsites/cpgsites.bed'
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
        df.to_csv(output.bed, sep='\t', index=False)

#
# SD
#

# variant_anno_caller_sd_max
#
# Intersect with SD annotations and determine the maximum identity
rule variant_anno_caller_sd_max:
    input:
        bed='results/variant/caller/{caller}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        sd_bed='data/anno/sd/sd-max-{match_type}.bed'
    output:
        tab='results/variant/caller/{caller}/anno/{sample}/all/{filter}/sd/sd-max-{match_type}_{vartype}_{svtype}.tab'
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    shell:
        """{{\n"""
        """    echo -e "ID\\tMATCH_MIN\\tMATCH_MAX";\n"""
        """    cut -f1-4 {input.bed} |\n"""
        """    bedtools map -a stdin -b {input.sd_bed} -c 4 -o min,max |\n"""
        """    cut -f4-6\n"""
        """}} > {output.tab}"""
