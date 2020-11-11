"""
Import external BED with variant calls for samples.
"""


###################
### Definitions ###
###################

def _variant_caller_extern_get_fa_input(wildcards):

    # Get FASTA file
    fasta_file = analib.variant_extern.get_fa(wildcards, config)

    # No input if no FASTA
    if fasta_file is None:
        return []

    # Return path to FASTA
    return [fasta_file]


#############
### Rules ###
#############

# variant_caller_extern_fa
#
# Get sequence FASTA.
rule variant_caller_extern_fa:
    input:
        bed='results/variant/caller/extern-{extern_source}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        fa=_variant_caller_extern_get_fa_input,
        ref=config['reference']
    output:
        fa='results/variant/caller/extern-{extern_source}/fasta/{sample}/all/{filter}/{vartype}_{svtype}.fa.gz'
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|rgn|snv'
    run:

        # Get FASTA file
        fasta_file = analib.variant_extern.get_fa(wildcards, config)
        bed_file = analib.variant_extern.get_bed(wildcards, config)

        # Re-write ensuring FASTA is bgzipped
        if fasta_file is not None and bed_file is not None:

            if os.stat(fasta_file).st_size > 0:
                with analib.seq.PlainOrGzReader(fasta_file) as in_file:
                    with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
                        SeqIO.write(SeqIO.parse(in_file, 'fasta'), out_file, 'fasta')

            else:
                # Write an empty file if source is empty
                with open(output.fa, 'w') as in_file:
                    pass

            return

        # FASTA was not given, create if it can be constructed from the reference.
        if wildcards.svtype in {'del', 'inv', 'dup', 'rgn'} and os.stat(input.bed).st_size > 0:

            # Read
            df = pd.read_csv(input.bed, sep='\t', header=0)

            # Create FASTA from the reference if there are SV BED records
            if df.shape[0] > 0:

                # Get sequences
                with pysam.FastaFile(input.ref) as in_file:
                    df['SEQ'] = df.apply(lambda row: in_file.fetch(row['#CHROM'], row['POS'], row['END']), axis=1)

                # Write FASTA
                with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
                    SeqIO.write(analib.seq.bed_to_seqrecord_iter(df), out_file, 'fasta')

                shell("""samtools faidx {output.fa}""")

                # Done writing
                return


        # Fall through to this point means no pre-defined FASTA and:
        #  a) SV sequence cannot be pulled from the reference
        #    - or -
        #  b) No SV calls in the BED file

        shell("""touch {output.fa}""")


# variant_caller_extern_filter_bed
#
# Apply region filter to variant BED.
rule variant_caller_extern_filter_bed:
    input:
        bed='temp/variant/caller/extern-{extern_source}/bed/{sample}/pre_filter/{vartype}_{svtype}.bed',
        filter=lambda wildcards: analib.variant.get_filter_bed(wildcards.filter, UCSC_REF_NAME, config, SVPOP_DIR)
        #filter=os.path.join(SVPOP_DIR, 'files/filter/{}/{{filter}}.bed'.format(UCSC_REF_NAME))
    output:
        bed='results/variant/caller/extern-{extern_source}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        bed_filt='results/variant/caller/extern-{extern_source}/bed/{sample}/all/{filter}/byref/removed/{vartype}_{svtype}.bed'
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|rgn|snv',
        filter='((?!all)[^\/]*)|all.+',
    run:

        if wildcards.filter != 'all':
            shell(
                """bedtools intersect -wa -v -sorted -a {input.bed} -b {input.filter} -header > {output.bed}; """
                """bedtools intersect -wa -u -sorted -a {input.bed} -b {input.filter} -header > {output.bed_filt}; """  # Capture filtered variants
            )
        else:
            shell(
                """cp {input.bed} {output.bed}; """
                """touch {output.bed_filt}; """
            )

# variant_caller_extern_filter_bed
#
# Apply region filter to variant BED.
rule variant_caller_extern_get_bed:
    input:
        bed=lambda wildcards: analib.variant_extern.get_bed(wildcards, config, [])
    output:
        bed=temp('temp/variant/caller/extern-{extern_source}/bed/{sample}/pre_filter/{vartype}_{svtype}.bed')
    wildcard_constraints:
        svtype='ins|del|inv|dup|sub|rgn|snv',
    run:

        bed_file_name = analib.variant_extern.get_bed(wildcards, config)

        if bed_file_name is None:
            print('Writing empty BED for extern {extern_source}: {vartype} - {svtype}'.format(**wildcards))

            if wildcards.vartype == 'snv':
                col_names = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'REF', 'ALT']
            else:
                col_names = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN']

            pd.DataFrame(
                [], columns=col_names
            ).to_csv(
                output.bed, sep='\t', index=False
            )

        else:

            # Get config
            extern_entry = analib.variant_extern.get_config(wildcards.extern_source, config)

            # Read
            df = pd.read_csv(input.bed, sep='\t', low_memory=False)

            if extern_entry['svtype_combined']:  # BED file contains a mix of svtypes, separate
                df = df.loc[df['SVTYPE'] == wildcards.svtype.upper()]

            # Sort
            if df.shape[0] > 0:
                df.sort_values(['#CHROM', 'POS'], inplace=True)

            # Write
            df.to_csv(output.bed, sep='\t', index=False)
