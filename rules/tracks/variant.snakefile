"""
Generate browser tracks.
"""

rule tracks_sv_bb:
    input:
        bed='temp/tracks/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/{vartype}_{svtype}.bed',
        asfile='temp/tracks/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/{vartype}_{svtype}.as',
        fai=config['reference_fai']
    output:
        bb='tracks/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/{vartype}_{svtype}.bb'
    shell:
        """bedToBigBed -tab -as={input.asfile} -type=bed9+ {input.bed} {input.fai} {output.bb}"""


# tracks_sv_bed
#
# Make track from variant BED.
rule tracks_sv_bed:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz',
        fai=config['reference_fai']
    output:
        bed=temp('temp/tracks/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/{vartype}_{svtype}.bed'),
        asfile=temp('temp/tracks/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/{vartype}_{svtype}.as')
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t', header=0)
        df_fai = svpoplib.ref.get_df_fai(input.fai)

        # Filter columns that have track annotations
        field_set = set(
            pd.read_csv(
                os.path.join(SVPOP_DIR, 'files/tracks/ucsc_track_fields.tsv'),
                sep='\t', header=0
            )['FIELD']
        )

        df = df.loc[:, [col for col in df.columns if col in field_set]]

        # Get track name and description for AS file
        track_name = 'VariantTable'
        track_description = '{sourcetype} - {sourcename} ({svset} / {filter}): {sample} - {vartype} {svtype}'.format(**wildcards)

        # Make track
        svpoplib.tracks.variant.make_bb_track(df, df_fai, output.bed, output.asfile, track_name, track_description)
