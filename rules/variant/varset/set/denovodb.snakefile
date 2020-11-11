"""
Variants from denovo-db. The denovo-db version is set in config['varset']['set']['denovodb']['version']
"""

def _varset_set_denovodb_full_tab(wildcards, config):
    """
    Get a path to the full table with all denovo-db annotations. Use "DNDB_ID" to link varset BED files back to records
    in this table.

    :param varset_name:
    :param wildcards:
    :param config:
    :return:
    """

    config_entry = analib.varset.get_config_entry('denovodb', wildcards, config)

    if 'full' not in config_entry:
        raise RuntimeError('Missing "full" in varset configuration entry: {}'.format(varset_name))

    if 'version' not in config_entry:
        raise RuntimeError('Missing "version" in varset configuration entry: {}'.format(varset_name))


    return config_entry['full'].format(
        version=config_entry['version']
    )


def _varset_set_denovodb_bed(wildcards, config):

    config_entry = analib.varset.get_config_entry('denovodb', wildcards, config)

    if 'bed' not in config_entry:
        raise RuntimeError('Missing "full" in varset configuration entry: denovodb')

    if 'version' not in config_entry:
        raise RuntimeError('Missing "version" in varset configuration entry: denovodb')

    return os.path.join(
        SVPOP_DIR,
        config_entry['bed'].format(
            version=config_entry['version'],
            vartype=wildcards.vartype,
            svtype=wildcards.svtype,
        )
    )

def _varset_set_denovodb_fa(wildcards, config):

    config_entry = analib.varset.get_config_entry('denovodb', wildcards, config)

    if 'fa' not in config_entry:
        raise RuntimeError('Missing "full" in varset configuration entry: denovodb')

    if 'version' not in config_entry:
        raise RuntimeError('Missing "version" in varset configuration entry: denovodb')


    return os.path.join(
        SVPOP_DIR,
        config_entry['fa'].format(
            version=config_entry['version'],
            vartype=wildcards.vartype,
            svtype=wildcards.svtype,
        )
    )

# varset_set_denovodb_bed
#
# Get denovo-db BED.
rule varset_set_denovodb_bed:
    input:
        bed=lambda wildcards: _varset_set_denovodb_bed(wildcards, config)
    output:
        bed=temp('temp/variant/varset/denovodb-{subset}/bed/{sample}/all/all/byref/{vartype}_{svtype}.bed')
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t')

        # Subset
        if wildcards.subset != 'all':
            df = df.loc[df['PHENO'].apply(lambda val: wildcards.subset in set(val.split(',')))]

        if wildcards.sample != 'all':
            df = df.loc[df['SAMPLE'].apply(lambda val: wildcards.sample in set(val.split(',')))]

        # Write
        df.to_csv(output.bed, sep='\t', index=False)
