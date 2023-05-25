"""
Read from PAV BED files (not from VCF).
"""

# variant_pavbed_bed
#
# Get variants from PAV BED files.
rule variant_pavbed_bed:
    input:
        bed=lambda wildcards: os.path.join(
            svpoplib.rules.sample_table_entry(wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='pavbed')['DATA'],
            '{vartype}_{svtype}.bed.gz'.format(**wildcards)
        ),
        fa=lambda wildcards: os.path.join(
            svpoplib.rules.sample_table_entry(wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='pavbed')['DATA'],
            'fa/{vartype}_{svtype}.fa.gz'.format(**wildcards)
        )
    output:
        bed=temp('temp/variant/caller/pavbed/{sourcename}/{sample}/all/all/bed/pre_filter/{vartype}_{svtype}.bed.gz'),
        fa=temp('temp/variant/caller/pavbed/{sourcename}/{sample}/all/all/bed/pre_filter/fa/{vartype}_{svtype}.fa.gz')
    wildcard_constraints:
        vartype='sv|indel'
    shell:
        """cp {input.bed} {output.bed}; """
        """cp {input.fa} {output.fa}"""

# variant_pavbed_bed_snv
#
# Get SNVs from PAV BED files.
rule variant_pavbed_bed_snv:
    input:
        bed=lambda wildcards: os.path.join(
            svpoplib.rules.sample_table_entry(wildcards.sourcename, SAMPLE_TABLE, wildcards=wildcards, caller_type='pavbed')['DATA'],
            '{vartype}_{svtype}.bed.gz'.format(**wildcards)
        )
    output:
        bed=temp('temp/variant/caller/pavbed/{sourcename}/{sample}/all/all/bed/pre_filter/{vartype}_{svtype}.bed.gz'),
        fa=temp('temp/variant/caller/pavbed/{sourcename}/{sample}/all/all/bed/pre_filter/fa/{vartype}_{svtype}.fa.gz')
    wildcard_constraints:
        svtype='snv',
        vartype='snv'
    shell:
        """cp {input.bed} {output.bed}; """
        """> {output.fa}"""
