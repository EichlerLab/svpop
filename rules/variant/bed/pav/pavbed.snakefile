"""
Read from PAV BED files (not from VCF).
"""

# variant_pavbed_bed
#
# Get variants from PAV BED files.
rule variant_pavbed_bed:
    input:
        bed=lambda wildcards: os.path.join(
            variant_global_sample_table_entry('pavbed', wildcards=wildcards)['DATA'],
            '{vartype}_{svtype}.bed.gz'.format(**wildcards)
        ),
        fa=lambda wildcards: os.path.join(
            variant_global_sample_table_entry('pavbed',wildcards=wildcards)['DATA'],
            'fa/{vartype}_{svtype}.fa.gz'.format(**wildcards)
        )
    output:
        bed=temp('temp/variant/caller/pavbed-{seq_set}/{sample}/all/all/bed/pre_filter/{vartype}_{svtype}.bed.gz'),
        fa=temp('temp/variant/caller/pavbed-{seq_set}/{sample}/all/all/bed/pre_filter/fa/{vartype}_{svtype}.fa.gz')
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
            variant_global_sample_table_entry('pavbed', wildcards=wildcards)['DATA'],
            '{vartype}_{svtype}.bed.gz'.format(**wildcards)
        )
    output:
        bed=temp('temp/variant/caller/pavbed-{seq_set}/bed/{sample}/all/{filter}/{vartype}_{svtype}.bed.gz')
    wildcard_constraints:
        vartype='snv'
    shell:
        """cp {input.bed} {output.bed}; """
        """cp {input.fa} {output.fa}"""
