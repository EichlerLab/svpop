# Flag rules. Generates a collection of data, but does not produce anything.
# Example: Generate a caller for a list of samples, but do not merge.


#
# Definitions
#

def _input_svpop_variant_flag_input_list(wildcards, file_pattern):

    # Get sample list
    sample_list = svpoplib.rules.get_sample_list(wildcards.sample, config)

    if sample_list is None:
        raise RuntimeError(f'No sample list with name: {wildcards.sample}')

    # Protect all wildcards except "{sample}" from formatting (i.e. "/path/{sample}/{vartype}_{svtype}.ext" to "/path/{sample}/{{vartype}}_{{svtype}}.ext")
    file_pattern = re.sub(r'\{([^}]+)\}(?<!{sample})', r'{{\1}}', file_pattern)

    return [file_pattern.format(sample=sample) for sample in sample_list]


#
# Rules
#

# svpop_variant_flag_bed
#
# Generate the callset for a collection of samples.
localrules: flag_variant_caller_bed

rule flag_variant_caller_bed:
    input:
        bed=lambda wildcards: _input_svpop_variant_flag_input_list(
            wildcards,
            'results/variant/caller/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz'
        )
    output:
        bed=touch('flag/variant/caller/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz')

# flag_variant_caller_anno
#
# Generate annotations for a collection of samples.
localrules: flag_variant_caller_anno

rule flag_variant_caller_anno:
    input:
        anno=lambda wildcards: _input_svpop_variant_flag_input_list(
            wildcards,
            'results/variant/caller/{sourcename}/{sample}/{filter}/{svset}/anno/{annodir}/{annotype}_{vartype}_{svtype}.{ext}.gz'
        )
    output:
        anno=touch('flag/variant/caller/{sourcename}/{sample}/{filter}/{svset}/anno/{annodir}/{annotype}_{vartype}_{svtype}.{ext}.gz')

# flag_variant_caller_bed_fa
#
# Get callset BED
localrules: flag_variant_caller_bed_fa

rule flag_variant_caller_bed_fa:
    input:
        bed=lambda wildcards: _input_svpop_variant_flag_input_list(
            wildcards,
            'results/variant/caller/{sourcename}/{sample}/{filter}/{svset}/bed/fa/{vartype}_{svtype}.fa.gz'
        )
    output:
        bed=touch('flag/variant/caller/{sourcename}/{sample}/{filter}/{svset}/bed/fa/{vartype}_{svtype}.fa.gz')
