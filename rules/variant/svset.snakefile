"""
Filter variants based on variant characteristics.
"""

# variant_svsetfilter_run
#
# Filter by SV set.
rule variant_svsetfilter_run:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        filter_files=lambda wildcards: analib.svset.get_filter_input_files(wildcards.svset, wildcards, config)
    output:
        bed='results/variant/{sourcetype}/{sourcename}/bed/{sample}/{svset}/{filter}/byref/{vartype}_{svtype}.bed'
    wildcard_constraints:
        svset='((?!all).*|all.+)'  # Do not allow "all"
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t')

        # Filter
        df = analib.svset.apply_svset_filter(
            pd.read_csv(input.bed, sep='\t'),
            wildcards.svset,
            wildcards,
            config
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False)
