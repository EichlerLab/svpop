"""
Global annotation rules.
"""

# variant_caller_anno_global_varset
#
# Filter annotations by variant class.
rule variant_caller_anno_global_varset:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/bed/{sample}/{varset}/{filter}/byref/{vartype}_{svtype}.bed',
        tab='results/variant/{sourcetype}/{sourcename}/anno/{sample}/all/{filter}/{annodir}/{annotype}_{vartype}_{svtype}.{ext}'
    output:
        tab='results/variant/{sourcetype}/{sourcename}/anno/{sample}/{varset}/{filter}/{annodir}/{annotype}_{vartype}_{svtype}.{ext}'
    wildcard_constraints:
        varset='((?!all).*)|all.+',
        svtype='ins|del|inv|dup|snv|rgn|sub'
    run:

        # Read
        df = pd.read_csv(input.tab, sep='\t')

        id_set = set(pd.read_csv(input.bed, sep='\t', usecols=('ID', ), squeeze=True))

        # Filter
        df = df.loc[df['ID'].apply(lambda val: val in id_set)]

        # Write
        df.to_csv(output.tab, sep='\t', index=False)
