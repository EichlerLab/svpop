# Write callsets as a VCF

import Bio.bgzf
import pandas as pd

import svpoplib

global shell

# vcf_write_vcf
#
# Make VCF headers.
rule vcf_write_vcf:
    input:
        bed='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz',
        fa='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/fa/{vartype}_{svtype}.fa.gz',
        ref_tsv='data/ref/contig_info.tsv.gz'
    output:
        vcf='results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/vcf/{vartype}_{svtype}_{altfmt}.vcf.gz'
    wildcard_constraints:
        altfmt='alt|sym|sym-noseq'
    run:

        # Read reference
        df_ref = pd.read_csv(input.ref_tsv, sep='\t')

        # Get VCF table and headers
        vcf_tab = svpoplib.vcf.VariantVcfTable(
            pd.read_csv(input.bed, sep='\t'),
            wildcards.sample,
            config['reference'],
            wildcards.altfmt
        )

        # Write VCF
        with Bio.bgzf.open(output.vcf, 'wt') as out_file:
            vcf_tab.write(df_ref, out_file)

        # Write tabix index if possible
        try:
            shell("""tabix {output.vcf} && touch -r {output.vcf} {output.vcf}.tbi""")
        except:
            pass
