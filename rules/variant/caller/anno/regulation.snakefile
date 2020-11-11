#
# ORegAnno
#

# variant_anno_reg_oreganno_intersect
#
# Intersect variants with ORegAnno.
rule variant_anno_reg_oreganno_intersect:
    input:
        bed='results/variant/caller/{caller}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        anno='data/anno/oreganno/oreganno.bed'
    output:
        tab='results/variant/caller/{caller}/anno/{sample}/all/{filter}/oreganno/oreganno_{vartype}_{svtype}.tab'
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    shell:
        """{{\n"""
        """    {{\n"""
        """        head -n 1 {input.bed};\n"""
        """        head -n 1 {input.anno} | sed -re 's/#|\s+/\\tOREG_/g' | sed -re 's/^\s+//'\n"""
        """    }} | awk -vORS="\t" '{{print}}' | sed -re 's/\\s+/\\t/g' | sed -re 's/\\s*$/\\n/';\n"""
        """    bedtools intersect -a {input.bed} -b {input.anno} -loj\n"""
        """}} | """
        """awk -vFS="\t" -vOFS="\t" '"""
            """(NR == 1) {{for (i = 1; i <= NF; ++i) {{ix[$i] = i}}; print $ix["ID"], $ix["OREG_ID"]}} """
            """(NR > 1 && $ix["OREG_ID"] != ".") {{print $ix["ID"], $ix["OREG_ID"]}} """
        """' """
        """> {output.tab}"""

# Transitioning to merge_variants method (Stashed using comments on 2019-07-31, finish or remove).
#
# rule variant_anno_reg_oreganno_intersect:
#     input:
#         bed='results/variant/caller/{caller}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
#         anno='data/anno/oreganno/oreganno.bed'
#     output:
#         tab='results/variant/caller/{caller}/anno/{sample}/all/{filter}/oreganno/oreganno_{merge_def}_{vartype}_{svtype}.tab'
#     wildcard_constraints:
#         svtype='ins|del|inv|dup|snv'
#     params:
#         cpu=4,
#         mem='1G',
#         rt='8:00:00'
#     run:
#
#         config_def = analib.svmerge.get_merge_def(wildcards.merge_def, config)
#
#         # Merge
#         df = analib.svmerge.merge_variants(
#             bed_list=[input.bed, input.anno],
#             sample_names=['ID', 'OREGANNO_ID'],
#             strategy=config_def,
#             threads=params.cpu
#         )
#
#         # Subset columns
#         df = df.loc[:, ['ID', 'MERGE_SAMPLES', 'MERGE_SRC', 'MERGE_VARIANTS']]
#
#         # Create an ID column for sample (empty string if the variant was not in that sample)
#         df.loc[df['MERGE_SRC'] == 'A', 'MERGE_VARIANTS'] =  df.loc[df['MERGE_SRC'] == 'A', 'MERGE_VARIANTS'] + ','
#         df.loc[df['MERGE_SRC'] == 'B', 'MERGE_VARIANTS'] =  ',' + df.loc[df['MERGE_SRC'] == 'B', 'MERGE_VARIANTS']
#
#         df['ID_A'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[0])
#         df['ID_B'] = df['MERGE_VARIANTS'].apply(lambda val: val.split(',')[1])
#
#         # Subset and write
#         df['SOURCE_SET'] = df['MERGE_SAMPLES']
#         df = df.loc[:, ['ID_A', 'ID_B', 'SOURCE_SET']]
#
#         df.to_csv(output.tab, sep='\t', index=False)

#
# Encode
#

# variant_anno_reg_encode_histone
#
# Intersect variants with histone encode marks (histone modifications)
rule variant_anno_reg_encode_histone:
    input:
        bed='results/variant/caller/{caller}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        anno='data/anno/encode/threshold/{mark}_all_{threshold}.bed'
    output:
        tab='results/variant/caller/{caller}/anno/{sample}/all/{filter}/encode/encode-{mark}-{threshold}-{overlap}_{vartype}_{svtype}.tab'
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    run:

        # Get overlap parameters
        if wildcards.overlap.lower() == 'any':
            overlap_param = ''
        else:
            overlap_param = '-f {wildcards.overlap}'

        # Intersect
        shell(
            """echo "ID" > {output.tab}; """
            """bedtools intersect -a {input.bed} -b {input.anno} -wa {overlap_param} -u | """
            """cut -f4 """
            """>> {output.tab}"""
        )

#
# Candidate cis-regulatory elements (CCRE2020) - ENCODE 2020
#

# variant_anno_reg_dnase_cluster
#
# Intersect variants with DNAse clusters.
rule variant_anno_reg_ccre2020:
    input:
        bed='results/variant/caller/{caller}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        anno='data/anno/ccre2020/ccre_{score}.bed.gz'
    output:
        tab='results/variant/caller/{caller}/anno/{sample}/all/{filter}/ccre/ccre2020-{score}_{vartype}_{svtype}.tab'
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    shell:
        """\n{{\n"""
        """    {{\n"""
        """        head -n 1 {input.bed};\n"""
        """        echo $(zcat {input.anno} | head -n 1)\n"""
        """    }} | awk -vORS="\t" '{{print}}' | sed -re 's/\\s+/\\t/g' | sed -re 's/\\s*$/\\n/';\n"""
        """    bedtools intersect -a {input.bed} -b {input.anno} -loj\n"""
        """}} | """
        """awk -vOFS="\t" '"""
            """(NR == 1) {{for (i = 1; i <= NF; ++i) {{ix[$i] = i}}; print $ix["ID"], $ix["CCRE2020"], $ix["EH38D"], $ix["EH38E"]}} """
            """(NR > 1 && $ix["EH38D"] != ".") {{print $ix["ID"], $ix["CCRE2020"], $ix["EH38D"], $ix["EH38E"]}} """
        """' """
        """> {output.tab}"""

#
# DNAse hypersensitivity clusters - ENCODE 2020
#

# variant_anno_reg_dnase_cluster
#
# Intersect variants with DNAse clusters.
rule variant_anno_reg_dnase2020:
    input:
        bed='results/variant/caller/{caller}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        anno='data/anno/dhs2020/dhs_{score}.bed.gz'
    output:
        tab='results/variant/caller/{caller}/anno/{sample}/all/{filter}/dhs/dhs2020-{score}_{vartype}_{svtype}.tab'
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    shell:
        """\n{{\n"""
        """    {{\n"""
        """        head -n 1 {input.bed};\n"""
        """        echo $(zcat {input.anno} | head -n 1)\n"""
        """    }} | awk -vORS="\t" '{{print}}' | sed -re 's/\\s+/\\t/g' | sed -re 's/\\s*$/\\n/';\n"""
        """    bedtools intersect -a {input.bed} -b {input.anno} -loj\n"""
        """}} | """
        """awk -vOFS="\t" '"""
            """(NR == 1) {{for (i = 1; i <= NF; ++i) {{ix[$i] = i}}; print $ix["ID"], $ix["DHS_ID"], $ix["DHS_SIG_MEAN"], $ix["MOTIF_CLUSTERS"]}} """
            """(NR > 1 && $ix["DHS_ID"] != ".") {{print $ix["ID"], $ix["DHS_ID"], $ix["DHS_SIG_MEAN"], $ix["MOTIF_CLUSTERS"]}} """
        """' """
        """> {output.tab}"""

#
# DNAse hypersensitivity clusters (Encode) - UCSC browser
#

# variant_anno_reg_dnase_cluster
#
# Intersect variants with DNAse clusters.
rule variant_anno_reg_dnase_cluster:
    input:
        bed='results/variant/caller/{caller}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        anno='data/anno/dhs_cluster/dhs_cluster_{score}.bed'
    output:
        tab='results/variant/caller/{caller}/anno/{sample}/all/{filter}/dhs/dhs-{score}_{vartype}_{svtype}.tab'
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    shell:
        """{{\n"""
        """    {{\n"""
        """        head -n 1 {input.bed};\n"""
        """        head -n 1 {input.anno}\n"""
        """    }} | awk -vORS="\t" '{{print}}' | sed -re 's/\\s+/\\t/g' | sed -re 's/\\s*$/\\n/';\n"""
        """    bedtools intersect -a {input.bed} -b {input.anno} -loj\n"""
        """}} | """
        """awk -vOFS="\t" '"""
            """(NR == 1) {{for (i = 1; i <= NF; ++i) {{ix[$i] = i}}; print $ix["ID"], $ix["DHS_NAME"], $ix["DHS_SOURCE_IDS"], $ix["DHS_SOURCE_SCORES"]}} """
            """(NR > 1 && $ix["DHS_NAME"] != ".") {{print $ix["ID"], $ix["DHS_NAME"], $ix["DHS_SOURCE_IDS"], $ix["DHS_SOURCE_SCORES"]}} """
        """' """
        """> {output.tab}"""


#
# CpG Islands
#

# variant_anno_reg_cpgisland
#
# Intersect variants with DNAse clusters.
rule variant_anno_reg_cpgisland:
    input:
        bed='results/variant/caller/{caller}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        anno='data/anno/cpgisland/cpgisland_{mask}.bed'
    output:
        tab='results/variant/caller/{caller}/anno/{sample}/all/{filter}/cpgisland/cpgisland-{mask}_{vartype}_{svtype}.tab'
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    shell:
        """{{\n"""
        """    {{\n"""
        """        head -n 1 {input.bed};\n"""
        """        head -n 1 {input.anno}\n"""
        """    }} | awk -vORS="\t" '{{print}}' | sed -re 's/\\s+/\\t/g' | sed -re 's/\\s*$/\\n/';\n"""
        """    bedtools intersect -a {input.bed} -b {input.anno} -loj\n"""
        """}} | """
        """awk -vOFS="\t" '"""
            """(NR == 1) {{for (i = 1; i <= NF; ++i) {{ix[$i] = i}}; print $ix["ID"], $ix["LENGTH"], $ix["CPG_NUM"], $ix["GC_NUM"], $ix["CPG_PCT"], $ix["GT_PCT"], $ix["OBS_EXP"]}} """
            """(NR > 1 && $ix["CPG_NUM"] != ".") {{print $ix["ID"], $ix["LENGTH"], $ix["CPG_NUM"], $ix["GC_NUM"], $ix["CPG_PCT"], $ix["GT_PCT"], $ix["OBS_EXP"]}} """
        """' """
        """> {output.tab}"""

#
# All regulatory
#

# variant_anno_reg_all_reg
#
# Merge all regulatory elements into one table.
rule variant_anno_reg_all_reg:
    input:
        bed='results/variant/caller/{caller}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        cpg='results/variant/caller/{caller}/anno/{sample}/all/{filter}/cpgisland/cpgisland-unmasked_{vartype}_{svtype}.tab',
        h3k27ac='results/variant/caller/{caller}/anno/{sample}/all/{filter}/encode/encode-H3K27Ac-50-any_{vartype}_{svtype}.tab',
        h3k4me3='results/variant/caller/{caller}/anno/{sample}/all/{filter}/encode/encode-H3K4Me3-50-any_{vartype}_{svtype}.tab',
        h3k4me1='results/variant/caller/{caller}/anno/{sample}/all/{filter}/encode/encode-H3K4Me1-50-any_{vartype}_{svtype}.tab',
        dhs='results/variant/caller/{caller}/anno/{sample}/all/{filter}/dhs/dhs-200_{vartype}_{svtype}.tab',
        oreganno='results/variant/caller/{caller}/anno/{sample}/all/{filter}/oreganno/oreganno_{vartype}_{svtype}.tab'
    output:
        tab='results/variant/caller/{caller}/anno/{sample}/all/{filter}/allreg/allreg_{vartype}_{svtype}.tab'
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    run:

        # Read SVs
        df = pd.read_csv(input.bed, sep='\t', header=0, usecols=('ID',), squeeze=False)

        df['CPG_ISLAND'] = False
        df['ENCODE_H3K27AC'] = False
        df['ENCODE_H3K4ME3'] = False
        df['ENCODE_H3K4ME1'] = False
        df['DHS'] = False
        df['OREGANNO'] = False

        df.set_index('ID', inplace=True)

        # Add annotations
        df.loc[set(pd.read_csv(input.cpg, sep='\t', header=0, usecols=('ID',), squeeze=True)), 'CPG_ISLAND'] = True
        df.loc[set(pd.read_csv(input.h3k27ac, sep='\t', header=0, usecols=('ID',), squeeze=True)), 'ENCODE_H3K27AC'] = True
        df.loc[set(pd.read_csv(input.h3k4me3, sep='\t', header=0, usecols=('ID',), squeeze=True)), 'ENCODE_H3K4ME3'] = True
        df.loc[set(pd.read_csv(input.h3k4me1, sep='\t', header=0, usecols=('ID',), squeeze=True)), 'ENCODE_H3K4ME1'] = True
        df.loc[set(pd.read_csv(input.dhs, sep='\t', header=0, usecols=('ID',), squeeze=True)), 'DHS'] = True
        df.loc[set(pd.read_csv(input.oreganno, sep='\t', header=0, usecols=('ID',), squeeze=True)), 'OREGANNO'] = True

        # Clear empty records
        df = df.loc[df.apply(any, axis=1)]

        # Write
        df.to_csv(output.tab, sep='\t', index=True)

# variant_anno_reg_all_reg_summary
#
# Collapse matrix of boolean values (SV x Annotation) into a list of annotations for each SV (where
# its field is True).
rule variant_anno_reg_all_reg_summary:
    input:
        tab='results/variant/caller/{caller}/anno/{sample}/all/{filter}/allreg/allreg_{vartype}_{svtype}.tab'
    output:
        tab='results/variant/caller/{caller}/anno/{sample}/all/{filter}/allreg/allreg-list_{vartype}_{svtype}.tab'
    run:

        # Setup translation table
        header_trans = {
            'CPG_ISLAND': 'CpG',
            'ENCODE_H3K27AC': 'H3K27Ac',
            'ENCODE_H3K4ME3': 'H3K4Me3',
            'ENCODE_H3K4ME1': 'H3K4Me1',
            'DHS': 'DHS',
            'DHS_STAM': 'DHS(Stam)',
            'OREGANNO': 'ORegAnno'
        }

        # Read variants
        df = pd.read_csv(input.tab, sep='\t', header=0, index_col='ID')

        # Translate to a comma-separated list (changing column names to formatted names)
        name_list = pd.Series([header_trans[col] for col in df.columns])

        reg_names = df.apply(lambda row: ','.join(name_list.loc[list(row)]), axis=1)
        reg_names.name = 'REG_ELEMENTS'

        # Write
        reg_names.to_csv(output.tab, sep='\t', index=True, header=True)


# variant_anno_reg_cns_dhs
#
# Intersect with fetal CNS DHS sites.
rule variant_anno_reg_cns_dhs:
    input:
        bed='results/variant/caller/{caller}/bed/{sample}/all/{filter}/byref/{vartype}_{svtype}.bed',
        dhs_bed=lambda wildcards: os.path.join(SVPOP_DIR, 'files/anno/fetal_cns_dhs/CNS_DHS_b38.bed.gz')
    output:
        tab='results/variant/caller/{caller}/anno/{sample}/all/{filter}/dhs/dhsfetal_{vartype}_{svtype}.tab'
    wildcard_constraints:
        svtype='ins|del|inv|dup|snv|rgn|sub'
    shell:
        """{{\n"""
        """    echo "ID"; """
        """    zcat {input.dhs_bed} | \n"""
        """    bedtools intersect -a {input.bed} -b stdin | \n"""
        """    cut -f4;\n"""
        """}} > {output.tab}"""

