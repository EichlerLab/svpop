"""
Flags for intersects.
"""

def _input_svpop_flag_intersect_input_list(input_pattern, wildcards, policy):

    wc_dict = dict(wildcards)

    input_list = list()

    for sample_a, sample_b in _input_svpop_flag_intersect_sample_tuples(wildcards, policy):
        wc_dict['sample_a'] = sample_a
        wc_dict['sample_b'] = sample_b

        input_list.append(input_pattern.format(**wc_dict))

    return input_list

def _input_svpop_flag_intersect_sample_tuples(wildcards, policy):

    # Get sample lists
    sample_list_a = svpoplib.rules.get_sample_list(wildcards.sample_a, config)
    sample_list_b = svpoplib.rules.get_sample_list(wildcards.sample_b, config)

    if sample_list_a is None:
        sample_list_a = [wildcards.sample_a]

    if sample_list_b is None:
        sample_list_b = [wildcards.sample_a]

    len_a = len(sample_list_a)
    len_b =  len(sample_list_b)

    if policy == 'match':
        return [
            (
                sample_list_a[i % len_a],
                sample_list_b[i % len_b]
            ) for i in range(max([len_a, len_b]))
        ]

    elif policy == 'comb':
        return [
            (sample_a, sample_b)
                for sample_a in sample_list_a for sample_b in sample_list_b if sample_a != sample_b
        ]

    elif policy == 'fullcomb':
        return [
            (sample_a, sample_b)
                for sample_a in sample_list_a for sample_b in sample_list_b
        ]

    else:
        raise RuntimeError(f'Unrecognized sample policy: {policy}')


#
# Rules
#

localrules: flag_intersect

# flag_intersect
#
# Generate intersects
rule flag_intersect:
    input:
        tsv=lambda wildcards: _input_svpop_flag_intersect_input_list(
            'results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/intersect.tsv.gz',
            wildcards, wildcards.policy
        )
    output:
        flag=touch('flag/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/{policy}')
    wildcard_constraints:
        policy='match|comb|fullcomb'
