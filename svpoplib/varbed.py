"""
Utilities for handling variants in BED format.
"""

import collections
import pandas as pd
import re


def bcftools_query_to_tsv(df, sample, strict_sample=False, filter_gt=True, multi_gt_override=False):
    """
    Process a table output from `bcftools query`. Correct column names and extract a target sample.

    :param df: DataFrame of records from bcftools query.
    :param sample: Sample to extract if VCF had multiple samples.
    :param strict_sample: If the VCF has only one sample, require the same name to match the `sample` field. Otherwise,
        accept any sample name (default).
    :param filter_gt: Filter by the GT column. Only variants with an alternate allele are accepted. If GT is missing and
        there was only one sample, filtering on GT is skipped. If `filter_gt` is False and there is more than one
        sample, an error is output unless `multi_gt_override` is `True`; this will output the same variant callset for
        all samples in the VCF and ignore genotypes.

    :return: Formatted dataframe containing variant calls for one sample.
    """

    # Remove column formatting
    df.columns = [re.sub(r'^(#?)\s*\[[^\]]+\](.*)$', r'\1\2', col) for col in df.columns]

    if len(set(df.columns)) != df.shape[1]:
        dup_cols = [col for col, count in collections.Counter(df.columns).items() if count > 1]

        raise RuntimeError('Fonud duplicate column names before sample filtering: {}'.format(', '.join(sorted(dup_cols))))

    # Get sample sets
    vcf_sample_set = {val.rsplit(':', 1)[0] for val in df.columns if ':' in val}
    vcf_sample_count = len(vcf_sample_set)

    # Check number of samples
    if not filter_gt and vcf_sample_count > 1 and not multi_gt_override:
        raise RuntimeError(f'Must filter by genotype if more than one sample is present in the VCF (found {vcf_sample_count}): Set "multi_gt_override=True" to override and output the same calls for all samples in the VCF')

    # Adjust column labels and select target sample
    if vcf_sample_count > 1:

        # Multiple samples in VCF
        if sample not in vcf_sample_set:
            vcf_sample_list_str = ', '.join(sorted(vcf_sample_set))

            raise RuntimeError(f'Found {vcf_sample_count} samples in VCF and none match target sample "{sample}": {vcf_sample_list_str}')

        target_sample = sample

    else:
        # zero or one sample

        # Fail if strict_samples and sample name does not match
        if strict_sample and sample not in vcf_sample_set:
            vcf_sample_list_str = ', '.join(sorted(vcf_sample_set))
            raise RuntimeError(f'Found {vcf_sample_count} samples in VCF and none match target sample "{sample}": {vcf_sample_list_str}')

        if vcf_sample_count == 1:
            target_sample = list(vcf_sample_set)[0]
        else:
            target_sample = None  # No sample columns

    # Subset to specific sample
    df = df[[col for col in df.columns if ':' not in col or col.rsplit(':', 1)[0] == target_sample]]

    df.columns = [col if ':' not in col else col.rsplit(':', 1)[1] for col in df.columns]

    if len(set(df.columns)) != df.shape[1]:
        dup_cols = [col for col, count in collections.Counter(df.columns).items() if count > 1]

        raise RuntimeError('Fonud duplicate column names after sample filtering: {}'.format(', '.join(sorted(dup_cols))))

    if filter_gt:
        if 'GT' in df.columns:
            df = df.loc[df['GT'].apply(gt_has_alt)]

        elif vcf_sample_count > 1:
            raise RuntimeError(f'Cannot filter on GT column for a VCF with {vcf_sample_count} samples: No GT columns to filter')

    # Return formatted DataFrame
    return df


def gt_has_alt(gt):
    """
    Parse a genotype (GT field) separating on "/" or "|" and output `True` if any  field is an integer greater than 1.

    :param gt: Genotype.

    :return: `True` if any allele is an alternate and `False` if all alleles are reference and/or missing.
    """

    if gt is None or pd.isnull(gt):
        return False

    for val in re.split(r'\||\/', gt):
        if val == '.':
            continue

        try:
            if int(val) > 0:
                return True

        except ValueError:
            pass  # Ignore

    return False

