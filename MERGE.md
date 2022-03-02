# SV-Pop merging and intersecting variants

SV-Pop has a variant merge and intersect system built into it for comparing variants.

1. Merging: Consolidating variants from multiple sources into one nonredundant callset
    * sampleset: Merge across samples.
    * callerset: Merge across callers.
1. Intersect: Generate a table of variants from two sources indicating which are found in both and which are exclusive
   to either source.

## Merging paradigm

SV-Pop has a flexible merging system and could support any number of merging strategies. Currently, it only has a
nonredundant strategy build in ("nr"). This document will describe parameters for the "nr" merge in SV-Pop.

We will likely support strategies provided by other tools Jasmine in the near future
(https://www.biorxiv.org/content/10.1101/2021.05.27.445886v1). 


## Nonredundant merging (nr)

The "nr" merge system built into SV-Pop is based on an order of input sources and merging is done iteratively
in that order. The first variant set (first caller or first sample) is copied into the merged set. The second
variant set is intersected with the merged set, matching calls are annotated with the intersected variants, and
new calls are added to the merged set. This process continues iteratively until all callers or samples are merged.

The result is that each variant in the merged set is represented by the call from the caller or sample where it was
first seen. Each merged variant contains enough information to trace back to the calls that support it; the list of
supporting callers/samples and variant IDs matching each are retained.

Example:

    results/variant/sampleset/freeze1/allsamples/all/all/bed/sv_ins.bed.gz

In this example, the configuration JSON will contain a definition for sampleset "freeze1", which includes variant
merging parameters. A set of samples is defined in a list called "allsamples" (defined in the `samplelist` section of
`config/config.json`).

### Comparing variants

Variant comparisons during merges and intersects occur in several tiers.

1. Exact: Variants matching exactly by size and position are considered first.
1. RO: Variants must overlap reciprocally by size and position.
1. Offset-size-RO: Variants within a set distance (offset) are compared by size whether or not they intersect.

RO is the classic way of comparing SVs before they were sequence-resolved and where breakpoint locations were less
certain. For two variants (A and B), A must overlap B by a certain percentage and vice-versa. It still works fairly
well for large SVs, but tends to under-merge small SVs and indels because small changes in the position cause RO to
decline rapidly. The typical cutoff is 50% RO. Insertions occur at a point, but the SVLEN is added to the position
to create a region for RO (only for the merge, the resulting variant BED file still has INS at a single point).

Offset-size was introduced to provide better merge performance for small SVs and indels. First, the offset is computed
(minimum of POS difference or END difference), and if the offset is in a certain distance (typically 200 bp), then
sizes are compared required a reciprocal size overlap (typically 50%). The size overlap is like RO if the variants were
shifted to maximum overlap (i.e. shift variant B POS to A's POS).

Variants passing one of these criteria are considered a match.


### Merge parameters

Merging is controlled by a string describing the parameters for the above merging process. The string begins with "nr:"
to indicate that SV-Pop's nonredundant merge system should be used. The remaining string a colon-separated list of
parameters and their values.

nr merge parameters:

| Parameter | Argument | Description |
| --- | --- | --- |
| ro | Percentage | Reciprocal overlap |
| szro | Percentage | Size reciprocal overlap |
| offset | Basepair distance | Variant offset |
| refalt | | Match REF and ALT (SNV) |
| ref | | Match REF (SNV) |
| alt | | Match ALT (SNV) |

Percentages (ro and szro) are a whole number (e.g. "50" means 50%).

The merge definitions we have used so far (Ebert 2021) are "nr:szro=50:offset=200" (SV and indel) or "nr:refalt" (SNVs). 

## Configuring merges

The configuration file has three sections that control merges.


### samplelist

Named lists of samples for sampleset merges. These names appear as "{sample}" wildcards.

Example:

    "samplelist": {
        "hgsvc2": [
            "HG00731", "HG00732",
            "HG00512", "HG00513",
            "NA19238", "NA19239",
            "NA24385", "HG03125",
            "NA12878", "HG03486", "HG02818",
            "HG03065", "HG03683", "HG02011",
            "HG03371", "NA12329", "HG00171",
            "NA18939", "HG03732", "HG00096",
            "NA20847", "HG03009", "NA20509",
            "HG00864", "HG01505", "NA18534",
            "NA19650", "HG02587", "HG01596",
            "HG01114", "NA19983", "HG02492",
            "HG00733", "HG00514", "NA19240"
        ],
        "PUR": ["HG00733", "HG00732", "HG00731"],
        "CHS": ["HG00514", "HG00513", "HG00512"],
        "YRI": ["NA19240", "NA19238", "NA19239"]
    }

### sampleset

Example sampleset section defining a sampleset merge ("pavhifi").

    "sampleset": {
        "pavhifi": {
            "sourcetype": "caller",
            "sourcename": "pav-hifi",
            "merge": {
                "sv,indel": "nr:szro=50:offset=200",
                "snv": "nr:refalt"
            },
            "name": "PAV HiFi",
            "description": "Pav HiFi"
        }

`sourcetype` and `sourcename` determine where samples are pulled from. The `merge` section can be a single merge
definition or a dictionary of merge strategies keyed by `svtype` (comma-separated list of svtypes is allowed). The
`name` will show up in figures generated for this merged callset. `description` is for documentation and is unused.

### callerset

Example callerset section.

    "callerset": {
        "longreads": {
            "callsets": [
                ["caller", "pav-hifi"],
                ["caller", "pav-clr"],
                ["caller", "pbsv-hifi"],
                ["caller", "pbsv-clr"]
            ],
            "name_list": ["PAVHIFI", "PAVCLR", "PBSVHIFI", "PBSVCLR"],
            "merge": {
                "sv,indel": "nr:szro=50:offset=200",
                "snv": "nr:refalt"
            },
            "name": "HGSVC2 LR",
            "description": "PAV/PBSV (HiFi & CLR)"
        }
    }

`callsets` contains a list variant sources (each is a souncetype/souncename pair). `name_list` must be as long as
`callsets` and contain short names for each source, which will be listed with each record in the merged variant table.
`merge` is the same as for `sampleset` (a single merge definition or a dictionary of definitions keyed by svtype(s)).
`name` is used for figures, and `description` is for documentation purposes and is unused.

## Merging annotations

Each merged variant has a list of samples and a list of variant IDs from each sample supporting the merged calls. This
allows annotations performed on individual callers and samples to be pulled through the merge process using the lead
variant for each merged record without re-running annotations.

For example, to get the RefSeq intersections for a merged set:

    results/variant/sampleset/freeze1/allsamples/all/all/anno/refseq/refseq-count_sv_ins.tsv.gz


## Variant intersects

Instead of generating an fully merged callset, it is often necessary to match variants across multiple sets.
For example, there might by many orthogonal callsets supporting a main callset, but generating a sampleset merge across
all the orthogonal callsets and the main callset is extremely messy and should be avoided. Instead, the main callset
can be intersected with other callsets with each intersect generating a table of variants from both showing which
were in both samples and which did not intersect the other sample.

### merge_def

The configuration file contains a `merge_def` section to define aliases for merge strings. The intersect file (see
below) contains the merge definition, which can be long and messy, and `merge_def` alias provide a way to clean it up.

For example:

    "merge_def": {
        "szro-50-200": "nr:szro=50:offset=200",
        "snvex": "nr:refalt",
    }

### Intersect file

Pattern:

    results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/intersect.tsv.gz

1. sourcetype_a, sourcename_a, sample_a: sourcetype, sourcename, and sample for the first variant set.
1. sourcetype_b, sourcename_b, sample_b: sourcetype, sourcename, and sample for the second variant set.
1. merge_def: Merge definition.
1. filter, svset, vartype, and svtype: Same wildcards for variants. 

Example:

    results/variant/intersect/caller+pav+HG00733/caller+pbsv+HG00733/szro-50-200/all/all/sv_ins/intersect.tsv.gz

This compares two callsets, a PAV callset and a pbsv callset on HG00733 with an RO/size-RO strategy
("szro-50-200", defined in `merge_def` in the config, see above). This format allows flexibly comparisons, for example,
comparing two samples or two callers. It is not restricted to callers, and samples could be compared to a whole merged
callset.
