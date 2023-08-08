# SV-Pop merging and intersecting variants

SV-Pop has a variant merge and intersect system built into it for comparing variants.

1. Merging: Consolidating variants from multiple sources into one nonredundant callset
    * sampleset: Merge across samples.
    * callerset: Merge across callers.
1. Intersect: Generate a table of variants from two sources indicating which are found in both and which are exclusive
   to either source.

For samplesets, SVs and indels are merged together to avoid artifacts on SV boundaries (i.e. missing support for a 50
bp SV from a 49 bp indel).


## Merging paradigm

SV-Pop has a flexible merging system and could support any number of merging strategies. Currently, it only has a
nonredundant strategy build in ("nr", or "nrsnv" for SNVs to enforce REF/ALT matches). This document will describe
parameters for the "nr" merge in SV-Pop.

We will likely support strategies provided by other tools, such as Jasmine
(https://www.biorxiv.org/content/10.1101/2021.05.27.445886v1) and Truvari
(https://www.biorxiv.org/content/10.1101/2022.02.21.481353v1.abstract).


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

### Comparison strategies

Several strategies are supported by SV-Pop's nonredundant merge process. Square brackets contain the merge strategy
identifier that is used in merge configuration (next section).

1. Exact [exact]: Size and position must match exactly.
1. Reciprocal overlap [ro]: Variants must overlap reciprocally by size and position. This is a classic merge strategy
   that works well for very large variants. For insertions, the length is added to the start position to obtain the
   end position for this strategy.
1. Offset-size reciprocal overlap [szro]: Variants within a set distance (offset) are compared by size whether or not
   they intersect. The offset can be restricted by either a set number of bases or as a multiplier of the variant size
   (or both). This strategy works better for smaller variants where small changes is position create large differences
   by ro.
1. Distance [distance]: A strict distance match that ignores the variant size. This strategy is useful for comparing
   callsets with extreme variant size inaccuracies, typcially legacy short-read SV callsets (modern short-read callsets
   should have more accurate sizes).

The merging process can use one or more of these strategies in order. For example, merging by "ro" then "szro" would
match everything by reciprcal overlep, and anything unmatched would then be matched by size-offset recprocal overlap.

### Sequence match similarity

The nonredundant merging process can compare variant sequences if sequence-resolved callests are being compared (i.e.
the sequence of each variant is known) using "match" parameters. Comparing sequences reduces over-merging by requiring
the sequence content of matched variants be similar (e.g. 80% similarity). Sequence-resolved callsets cannot be
generated from a VCF with symbolic ALTs (e.g. "<INS>") unless there is a SEQ INFO field for each containing the full
variant sequence. Many modern tools, especially from long reads, provide appropriate callsets for match similarity

The match similarity can be performed in two ways:
1. Similarity by alignment
1. Similarity by k-mer Jaccard

Alignment similarity is achieved by aligning both sequences with a Smith-Waterman alignment with special code that
allows it to remain sensitive when shifted variants alter the sequence inside the SV.

Jaccard similarity is achieved by taking k-mers (default 9-mer, from Jasmine's default) from both sequences and
comparing: A&B / (A!B + A&B + !AB).

Similarity by alignment is the most accurate comparison method, but it becomes very slow for large variants. For
HiFi callsets, k-mer Jaccard and alignment produce equivalent similarity scores above 4 kbp.


### Merge parameter syntax

Merging is controlled by a grammar starting with "nr::" followed by a colon-separated list of comparison strategies
each with optional arguments in parenthesis. If no merge strategies are given, it defaults to exact merging without
sequence matches.

The "match" parameter may appear in the configuration string where it sets a default match for all strategies. It may
also appear as an argument within a strategy where it overrides the default.

Examples:
1. "nr::" or "nr": Exact match all variants by size and positon, no match.
1. "nr::exact": Exact match all variants by size and positon with match.
1. "nr::ro:szro:exact": Match by ro, then by szro. Apply match parameters to both ro and szro.

The default values may be tuned for each step.

Example:
1. "nr::ro(0.8):match": Match by reciprocal overlap including default sequence match parameters.

Each parameter has a list of parameters with default values for each. Arguments may be supplied in the same order as
the list, with keywords, or both. Like Python, once a keyword argument is found, all subsequent arguments for that
strategy must be keyword arguments. Both empty parethesis and no parenthesis indicate that all default parameters are
used.

Example:
1. "nr::ro(ro=0.8):szro(0.2,szdist=2):match()"

Match may appear inside or outside of strategies. This is typically not needed, but SV-Pop allows the flexibility.

Examples:
1. "nr::ro(0.8,match(0.75)):szro:match"

This example uses the default match parameters for szro (and other strategies if they were present) and overrides the
default to require 75% match for ro.

#### Variant offsets

Several merges can be restricted by breakpoint offsets. SV-Pop defines offset as the minimum distance between
start positions and end positions. By taking the minimum of start and end positions, the distance parameter does not
over-penalize distance for variants of differing sizes.

#### Numeric types

Each parameter has a set of types it can accept.

Floating point (float):
1. Any unqualified number with a dot
   * 500.0, -5.0, 10.05, 0.002
1. Any exponential number with "e" or "E":
   * 5.0e2, -5e0, 1.005E1, 2e-3

Integer (int):
1. Any number without a dot
   * 500, -5, 110000000
1. Any number with a multiplier: k, m, g (not case sensitive)
   * 0.5k, 110m, 3.1G

The "unlimited" keyword may be substituted for a value to indicate that there is no limit. Not all paramers accept
the unlimited keyword.

Several parameter types may accept a mix of int, float, and unlimited. These are used in the parameter documentation
below.

1. num: float, int,
1. num_u: float, int, unlimited
1. int_u: int, unlimited
1. float_u: float, unlimited


#### Merge parameter: ro

Match by reciprocal-overlap (ro). For two variants to pass, their size and position must overlap by at least some
percentage. Since insertions are represented as a single point in the reference, the insertion length is added to the
start position for the RO calculation (makes it equivalent to deletions). This is the classic intersect strategy
(typically 50% RO), and it works well for very large variants. Small variants are over-penalized by RO alone since small
breakpoint differences cause larger movement in RO for small variants than large variants.

| Parameter | Type | Default | Range | Description |
| --- | --- | --- | --- | --- |
| ro | num | 0.5 | 0.0 < ro ≤ 1.0 | Proportion of overlap (defaults to 50% RO) |
| dist | int_u | unlimited | 0 ≤ dist ≤ unlimited | Breakpoints must be within dist bp |


#### Merge parameter: szro

Size-reciprocal-overlap (szro) is more flexible than ro. Like ro, it requires that variant sizes match within some
percentage by size, but it relaxes the restrictions on placement.

| Parameter | Type | Default | Range | Description |
| --- | --- | --- | --- | --- |
| szro | num | 0.5 | 0.0 < ro ≤ 1.0 | Proportion of size overlap (defaults to 50% size match) |
| dist | int_u | 200 | 0 ≤ dist ≤ unlimited | Breakpoints must be within dist bp |
| szdist | num_u | unlimited | 0 ≤ dist ≤ unlimited | Breakpoints must be offset ≤ szdist * size |

One of dist or szdist must not be unlimited.

#### Merge parameter: distance

Merge by distance. Can be restricted by szro and szdist, but it's better to use szro for those parameters. By default,
the szro parameter is unset and ignore (size is not requried to match). Use this when comparing against variants with
unknown or extremely inaccurate sizes.

| Parameter | Type | Default | Range | Description |
| --- | --- | --- | --- | --- |
| szro | num | NA | 0.0 < ro ≤ 1.0 | Proportion of size overlap (defaults to 50% size match) |
| dist | int_u | 500 | 0 ≤ dist ≤ unlimited | Breakpoints must be within dist bp |
| szdist | num_u | unlimited | 0 ≤ dist ≤ unlimited | Breakpoints must be offset ≤ szdist * size |

#### Merge parameter: exact

Requires exact position and size matches. There are no parameters to control exact matches. If match parameters are
also given, matches are further restricted by sequence content, which may or may not exactly match (even though the
size and position match exactly).

#### Match parameters

Sequences are compared by alignment similarity or Jaccard.

Alignment similarity is computed by taking the alignment score based on a Smith-Waterman alignment and dividing by
the maximum score if all bases match based on the match score and the length of the smaller variant
(match * min(len A, len B)). Alignment arguments (match, mismatch, open, and extend) parameterize an affine
Smith-Waterman alignment over the sequences (affine allows separate gap-open and gap-extend scores to better model true
biological indels and SVs). This method is far more accurate than Jaccard, especially for smaller variants.

For larger variants, comparison by alignment consumes excessive CPU time, so the match can fall-back to Jaccard
above a set size (limit). Each sequence is k-merized using a set k-mer size (ksize) and Jaccard similarity is
computed, Jaccard = A&B / (!AB + A&B + A!B). If the limit is set to unlimited, all variant are compared by alignment.

Whether computed by alignment or Jaccard, the similarity score must meet or exceed a set threshold (score).

| Parameter | Type | Default | Range | Description |
| --- | --- | --- | --- | --- |
| score | num | 0.8 | 0 < score ≤ 1.0 | Minimum sequence match proportion |
| match | num | 2.0 | 0 < match | Alignment base match score |
| mismatch | num | -1.0 | mismatch < 0.0 | Alignment base mismatch score |
| open | num | -1.0 | open ≤ 0.0 | Gap open score |
| extend | num | -0.25 | extend ≤ 0.0 | Gap extend score |
| limit | int_u | 4000 | 0 ≤ limit ≤ unlimited | Use Jaccard for length > limit |
| ksize | int | 9 | 0 < ksize | K-mer size for Jaccard |


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
                "svindel": "nr::szro(0.8,,4):match(0.8)",
                "sv:inv": "nr::ro(0.8)",
                "snv": "nrsnv::exact"
            },
            "name": "PAV HiFi",
            "description": "Pav HiFi"
        }

`sourcetype` and `sourcename` determine where samples are pulled from. The `merge` section can be a single merge
definition or a dictionary of merge strategies keyed by `svtype` (comma-separated list of svtypes is allowed). The
`name` will show up in figures generated for this merged callset. `description` is for documentation and is unused.

To merge SVs or indels, the entry vartype must be "svindel" since they are merged together.

The svtype (after the ":") is optional, but if it is present, it must match the svtype being merged. Entries with a
matching svtype are prioritized over entries without a matching svtype. If more than one entry can match, then the
first one in the list has the higher precedence. A key of "DEFAULT" matches everything with the lowest precedence.


### callerset

Example callerset section.

    "callerset": {
        "longreads": {
            "callsets": [
                ["caller", "pav-hifi", "PAVHIFI"],
                ["caller", "pav-clr", "PAVCLR"],
                ["caller", "pbsv-hifi", "PBSVHIFI"],
                ["caller", "pbsv-clr", "PBSVCLR"]
            ],
            "merge": {
                "sv:ins,del": "nr::szro(0.8,,4):match(0.8)",
                "sv:inv": "nr::ro(0.8)",
                "indel": "nr::szro(0.8,,4):match(0.8)",
                "snv": "nrsnv::exact"
            },
            "name": "HGSVC2 LR",
            "description": "PAV/PBSV (HiFi & CLR)"
        }
    }

`callsets` contains a list variant sources (each is a souncetype/souncename/name trio). In an older version of SV-Pop,
a separate `namelist` element contained the names of each source, but that information has been moved to the
`callsets` list as a third element for each source. `merge` is the same as for `sampleset` (a single merge definition
or a dictionary of definitions keyed by svtype(s)). `name` is used for figures, and `description` is for documentation
purposes and is unused by the pipeline.


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
        "szro-80": "nr::szro(0.8,,4):match(0.8)",
        "snv-exact": "nrsnv::exact",
        "ro-80-nm": "nr::ro(0.8)"
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
