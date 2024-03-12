# SV-Pop targets

The SV-Pop pipeline is run by requesting one or more target files. Each file has a specific path with wildcards that
guide SV-Pop. Targets may be callsets (parsed from input or merged) or annotations on callsets. This guide lists
available targets.

## Variant calls

Base variant BED files (bed 6+).

Columns are:
1. #CHROM: Chromosome
1. POS: Position
1. END: End position
1. ID: Unique variant ID (only unique to a single file, other files may have another variant with the same ID)
1. SVTYPE: SV type (ins, del, inv, etc)
1. SVLEN: Length of the SV
1. REF: Reference base (SNV)
1. ALT: Alternate base (SVV)

### BED and variant FASTA

BED file:

    results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/bed/{vartype}_{svtype}.bed.gz

Sequence-resolved variants have a FASTA file with sequences. FASTA record IDs match BED variant IDs.

Variant sequence FASTA:

    results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/all/bed/fa/{vartype}_{svtype}.fa.gz

### Variant summary tables

    results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/tables/variant_summary/{vartype}_{svtype}.tsv.gz

Can also generate `.xlsx` after `.tsv.gz`.


## Intersect

Intersect variants:

    results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/intersect.tsv.gz

1. sourcetype_a, sourcename_a, sample_a: sourcetype, sourcename, and sample for the first variant set.
1. sourcetype_b, sourcename_b, sample_b: sourcetype, sourcename, and sample for the second variant set.
1. merge_def: Merge definition.
1. filter, svset, vartype, and svtype: Same wildcards for variants. 

A Venn of intersects can be generated (requires matplotlib_venn):

    results/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/venn/variant_venn.{ext}

ext may be "pdf" or "png". Both will be generated regardless of which is requested.

Nearest intersects (distance to nearest variant of the same type):

    results/variant/intersect_nearest/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{filter}/{svset}/{vartype}_{svtype}/intersect.tsv.gz

## Annotations

### Alt-map

Insertion sequences can be re-mapped. Currently the only `mapper` supported is `minimap2`. 

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/altmap/altmap-{mapper}_{vartype}_{svtype}.bed.gz

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/altmap/altmap-{mapper}{vartype}_{svtype}.bam

Find Alt-map hits near the original insertion sites (finds putative tandem duplications):

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/altmap/altmap-{mapper}-distance-{dist_prop}-{ro}-{map_prop}_sv_ins.bed.gz

1. dist_prop: Allowed distance from the insertion site as a proportion of the SV length (e.g. 2 is 2x SVLEN from POS).
1. ro: Size reciprocal overlap between SV and mapped site ("any" means do not restrict).
1. map_prop: Require at least this proportion of the SV to map (e.g. 90 means 90% of the SV insertion sequence was mapped).


### Homopolymer and dinucleotides

Homopolymer and dinucleotide repeats are a source of error and accelerated variation.

`rep_len` is "homopolymer" or "dinucleotide"

Direct intersect:

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/homopolymer/{rep_len}_intersect_all_any_{vartype}_{svtype}.tsv.gz

Find nearest non-intersecting site:

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/homopolymer/{rep_len}_nearest_all_{flank}_{vartype}_{svtype}.tsv.gz

`flank` is "up" (upstream in reference space) or "dn" (downstream in reference space).


### RefSeq

For each RefSeq hit, count the number of gene features affected (CDS, UTR5, UTR3, INTRON):

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/refseq/refseq-count_{vartype}_{svtype}.tsv.gz

Find all RefSeq genes proximal to the variant

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/refseq-prox/refseq-prox-{direction}-{flank}_{vartype}_{svtype}.tsv.gz


1. direction: "up" (upstream of gene relative to gene orientation) or "dn" (downstream of the gene relative to gene orientation)
1. flank: Max distance in basepairs away from gene.


### Reference annotated regions

Intersect with reference annotations pulled from the UCSC genome browser. Many of the regions are redundant and are
merged. For example, segmental duplications (SDs) and tandem repeats (TRs) have many overlapping records. These
annotations merge those into a flat region.

1. distance: Merge records if they are within this many bases.
1. flank: Add this number of bases to each merged record (expands regions).
1. overlap: Percentage of variant that must intersect with a region (e.g. 50 is 50% of the variant is in a region).


    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/{annotype}/{annoname}_regions_{distance}_{flank}_{overlap}_{vartype}_{svtype}.tsv.gz


| annotype | annoname | Description |
| --- | --- | --- |
| cen | cen | Centromeres |
| gap | gap | Reference gaps |
| agp | agp | Reference source contig |
| agp | agp_switch_{flank} | Within `flank` bp of a reference contig switch |
| bands | bands | Chromosome bands |
| refseq | refseq | Whole refseq region (refseq-count is recommended) |
| trf | trf | Simple repeats (STR/VNTR) |
| sd | sd-max-{match_type} | Get max SD identity for each SD intersect |
| sd | sd | SD intersect (True/False) |
| dhs_cluster | dhs_cluster_{score} | ENCODE DHS clusters with min score (or "all") |
| dhs2020 | dhs_{score} | Vierstra 2020 DHS sites with min score (or "all") |
| ccre | ccre2020 | Candidate cis regulatory elements |
| oreganno | oreganno | Curated regulatory sites |

Notes:
1. agp_switch_{flank}: "flank" in "agp_switch_{flank}" is not the same as the "flank" that comes after regions in the
   path ("..._{annoname}_regions_{distance}_{flank}_{overlap}_{vartype}_{svtype}.tsv.gz").
1. refseq: refseq-count (see above) is a better annotation. The refseq here is just the whole region, but would not
   annotate which genes were intersected or how.
1. sd: Does not identify the SD identity. sd-max is a better annotation (see below).

SD: Segmental duplication<br/>
DHS: DNAase hypersensitivity

### Chromosome bands

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/bands/bands_{vartype}_{svtype}.tsv.gz


### Max SD

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/sd/sd-max-{match_type}_{vartype}_{svtype}.tsv.gz

match_type:

1. frac: Identity not including indels.
1. fracindel: Identity including indels.

frac and fracindel are very close, we tend to use frac.


### RepeatMasker

Intersect with reference RepeatMasker annotations.

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/rmsk/rmsk-{filter_spec}-{ident}_intersect_{vartype}_{svtype}.tsv.gz

### ORegAnno

Intersect with curated regulatory elements, ORegAnno.

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/oreganno/oreganno_{vartype}_{svtype}.tsv.gz

Tags variants with ORegAnno ID. A separate table of annotations for each ID, such as the candidate regulatory element
and the gene it likely affects, is available from ORegAnno.

### ENCODE marks


    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/encode/encode-{mark}-{threshold}-{overlap}_{vartype}_{svtype}.tsv.gz

1. mark: H3K27Ac H3K4Me3 H3K4Me1
1. threshold: Min annotation score
1. overlap: Percentage of variant that must intersect with a region (e.g. 50 is 50% of the variant is in a region).


### cCRE

Candidate cis-regulatory elements (ENCODE 2020).

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/ccre/ccre2020-{score}_{vartype}_{svtype}.tsv.gz

1. score: Minimum annotation score.


### DNAse 2020

ENCODE DHS marks updated in 2020.

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/dhs/dhs2020-{score}_{vartype}_{svtype}.tsv.gz

1. score: Minimum annotation score.


### RepeatMasker and TRF on SV sequences

Derived from reference annotations on UCSC.

RepeatMasker:

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/rmsk/rmsk-table_{vartype}_{svtype}.tsv.gz
    
TRF (Tandem Repeats Finder)

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/trf/trf-table_{vartype}_{svtype}.tsv.gz


### Sequence content

Annotations run on variant sequences. Requires sequence-resolved variants.

GC content:

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/gc/gc_content_{vartype}_{svtype}.tsv.gz

Breakpoint homology with the reference. Many variants are driven by homology (recombination, replication based repair)
or generate homology (target site duplications, TSDs, for mobile elements). This annotation looks for perfect homology
between the SV sequence and the reference. It relies on exact breakpoint placement and stops at the first mismatch.
PAV generates this homology annotation, but SV-Pop provides this for variant called by other sources.
    
    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/seq/break_hom_ref_{vartype}_{svtype}.tsv.gz

Table of variant IDs and a SEQ column with the variant sequence. Used to pull sequences into a table if needed. Requires
sequence-resolved variant input.

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/seq/seq_{vartype}_{svtype}.tsv.gz

## Flag rules

Flag rules use lists of sample names (`samplelist` section of `config.json`) to run a collection of targets for multiple
samples. Unlike other targets, the output file of a flag rule is an empty flag file once all the targets are run, and
this flag file can be immediately removed if desired. Running a flag rule will check all the target files and re-run
them if any are missing or data they used has changed. These flag targets are a convenience, all the targets they run
can be executed without them, but the command-line process is greatly simplified when there are many samples to be run.

Two types of flag targets exist:
1) Callset flag rules: Retrieve and prepare callset variants and annotations across samples.
1) Merge flag rules: Run an intersect across multiple samples. May run sample-to-sample or sample combinations.

### Callset flag rules

Variant calls, annotations, and variant sequence FASTA files:

```
flag/variant/caller/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz
flag/variant/caller/{sourcename}/{sample}/{filter}/{svset}/anno/{annodir}/{annotype}_{vartype}_{svtype}.{ext}.gz
flag/variant/caller/{sourcename}/{sample}/{filter}/{svset}/bed/fa/{vartype}_{svtype}.fa.gz
```

The sample wildcard is the name of a sample list. The rest of the wildcards are consistent with the per-sample variant
callest and annotation wildcards.

### Intersect flag rules

```
flag/variant/intersect/{sourcetype_a}+{sourcename_a}+{sample_a}/{sourcetype_b}+{sourcename_b}+{sample_b}/{merge_def}/{filter}/{svset}/{vartype}_{svtype}/{policy}
```

The two samples, `sample_a` and `sample_b`, can be defined sample lists, but may be individual sample names (e.g.
compare many samples to one).

The `policy` wildcard has three options:
1) match: One-to-one match for the `sample_a` list and the `sample_b` list. For example, if one list contains "a", "b",
  and "c", and the other list list is "x", "y", "z", then "a" is compared to "x", "b" is compared to "y", and "c" is
  compared to "z".
1) comb: Run all combinations of the sample lists except where the sample names match.
1) fullcomb: Run all combinations of the sample lists including where the sample names match.

The sample list may have ":hap" or ":hapmatch" appended to it to accomodate callsets with and without haplotype
assignments.

When ":hap" is appended to a sample list, it expands the list by adding "_h1" and "_h2" to each sample
name. For example, if list "s1" is ["a", "b"], then setting a sample wildcard to "s1:hap" expands it to ["a_h1", "a_h2",
"b_h1", "b_h2"]. Note the order of the original list is preserved (i.e. "a" still comes before "b"), so the list with
haplotypes is predictable. If both `sample_a` and `sample_b` are "s1:hap", then "a_h1" is compared with "a_h1",
"b_h1" is compared with "b_h1", and so on. This mode accomodates PAV callsets imported as "pavbedhap" in the sample
table.

When ":hapmatch" is appended, the list is expanded, but the haplotype assignments are not made. For example, if list
"s1" is ["a", "b"], then setting a sample wildcard to "s1:hapmatch" expands it to ["a", "a", "b", "b"]. This is
meant to accommodate cases where phased callsets are compared with unphased callsets (e.g. PAV imported as "pavbedhap"
with PBSV calls for the same sample).
