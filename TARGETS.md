# SV-Pop targets

The SV-Pop pipeline is run by requesting one or more target files. Each file has a specific path with wildcards that
guide SV-Pop. Targets may be callsets (parsed from input or merged) or annotations on callsets. This guide lists
available targets.

## Variant calls

Base variant BED files (bed 6+).

Columns are:
1. \#CHROM: Chromosome
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
| trf | trf | Simple repeats (STR/VNTR) |
| sd | sd | Any SD region |
| sd | sd-max-{match_type} | Maximum SD identity |
| rmsk | rmsk-{filter_spec}-{ident} | RepeatMasker |
| oreganno | oreganno | Curated regulatory sites |
| dhs2020 | dhs_{score} | Vierstra 2020 DHS sites |
| dhs_cluster | dhs_cluster_{score} | Original ENCODE DHS clusters |
| ccre | ccre2020 | Candidate cis regulatory elements |
| encode | encode-{mark}-{threshold} | ENCODE mark |
| cpgisland | cpgisland_unmasked | CpG island including RepeatMasked regions |
| bands | bands | Chromosome bands |
| agp | agp | Reference scaffolds |
| gap | gap | Reference gaps |
| cen | cen | Centromeres |

SD: Segmental duplication<br/>
DHS: DNAase hypersensitivity

#### dhs2020

URL: 
https://resources.altius.org/~jvierstra/projects/footprinting.2020/consensus.index/consensus_footprints_and_motifs_hg38.bed.gz

    Vierstra, J., Lazar, J., Sandstrom, R. et al. Global reference mapping of human transcription factor footprints. Nature 583, 729–736 (2020).

#### sd-max

match_type:

1. frac: Identity not including indels.
1. fracindel: Identity including indels.

frac and fracindel are very close, we tend to use frac.

#### dhs2020

1. score: Minimum score (mean DHS signal) or "all".


### RepeatMasker

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/rmsk/rmsk-{filter_spec}-{ident}_intersect_{vartype}_{svtype}.tsv.gz


### ENCODE marks

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/encode/encode-{mark}-{threshold}-{overlap}_{vartype}_{svtype}.tsv.gz
    
    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/ccre/ccre2020-{score}_{vartype}_{svtype}.tsv.gz
    
    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/dhs/dhs2020-{score}_{vartype}_{svtype}.tsv.gz
    
    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/dhs/dhs-{score}_{vartype}_{svtype}.tsv.gz
    
    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/cpgisland/cpgisland-unmasked_{vartype}_{svtype}.tsv.gz


1. mark: H3K27Ac H3K4Me3 H3K4Me1
1. threshold: Min score
1. overlap: Percentage of variant that must intersect with a region (e.g. 50 is 50% of the variant is in a region).

Note: cpgisland-unmasked include RepeatMasked sites.


### RepeatMasker and TRF on SV sequences

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/rmsk/rmsk-table_{vartype}_{svtype,ins|del|inv|dup}.tsv.gz
    
    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/trf/trf-table_{vartype}_{svtype}.tsv.gz


### Sequence content

    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/gc/gc_content_{vartype}_{svtype}.tsv.gz
    
    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/seq/seq_{vartype}_{svtype}.tsv.gz
    
    results/variant/caller/{sourcename}/{sample}/{filter}/all/anno/seq/break_hom_ref_{vartype}_{svtype}.tsv.gz

