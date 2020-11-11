# SV-Pop Pipeline

This pipeline performs many tasks related to variant parsing, merging, and annotation. It was initially developed for
"Characterizing the Major Structural Variant Alleles of the Human Genome" (Audano, Sulovari, et al. 2019. Cell), but
has since evolved significantly.

The pipeline takes variants from different callers. The pipeline has built-in parsers for some callers, such as PBSV and
Sniffles, and so variants can be read directly from their output VCFs. Variants from any other source can
be read from prepared BED files.

This is not the same as "SV-Pop: population-based structural variant analysis and visualization"
(Ravenhall et al. 2019. BMC Bioinformatics). It was named before that paper came out, and it will probably get a new
name eventually.


## Configuration

The main configuration file `config/config.json`, and is in JSON format. It defines paths to the reference, input,
rules for merging, and links to external variant sources.

For variants with built-in parsers (PBSV, Sniffles, SMRT-SV, and Phased-SV), `config/samples.tab` defines where to find these
samples.

Some genotyping analysis steps need to know sex and ethnicity of samples, and that is stored in `config/sample_info.tab`,
although this files is not needed for most analyses. 


## Running the pipeline

The piplenie should be run from a clean analysis directory specific for the project. Pipeline code is installed in
another directory.

There is no defined endpoint for this pipeline. Desired output files are requested by running a Snakemake command,
and the output is generated. That output may be a set of variant calls, merged variants for multiple samples and/or
callers, annotations, or plots.

Example: `results/variant/caller/pbsv-hifi/bed/HG00733/all/lc/byref/sv_ins.bed`

This file requests variants from caller PBSV (from a HiFi run) in BED format for HG00733. It applies no "variant set"
filters (e.g. TRF or SD; the `all` directory), does apply a low-confidence filter around centromeres and
peri-centromeres (`lc`), and requests BED files in reference context (`byref`). The variants extracted for this set are
SV insertions (`sv_ins.bed`).

Most output will be in `results/variant`. The next two directories are the "sourcetype" and "sourcename". If there is no
merging (see below), then "sourcetype" is "caller". If variants are merged, then "sourcetype" will be "sampleset" or
"callerset", and "sourcename" will be the name of a merge configuration parameter set specified in `config/config.json`.
See below for more information about merging.

## Annotations

The pipeline does many annotations including intersecting with UCSC (see below), running TRF and RepeatMasker on
inserted or deleted SV sequence, reference mapping location for SV insertions (finding SV donor sites for duplications),
and homopolymer run intersections. 

* UCSC tracks: GRC patches, centromeres, gaps, AGP, chromosome band, RefSeq (CDS, ncRNA, intron, upstream/downstream
flank), TRF, segmental duplications (SD), RepeatMasker, CpG islands, ENCODE histone marks, and ORegAnno.

The pipeline runs these annotations only original variant callsets before subsetting or merging. This way, variants can
annotated once, then merged and subset in any number of ways without re-annotating.

Example: `results/variant/caller/extern-asmhifi/anno/HG00733/all/lc/refseq/refseq-count_sv_ins.tab`

This file requests RefSeq intersections for PBSV HiFi HG00733 SV insertions. The table counts the number of bases
affected within coding regions, UTRs, introns, and ncRNA exons.


## Merging variants

The pipeline can merge variants in two ways:
  1. Callerset: Merge variants from different callers for the same sample.
  1. Sampleset: Merge variants from different samples.

### Merge process overview

Merging is done by specifying an order for callers (callerset) or samples (sampleset), and merging is done iteratively
in that order. The first variant set is copied into the merged set. The second variant set is intersected with the
merged set, matching calls are annotated with the intersected variants, and new calls are added to the merged set. This
process continues iteratively until all callers or samples are merged into the merged set.

The result is that each variant in the merged set is represented by the call from the caller or sample where it was
first seen. Each merged variant contains enough information to trace back to the calls that support it; the list of
supporting callers/samples and variant IDs matching each are retained.

Example: `results/variant/sampleset/pansv/bed/all_samples/all/lc/byref/sv_ins.bed`

In this example, the configuration JSON will contain a definition for sampleset "pansv", which includes variant merging
parameters. A set of samples is defined in a list called "all_samples" (also defined in the configuration JSON) that
defines which samples are to be merged. Notice that a merging process can be configured once ("pansv") and used for
any number ofsample sets.

### Merging annotations

Because the information is saved, it allows annotations done on individual callers and samples to be pulled through
the annotation process.

Example: `results/variant/sampleset/pansv/anno/all_samples/all/lc/refseq/refseq-count_sv_ins.tab`

In this example, it takes RefSeq annotations for all samples defined as part of the "all_samples" group and merges them
into a new table with identical format containing annotations for the whole merged set. The result is the same as taking
the merged set and annotating it except it can be re-merged and subset without re-running the annotations.


## Comparing variants

Example: `results/variant/intersect/caller+pbsv-hifi+HG00733_vs_caller+pbsv-clr+HG00733/all/lc/roor-50-200/sv_ins/venn/variant_venn.pdf`

This compares two callsets, a PBSV HiFi set (`caller+pbsv-hifi+HG00733`) against a PBSV CLR set
(`caller+pbsv-clr+HG00733`) for sample HG00733. The intersect definition, `roor-50-200` is defined in the configuration
JSON.

The format of the directory name specifying samples is flexible enough to do many types of comparisons. Comparisons can
be done within samples, across samples, against merged variants, etc.

Format: `{sourcetype_a}+{sourcename_a}+{sample_a}_vs_{sourcetype_b}+{sourcename_b}+{sample_b}`

If the variant set field contains "_vs_" (is `all` in the above eample), then comparisons can be done with differently-filtered
variants. For example, `notr_vs_all` would compare variants outside tandem repeats (`notr`) with all variants `all`.

