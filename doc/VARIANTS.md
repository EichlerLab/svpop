# Variant Rules

Rules for extracting variants from callsets and merging them are in `rules/variant`. These may be variants
directly from the callers (`rules/variant/caller`), variants merged from multiple callers per sample
(`rules/variant/callerset`), or variants merged from multiple samples (`rules/variant/mergeset`). Pipelines
will likely combine these, for example, merge variants (from callers) per sample (callerset), and then merge
a family or create a nonredundant set of variant from mulitple samples from those (mergeset).

## Callers

