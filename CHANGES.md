# Changes

## 2.2.0

* Added "match" option to nr strategy. Value is a comma-separated string where fields
  may be missing or empty to accept the default parameter (e.g. "match=0.75"). Fields:
  * SCORE_PROP [0.8]: Alignment score must be this proportion of the max where max
    is the score if all bases match (MATCH * max_length(A, B))
  * MATCH [2.0]: Alignment match score.
  * MISMATCH [-1.0]: Alignment mismatch score.
  * GAP-OPEN [-5]: Alignment gap-open score.
  * GAP-EXTEND [-0.5]: Alignment gap-extend score. 
  * MAP-LIMIT [500000]: Use Jaccard distance if the largest of the two sequences is
    larger that this size. Similarity is still SCORE_PROP.
  * JACCARD-KMER [9]: Jaccard k-mer size.

* Merge strategy "nrid" was removed. Replace with "nr:refalt" to do exact matches
  over position requiring REF and ALT to match.

* Merge strategy "fam" was removed. This was an NR strategy that added columns for
  inheritance. This annotation should be done outside SV-Pop.


# 2.3.3

Merging
* Dropped "expand" support for merging
* Dropped "MERGE_AC" and "MERGE_AF" columns. These are not true AC and AF calculations without confident
  genotypes, which depends heavily on the input callset. A future version might consider GT if present.
* Added "MERGE_N". Counts the number of samples supporting a variant.
