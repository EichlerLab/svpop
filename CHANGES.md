# Changes

## 3.4.1
* Added dvpepper (PEPPER-Margin-DeepVariant) and SVIM input VCF parsers.
* Added min_svlen and max_svlen options for VCF input parsers.

## 3.4.0
* Merging and intersects bring SVs and indels together for the intersect and separate back out after.
  * Avoids mis-merging at SV/indel boundaries.
* Added flag rules to run many related operations at once (for each sample in a sample list).
* Subset chromosome option supports multiple chromosomes.

## 3.3.7
* Added "pavbedhap" variant source for PAV input from each haplotype (not merged at the sample level).
* Added filter pass option for parsing VCFs and controlling which variants are accepted based on the FILTER column.
* Delly input VCF parser.
* Fixed bug with upstream deletions (ALT=* was not ignored)

## 3.3.6
* SNV format changed to "CHROM-POS-SVTYPE-REFALT" (eliminated dash between REF and ALT). All variant IDs can be
  split on "-" to obtain 4 fields (note that a ".n" may also be appended to distinguish multiple calls of the same
  type and at the same location).
* Improved multi-allele VCF handling. Resulting tables now include "VCF_ALT_IDX" to match the alt genotypes each
  record was retrieved from. Each alternate allele will be separated into a separate record and assigned to samples
  where VCF_ALT_IDX is in GT.

## 3.3.5
* Added threads to rules (supports --cores)
* Fixed output format in rule hpref_merge_bed (was not gzipped)
* Moved VCF parsing code to svpoplib (supports svpoplib use as an external library)
* Input VCF parser assumes GT "1/." if the GT field is missing
* Removed GT field for CuteSV (writes "./." causing the parser to drop variants)
* Better support for writing VCF files with no variants (was causing crashes)
* Can input dataframes to merging/intersect svpoplib routines (supports svpoplib use as an external library)

## 3.3.4
* Fixed parsing Sniffles2 VCFs that are incorrectly formatted with "N" and SEQ in REF/ALT instead of symbolic ALTs 

## 3.3.3
* Bug in ID de-duplication code
* Sorting RMSK annotations

## 3.3.2
* Fixed bugs processing chromosome names containing a "." when de-duplicating variant IDs

## 3.3.1
* Removed verbose output from merging (now optional)

## 3.3.0
* Fixed pavlib functions PAV uses that cause numeric chromosome problems.

## 3.2.0
* Set ply submodule version to 3.10
* Alt-map retains CIGAR operations
* Alt-DUP calls inner-variants (SV/indel INS/DEL & SNVs) from duplications using the alt-map CIGAR

## 3.1.0
* CCRE annotations needed to be sorted.
* VCF input support for DeepVariant and multiallelic sites.
* Moved Sniffles and SVIM-asm input parsers to the VCF parser framework (no longer custom parsers, not needed).
* Variant FASTA files get FAI.

## 3.0.0
* Added altdup for remapping INS as DUPs (allows DUP version to be treated as a separate callset - i.e. merging, intersects, annotations)
* Merging handles multi-allelic sites better.
* Flexible merging parameter backend.
* Revived VCF writer for standard BED files.
* VCF: SVLEN header had the incorrect data type. 

## 2.3.3
* Dropped "expand" support for merging
* Dropped "MERGE_AC" and "MERGE_AF" columns. These are not true AC and AF calculations without confident
  genotypes, which depends heavily on the input callset. A future version might consider GT if present.
* Added "MERGE_N". Counts the number of samples supporting a variant.

## 2.3.2
* Changed default merging parameters

## 2.3.1
* Updated pipeline documentation

## 2.3.0
* Added kanapy as a submodule

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

