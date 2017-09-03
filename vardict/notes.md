Publication
-----------
* Simultaneously calls SNV, MNV, indels, complex and structural varints
* Performs local realignments on the fly for more accurate allele freq
  estimation (e.g. for sub-clones, circulating tumor DNA)
* Can handle ultra-deep sequenced samples for detection of low AF mutations
* Detects differences in somatic and LOH variants between paired samples 

### Local realignment and indel calling
* Supervised and unsupervised local realignment
    * **supervised**: identify mismatched
      alignment of 3' and 5' read ends flanking short indels
      and add in support of indel, increasing allele frequency.
    * **Unsupervised**: scan local sequences near soft-clippings
      to look for larger indels.

### Complex variants
* Combinations of insertions and deletions
* Represents them as a single variant

### Structural variants
* Two-step approach:
    * Build a consensus sequence from clipped sequences 
    * Search whether this consensus can be uniquely aligned withing 5kb from
      of the given region

### Tumor-Normal analysis
* Extract read counts for REF and ALT alleles
* Perform Fisher's exact test to see if variant AF is significantly
  different between Tumor and Normal. Variant is classified as:
    * Somatic: present in tumor sample
    * Germline: present in both samples
    * LOH: heterozygous in normal but homozygous or lost in tumor
    * Deleted: present in normal but no coverage in tumor

### De-duplication
* Built-in option to perform de-duplication.
* Read pairs with same alignment positions for both reads are
  deemed duplicates
    * Only first one encountered is used for varcalling
