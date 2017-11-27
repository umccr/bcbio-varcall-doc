Summarising Somatic Variant Caller Filters
==========================================

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [VarDict-java Filters](#vardict-java-filters)
    * [`var2vcf_paired`](#var2vcf_paired)
    * [`bcftools`](#bcftools)
    * [`depth-freq-filter`](#depth-freq-filter)
    * [`freebayes-call-somatic`](#freebayes-call-somatic)

<!-- vim-markdown-toc -->

Introduction
------------
What sort of stuff can we look at when filtering (somatic) variants?

* Base Quality Score decreases towards the end of the read; error likelier than variant
    - Min BQ threshold for variant alleles
    - Avg Read Position of variant

* Duplicate reads from PCR amplification
    - Mark duplicates post-alignment

* Strand bias: variant supported by reads aligning to mostly one strand might be artefact
    - Discard variants within reads aligning in one direction only

* GC bias: lower coverage in low-GC content regions. Germline SNVs might be mistaken for somatic variants
  due to the lower coverage.
    - Require minimum read depth at variant locus especially in normal sample

* Misalignment due to missing reference genome sequence
    - Local indel realignment post-alignment
    - Discard variants within/close to germline indels
    - Discard variants positioned always toward start/end of alignment

* Misalignment due to repetitive regions: MQ0 for multi-aligned reads; low-mappability regions out of
  bounds for short-read sequencing
    - Minimum Mapping Quality MQ
    - Discard variants in low-complexity/low-mappability regions

VarDict-java Filters
--------------------

* Read misalignment/duplicate (`F 0x700`)
    - not primary alignment
    - read fails platform/vendor qc
    - read is PCR or optical duplicate
* Min allele frequency threshold (`-f 0.1`) (also in `var2vcf_paired`)
* Min allele frequency in a normal sample allowed for a putative somatic mutation (`-V 0.05`)
* Min Read Mapping Quality (`-Q 10`)
* Min Base Quality (`-q 25`)
* Local realignment (`-k 1`)
* Indel size in bp (`-I 120`)
* Min num of reads supporting variant (`-r 2`)
* Min num of reads to determine Strand Bias (`-B 2`)
* Max num mismatches allowed in read (`-m 8`)
* Min read position - if the mean read position of the variant is less than this, it's filtered out (`-P 5`)

### `var2vcf_paired`

* Max p-value for somatic evidence (`-P 0.9`)
* Max mean mismatches allowed (`-m 4.25`)
* Min allele frequency (`-f 0.1`)
* Output only candidate somatic variants (`-M`)
* Min total depth (`-d 5`)
* Min variant depth (`-v 3`)
* Min Mapping Quality (`-Q 0`)
* Min mean position of variant in reads (`-p 5`)
* Min Base Quality (`-q 22.5`)

### `bcftools`
* Add 'REJECT' soft filter to variants where the INFO-STATUS field doesn't contain the string 'Somatic'.
* Filter out variants with `QUAL` < 0

### `depth-freq-filter`
Command line to filter VarDict calls based on depth, frequency and quality.
Looks at regions with low depth for allele frequency (AF * DP < 6, the equivalent
of < 13bp for heterozygote calls, but generalised.

Within these calls filters if a call has:

- Low mapping quality and multiple mismatches in a read (NM)
    - For bwa only: MQ < 55.0 and NM > 1.0 or MQ < 60.0 and NM > 2.0
- Low depth (DP < 10)
- Low QUAL (QUAL < 45)

Also filters in low allele frequency regions with poor quality, if all of these are
true:

- Allele frequency < 0.2
- Quality < 55
- P-value (SSF) > 0.06

### `freebayes-call-somatic`
* For VarDict I believe it simply checks the AF of the tumor compared to the
  normal based on a threshold. If there is supporting evidence, adds 'SOMATIC'
  to the INFO column. Else, adds 'REJECT' to the FILTER column.

