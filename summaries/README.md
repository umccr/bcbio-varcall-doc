Summarising Somatic Variant Caller Filters
==========================================

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [VarDict-java Filters](#vardict-java-filters)
    * [`var2vcf_paired`](#var2vcf_paired)
    * [`bcftools`](#bcftools)
    * [`depth-freq-filter`](#depth-freq-filter)
    * [`freebayes-call-somatic`](#freebayes-call-somatic)
* [FreeBayes Filters](#freebayes-filters)
    * [`freebayes-call-somatic`](#freebayes-call-somatic-1)
* [MuTect2 Filters](#mutect2-filters)

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
* Min allele frequency (`-f 0.1`) (also in `VarDict-java`)
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

FreeBayes Filters
-----------------

Only one default filter changed in bcbio command line:

* Min alternate fraction (`--min-alternate-fraction 0.1`)
    - require at least 10% of observations supporting an alternate allele
      within a single individual in order to evaluate the position.

Next are default (several others default to 0 - can see with `freebayes -h`):

* Min alternate count (`--min-alternate-count 2`)
    - require at least 2 observations supporting an alternate allele
      within a single individual in order to evaluate the position.
* Read duplicate: exclude duplicates marked as such in alignments
* Min Read Mapping Quality (`-m 1`)
    - exclude reads from analysis if the have mapping quality < m
* Min Base Quality (`-q 0`)
    - exclude alleles from analysis if their supporting base quality
      is < q

### `freebayes-call-somatic`

Unsure about this:

Function in [bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/variation/freebayes.py)

* Call SOMATIC variants from tumor/normal calls, adding REJECT filters and SOMATIC flag.
* Works from stdin and writes to stdout, finding positions of tumor and normal samples.
* Extracts the genotype likelihoods (GLs) from FreeBayes, which are like phred scores
  except not multiplied by 10.
    * for tumors, we retrieve the best likelihood to not be reference (the first GL)
    * for normal, the best likelhood to be reference.
* After calculating the likelihoods, we compare these to thresholds to pass variants
  at tuned sensitivity/precision. Tuning done on DREAM synthetic 3 dataset evaluations.
* We also check that the frequency of the tumor exceeds the frequency of the normal by
  a threshold to avoid calls that are low frequency in both tumor and normal. This supports
  both FreeBayes and VarDict output frequencies.

bcbio has several steps after this, mostly cleaning up fields and normalising/removing dup variants with `vt`.


MuTect2 Filters
---------------

* `tumor_lod`: minimum likelihood of an allele as determined by the somatic
  likelihoods model required to pass. `TUMOR_LOD_THRESHOLD = 5.3`.
    * LOD threshold for calling tumor variant. Only variants with tumor LODs
      exceeding the threshold can pass filtering
* `maxEventsInHaplotype`: the maximum allowable number of called variants
  co-occurring in a single assembly region. If the number of called variants
  exceeds this they will all be filtered. Note that this filter is misnamed
  because it counts the total number of events over all haplotypes in an
  assembly region. `maxEventsInHaplotype = 2`.
    * Variants coming from a haplotype with more than this many events are filtered
* `uniqueAltReadCount`: the minimum number of unique (start position, fragment
  length) pairs required to make a call. This count is a proxy for the number
  of unique molecules (as opposed to PCR duplicates) supporting an allele.
  Normally PCR duplicates are marked and filtered by the GATK engine, but in
  UMI-aware calling this may not be the case, hence the need for this filter.
  `uniqueAltReadCount = 0`.
    * Filter a variant if a site contains fewer than this many unique
      (i.e. deduplicated) reads supporting the alternate allele
* `maxAltAllelesThreshold`: the maximum allowable number of alt alleles at a
  site. By default only biallelic variants pass the filter. `numAltAllelesThreshold = 1`
    * Filter variants with too many alt alleles.
* `max_germline_posterior`: the maximum posterior probability, as determined by
  the above germline probability model, that a variant is a germline event.
  `maxGermlinePosterior = 0.025`.
    * Maximum posterior probability that an allele is a germline variant.
      The default value of this argument was chosen to achieve roughly one false
      positive per covered megabase according to a back-of-the-envelope
      calculation assuming the neutral model of allele selection, with
      parameters estimated from ExAC.
* `normal_artifact_lod`: the maximum acceptable likelihood of an allele in the
  normal _by the somatic likelihoods model_. This is different from the normal
  likelihood that goes into the germline model, which makes a diploid assumption.
  Here we compute the normal likelihood as if it were a tumor in order to detect
  artifacts. `NORMAL_ARTIFACT_LOD_THRESHOLD = 0.0`.
    * LOD threshold for calling normal artifacts. Measure of the minimum
      evidence to support that a variant observed in the tumor is not
      also present in the normal as an artifact i.e. not as a germline event.
* `strandArtifactPosteriorProbability`: the posterior probability of a strand
  artifact, as determined by the model described above, required to apply the
  strand artifact filter. This is necessary but not sufficient - we also
  require the estimated max a posteriori allele fraction to be less than
  `strandArtifactAlleleFraction`. The second condition prevents filtering
  real variants that also have significant strand bias, i.e. a true
  variant that _also_ has some artifactual reads. `STRAND_ARTIFACT_POSTERIOR_PROB_THRESHOLD = 0.99`.
    * Filter a variant if the probability of strand artifact exceeds this
      number.
* `strandArtifactAlleleFraction`: `STRAND_ARTIFACT_ALLELE_FRACTION_THRESHOLD = 0.01`.
    * Filter a variant if the MAP estimate of allele fraction given artifact is
      below this number.
* `minMedianBaseQuality`: the minimum median base quality of bases supporting a SNV.
  `minMedianBaseQuality = 20`.
    * Filter variants for which the median base quality of alt reads is too
      low.
* `minMedianMappingQuality`: the minimum median mapping quality of reads supporting an allele.
  `minMedianMappingQuality = 30`.
    * Filter variants for which the median mapping quality of alt reads is too
      low.
* `maxMedianFragmentLengthDifference`: the maximum difference between the median fragment
  lengths reads supporting alt and reference alleles. Note that fragment
  length is based on where paired reads are mapped, not the actual physical fragment length.
  `maxMedianFragmentLengthDifference = 10000`.
    * Filter variants for which the median fragment length for alt reads is very
      different from the median for ref reads.
* `minMedianReadPosition`: the minimum median length of bases supporting an
  allele from the closest end of the read. Indels positions are measured by
  the end farthest from the end of the read. `minMedianReadPosition = 5`.
    * Filter variants for which the median position of alt alleles within reads
      is too near the end of reads.
* `contaminationTable`: if `FilterMutectCalls` is passed a `contaminationTable`
  from `CalculateContamination` it will filter alleles with allele fraction less
  than the whole-bam contamination in the table. `contaminationTable = null`.
    * Table containing contamination information.
