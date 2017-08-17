bcbio blog post notes
=====================

Ensemble method (2013-02-06)
----------------------------

* [Blog post](<https://bcbio.wordpress.com/2013/02/06/an-automated-ensemble-method-for-combining-and-evaluating-genomic-variants-from-multiple-callers/>)

### Problem
Say we use tool1 and tool2 to call variants in a sample. We can
check the number of variants that were called in both tools (concordant), those
that were called in only tool1 (discordant1) and only tool2 (discordant2).
How can we treat these sets?
* If we only focus on the concordant variants we might be missing out on
  true variants in the discordant groups, getting a lot of false negatives.
* If we focus on the concordant and discordant variants we will get a lot of
  false positives.

### Method

#### Pipeline
* Get two NA12878 exome replicates from EdgeBio. These were independently
  prepared using a Nimblegen capture kit and sequenced on an Illumina HiSeq.
* Align reads with Novoalign, deduplicate, base recalibrate and realign with
  GATK.
* Call variants with five tools.

#### Combine and annotate
* Get **union** of all variant calls
* Annotate each variant with metrics (e.g. strand bias, allele balance, regional
  sequence entropy, position of calls within reads, regional base quality,
  overall genotype likelihoods)

#### Filter
* First keep variants supported by `>=` N callers ('trusted variants')
* For variants supported by `<` N callers, use an SVM to distinguish between
  true and false positives:
    * The annotated metrics above are the input parameters
    * True positives are variants found in all callers
    * False positives are those found in a single caller
    * Use above training variants to get an initial set of below-cutoff variants
      to include and exclude, then re-train multiple classifiers stratified
      based on variant characteristics: variant type (indels vs. SNPs), zygosity
      (hom vs. het) and regional sequence complexity.
    * Use these final classifiers to identify included and excluded variants
      falling below the trusted calling support cutoff.
* Final set of variants includes the trusted ones and those that pass the SVM
  filtering. So at this stage we have two final sets of variants - one for each
  NA12878 replicate.

### Results
Compare concordant and discordant variant calls between the two replicates. Use
tools 1-5 as a baseline, then compare with ensemble approach.
