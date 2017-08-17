bcbio blog post notes
=====================


<!-- vim-markdown-toc GFM -->
* [Ensemble method (2013-02-06)](#ensemble-method-2013-02-06)
    * [Problem](#problem)
    * [Method](#method)
        * [Pipeline](#pipeline)
        * [Combine and annotate](#combine-and-annotate)
        * [Filter](#filter)
    * [Results](#results)
* [Evaluation framework (2013-05-06)](#evaluation-framework-2013-05-06)
    * [Problem](#problem-1)
    * [Method](#method-1)
        * [Aligners](#aligners)
        * [Post-alignment process](#post-alignment-process)
        * [Varcallers](#varcallers)
* [Updated Evaluation framework (2013-05-06)](#updated-evaluation-framework-2013-05-06)
* [WGS trio varcall evaluation (2014-05-12)](#wgs-trio-varcall-evaluation-2014-05-12)
    * [Method](#method-2)
    * [Filters](#filters)
        * [Low complexity regions](#low-complexity-regions)
        * [High depth and low quality](#high-depth-and-low-quality)
        * [GATK](#gatk)
* [Cancer variant callers (2015-03-05)](#cancer-variant-callers-2015-03-05)
    * [Problem](#problem-2)
    * [Variant Callers](#variant-callers)
    * [Tumor-only prioritisation](#tumor-only-prioritisation)

<!-- vim-markdown-toc -->

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

Evaluation framework (2013-05-06)
---------------------------------
* [Blog post](<https://bcbio.wordpress.com/2013/05/06/framework-for-evaluating-variant-detection-methods-comparison-of-aligners-and-callers/>)

### Problem
Relatively low concordance between variant calling methods. How can we evaluate
variant calls?

### Method
Compare:
* 2 aligners (bwa mem, novoalign)
* 2 post-alignment prep methods (GATK, gkno)
* 3 varcallers (GATK ug, GATK hc, FreeBayes)

* Call variants on NA12878 exome from EdgeBio (evaluation variants)
* Assess against NA12878 from GIAB (reference)

* Discordant variants where reference and evaluation variants differ are:
    * Extras: False positives (called in eval data but not in ref)
    * Missing: False negatives (found in ref but not called in eval data)
    * Shared: shared variants (found in ref and eval but represented
      differently e.g. allele differences (hom/het), indel start/end diff etc.)

* Take false positives and false negatives and subdivide using annotations from
  GEMINI (low coverage (4-9 reads), repetitive (RepeatMasker), error prone
  (motifs inducing seq errors))
* Compare SNPs and indels separately.
    * Indels have lower total counts but higher
      error rates.
* Distinguish between low coverage variants and actual diff in variant assessment
* Evaluate only where coverage > 4.

#### Aligners
* novoalign vs. bwamem (GATK post-align processing + GATK ug varcaller)
    * 1389 more concordant SNPs and 145 indels not seen with novoalign.
        * 1024 of these in novoalign low coverage regions
            * 941 of these had coverage `<` 10 with bwamem
* bwamem equally as good as novoalign in calling, and faster

#### Post-alignment process
* GATK vs. gkno
    * GATK: dedup picard, GATK BQSR, GATK indel realign
    * gkno: dedup samtools, ogap indel realign. No BQSR.
* GATK better than gkno:
    * better SNP calling with BQSR
    * 1% of SNPs missed due to poor quality calculations

#### Varcallers
* FreeBayes vs. GATK ug vs. GATK hc
* Biggest impact on final set of variants
* Highest difference of concordant SNPs between ug and hc.
    * ug best at detecting SNPs
    * hc best at detecting indels

Updated Evaluation framework (2013-05-06)
---------------------------------
* [Blog post](<https://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/>)

* FreeBayes detects more concordant SNPs and indels compared to GATK ug/hc
* BQSR and indel realignment have little impact
* The Ensemble calling method provides the best variant detection by combining
  inputs from GATK ug, hc and FreeBayes

WGS trio varcall evaluation (2014-05-12)
---------------------------------
* [Blog post](<https://bcbio.wordpress.com/2014/05/12/wgs-trio-variant-evaluation/>)

### Method
* Used NA12878/NA12891/NA12892 trio WGS 50x Illumina platinum genomes
* Alignment with bwa mem, varcalls with FreeBayes/GATK hc
* Evaluated calls using GIAB reference for NA12878

### Filters

#### Low complexity regions
* LCRs cover 2% of the genome
* Consist of locally repetitive genomic sections
* LCRs contribute:
    * 2% of SNPs
    * 45% of indels
* So: exclude LCRs in variant comparisons

#### High depth and low quality
* Reduce false positives with QUAL `<` 500 and DP `>` 208

#### GATK
* Overall VQSR provides good filtering
* Hard filters are also ok

Cancer variant callers (2015-03-05)
---------------------------------
* [Blog post](<https://bcbio.wordpress.com/2015/03/05/cancerval/>)

### Problem
* Cancer variant calling is challenging
* Tumor samples have mixed cellularity (contaminating normal sample)
* Multiple sub-clones with different somatic variants
* Low-frequency sub-clonal variants are critical but hard to detect
* Reference sets don't exist for cancer calling in public (e.g. NA12878 GIAB),
  but:
    * DREAM synthetic datasets exist with cellularity and multiple sub-clones
    * DREAM real non-simulated data also exist
    * ICGC, Bina


* Use DREAM synthetic dataset 3 to evaluate tumor/normal varcalling
* Ensemble callset generated with good sensitivity and precision

### Variant Callers
* MuTect (Broad)
* VarDict (AstraZeneca)
* FreeBayes (Eric Garrison - Marth lab)
* VarScan (Dan Koboldt)
* Scalpel (Schatz lab)
* LUMPY (Ryan Layer - Quinlan lab)
* DELLY (Tobias Rausch - EMBL)
* WHAM (Zev Kronenberg - Yandell lab)

* Ensemble method combines SNPs, indels and SVs
    * Simplified approach performs well and faster than previous SVM approach

### Tumor-only prioritisation
* Sometimes there aren't matched normal samples for filtration
    * E.g. FFPE samples and tumor cell lines

* Prioritise likely tumor specific variants using publicly available resources:
  COSMIC, ClinVar, 1KG, ESP, ExAC

* High/medium impact, MAF `<` 1% in 1KG/ExAC, in COSMIC/ClinVar

* Can't expect to filter private germline mutations that aren't in public
  databases
