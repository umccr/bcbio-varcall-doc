Publication
-----------
* Applies a Bayesian classifier to detect somatic mutations with very low AF
* Carefully tuned filters that ensure high specificity
* Has higher sensitivity with similar specificity, especially for variants with
  `AF < 0.1`


### Intro
* SNVs are important for altering gene function in cancer, but hard to identify:
    * Very low frequency in genome: 0.1 - 100 mutations per Mb
    * Present in small fraction of targeted DNA molecules due to:
        * contaminating normal cells in tumor sample
        * CNV in cancer genome
        * Subclonality (presence in a subpopulation of tumor cells)

* AF: allelic fraction (fraction of DNA molecules harboring variant)

* Sensitivity and specificity depend on:
    * depth of coverage in tumor and normal sample
    * sequencing error rate
    * allelic fraction
    * thresholds used to declare a mutation

* Two benchmarking approaches:
    * Downsampling
    * Virtual tumors
