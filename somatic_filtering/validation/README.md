Small variant validation
========================

This document summarizes tools and datasets that can be used for variant calling comparison with truth set and visualization.

## Tools

* [rtgeval](https://github.com/lh3/rtgeval) - wrapper for RTG's vcfeval.

* [Concordance between variant callers](https://github.com/vladsaveliev/concordance)

* [JS-based Venn diagrams for BED and VCF files](https://github.com/vladsaveliev/Venn)

* Bcbio:
  * [Examples of validation](https://github.com/bcbio/bcbio_validations)
  * [Somatic variant callers validation for tumor-only samples blog post](http://bcb.io/2015/03/05/cancerval/)

## Somatic validation datasets

### GiaB NA12878/NA24385 somatic-like mixture
Source: https://github.com/hbc/projects/tree/master/giab_somatic
Paired
WGS
Depth:
Tumor purity:
Samples:
Location in the FS:

### ICGC-TCGA DREAM challange
Text: https://www.synapse.org/#!Synapse:syn312572
Datasets:
* Set 3
* Set 4
Paired
WGS
Depth:
Tumor purity:
Samples:
Location in the FS:

### COLO829
Metastatic melanoma cell line
Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4837349/, 2016

AZ:
Justin: `We also have COLO829 (T/N) on the IDT Pan Cancer Panel along with the Horizon Standards.  I will share this as well in case it is of any use to you all.`
Miika: `/ngs/oncology/analysis/translation/TS_UK_0005__WGS_HD701_Colo829/bcbio/work/post_processing`



### ICGC MB
Paper: https://www.nature.com/articles/ncomms10001, 2015
Paired
WGS
Depth: T: 314x, N: 272x
Tumor purity: 95–98%
Samples:
Location in the FS:

### 
99X for COLO829 and 103X for the paired normal 


































