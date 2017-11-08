Somatic variant calling
=======================

Matt Eldridge made [slides](https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2017/Day3/somatic_snv_filtering.html) on problems around small somatic variant calling and filtering. We summarized the filtering ideas into a [spreadsheet](https://docs.google.com/spreadsheets/d/1Xbz4nW76mofKb9ym3C725W035qkA7JuUu8_FvYSCOT0/edit#gid=0), where for problem and idea we provide a corresponding solution (if available) used in [CaVEMan](https://github.com/cancerit/cgpCaVEManWrapper), [bcbio](http://bcb.io/2015/03/05/cancerval/) and [VarDict](somatic_filtering/vardict_filtering.md).


Callers used in bcbio
=====================

Somatic variant caller tools used by bcbio:

1. [VarDict](https://github.com/AstraZeneca-NGS/VarDict)
2. [MuTect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php)
3. [VarScan](http://dkoboldt.github.io/varscan/)
4. [FreeBayes](https://github.com/ekg/freebayes)
5. [Strelka2](https://github.com/Illumina/strelka)