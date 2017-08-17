bcbio general notes
===================

<!-- vim-markdown-toc GFM -->
* [Installation](#installation)
* [Upgrade](#upgrade)
* [Data](#data)
* [Pipelines](#pipelines)
    * [Germline variant calling](#germline-variant-calling)
        * [Basic germline calling](#basic-germline-calling)
        * [Population calling](#population-calling)
    * [Cancer variant calling](#cancer-variant-calling)
* [Project Structure](#project-structure)
* [Miscellaneous](#miscellaneous)

<!-- vim-markdown-toc -->

Installation
------------

First download the install script:

```bash
wget https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
```

Then run it:

```bash
python bcbio_nextgen_install.py \
  /usr/local/share/bcbio \
  --tooldir=/usr/local \
  --genomes GRCh37 \
  --aligners bwa --aligners bowtie2
```

* Data are installed in `/usr/local/share/bcbio`
* Tools are installed in `/usr/local`
    * If you don't have permission, use `~/local`
* You can install additional data and tools later using
  `bcbio_nextgen.py upgrade`
* Use `--nodata` to avoid installing genome data.
* After installation edit the config file in
  `/usr/local/share/bcbio/galaxy/bcbio_system.yaml`.

Upgrade
-------
You can upgrade the code, tools and data with:

```bash
bcbio_nextgen.py upgrade -u stable --tools --data
```

* `-u`: stable or development code from GitHub
* `--tools`: upgrade third party tools
* `--toolplus`: additional tools (GATK, MuTect)
* `--tooldir`: install in different tool directory
* `--data`: upgrade genome data
* `--datatarget`: customised installed data
* `--genomes` and `--aligners`: additional aligner indexes to download
  and prepare
* `--cores`: number of cores used for indexing

Data
----
By default `bcbio` will install data files for `variation`, `rnaseq` and
`smallrna`. You can sub-select one of these.

* Available data targets are: `variation`, `rnaseq`, `smallrna`,
  `gemini`, `cadd`, `vep`, `dbnsfp`, `dbscsnv`, `battenberg`, `kraken`.

Pipelines
---------
Several pipelines can be run using `bcbio`.

### Germline variant calling
* SNPs, indels and structural variants can be called for germline populations.
* Annotation with `snpEff` or `VEP`. `GEMINI` database prepared.

#### Basic germline calling
* You can use the `FreeBayes` or the `GATK HaplotypeCaller` template.
* The `GATK` template follows best practices, including BQSR, local realignment
  and HaplotypeCaller variant calling.
* You can enable SV calling.


#### Population calling
You can either use the joint calling or the batch calling method.
The former calls samples independently, then combines them together into
a single callset. The latter calls all samples together, but obviously will
have extremely large memory requirements for more than 100 samples.
Can work with `GATK HaplotypeCaller` or `FreeBayes`.

### Cancer variant calling
Summary:
- _tumor-only_ samples: filter out likely germline variants present in
  1KG/ExAC and not present in COSMIC.
  Likely germline variants are marked with a `LowPriority` filter.
  Two variant outputs are generated:
    - `sample-caller.vcf.gz` contains the somatic calls.
    - `sample-caller-germline.vcf.gz` contains likely germline mutations.
- _tumor + normal_ samples: call somatic (tumor-specific) and germline
  (pre-existing) variants. For example you can use `VarDict` (somatic) and
  `FreeBayes` (germline).
- The ensemble approach combines calls from multiple SNP and indel
  callers, and also flattens structural variant calls into a combined
  representation.

For example, to use 3 tools for somatic calling, 3 tools for germline calling,
create ensemble calls for both and include germline and somatic events from 2
structural variant callers you need to specify:

```
variantcaller:
    somatic: [vardict, varscan, mutect2]
    germline: [freebayes, gatk-haplotype, platypus]
ensemble:
    numpass: 2
svcaller: [manta, cnvkit]
```

Project Structure
-----------------
The structure of each project is generally recommended to be:

```
my-project/
    |-- config/
    |-- final/
    |-- work/
        |-- log/
```

* `config`: configuration files for input samples in here
* `final`: output from pipeline in here
* `work`: processing done here
* `log`: contains log files:
    * `bcbio-nextgen.log`: overview
    * `bcbio-nextgen-debug.log`: details
    * `bcbio-nextgen-commands.log`: full command lines used

Miscellaneous
-------------
* GIAB high-confidence set of reference calls: <http://arxiv.org/abs/1307.4661>
* Installation directory on Spartan: `/data/projects/punim0010/local/share/bcbio`
