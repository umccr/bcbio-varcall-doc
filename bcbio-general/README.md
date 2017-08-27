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
* [Running](#running)
* [Project Structure](#project-structure)
* [Configuration Notes](#configuration-notes)
    * [Variant calling](#variant-calling)
* [Code](#code)
* [Miscellaneous](#miscellaneous)

<!-- vim-markdown-toc -->

Installation
------------

First download the install script:

```
wget https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
```

Then run it:

```
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

```
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
- **tumor-only** samples: filter out likely germline variants present in
  1KG/ExAC and not present in COSMIC.
  Likely germline variants are marked with a `LowPriority` filter.
  Two variant outputs are generated:
    - `sample-caller.vcf.gz` contains the somatic calls.
    - `sample-caller-germline.vcf.gz` contains likely germline mutations.
- **tumor + normal** samples: call somatic (tumor-specific) and germline
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

Running
-------

Output from `bcbio_nextgen.py -h`:

```
usage: bcbio_nextgen.py [-h] [-n NUMCORES] [-t {local,ipython}]
                        [-s {lsf,sge,torque,slurm,pbspro}]
                        [--local_controller] [-q QUEUE] [-r RESOURCES]
                        [--timeout TIMEOUT] [--retries RETRIES] [-p TAG]
                        [-w WORKFLOW] [--workdir WORKDIR] [-v]
                        [--force-single]
                        [global_config] [fc_dir] [run_config [run_config ...]]

Community developed high throughput sequencing analysis.

positional arguments:
  global_config         Global YAML configuration file specifying details
                        about the system (optional, defaults to installed
                        bcbio_system.yaml)
  fc_dir                A directory of Illumina output or fastq files to
                        process (optional)
  run_config            YAML file with details about samples to process
                        (required, unless using Galaxy LIMS as input)

optional arguments:
  -h, --help            show this help message and exit
  -n NUMCORES, --numcores NUMCORES
                        Total cores to use for processing
  -t {local,ipython}, --paralleltype {local,ipython}
                        Approach to parallelization
  -s {lsf,sge,torque,slurm,pbspro}, --scheduler {lsf,sge,torque,slurm,pbspro}
                        Scheduler to use for ipython parallel
  --local_controller    run controller locally
  -q QUEUE, --queue QUEUE
                        Scheduler queue to run jobs on, for ipython parallel
  -r RESOURCES, --resources RESOURCES
                        Cluster specific resources specifications. Can be
                        specified multiple times. Supports SGE, Torque, LSF
                        and SLURM parameters.
  --timeout TIMEOUT     Number of minutes before cluster startup times out.
                        Defaults to 15
  --retries RETRIES     Number of retries of failed tasks during distributed
                        processing. Default 0 (no retries)
  -p TAG, --tag TAG     Tag name to label jobs on the cluster
  -w WORKFLOW, --workflow WORKFLOW
                        Run a workflow with the given commandline arguments
  --workdir WORKDIR     Directory to process in. Defaults to current working
                        directory
  -v, --version         Print current version
  --force-single        Treat all files as single reads
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

* `config`: configuration files for input samples
* `final`: output from pipeline
* `work`: processing
* `log`: contains log files:
    * `bcbio-nextgen.log`: overview
    * `bcbio-nextgen-debug.log`: details
    * `bcbio-nextgen-commands.log`: full command lines used


Configuration Notes
-------------------

You can create a sample configuration file with e.g.:

```
bcbio_nextgen.py -w template gatk-variant project1.csv sample1.bam sample2_1.fq sample2_2.fq
```

where:

* -w template: generate a template
* `gatk-variant`: name of yaml in `bcbio-nextgen/config/templates` or path
  to custom file.
* `project1.csv`: file with sample/algorithm metadata
* `bam`/`fq`: input `BAM`/`FASTQ` files.

**OR**

You can also do the above in a more simple manner:

```
bcbio_nextgen.py -w template tumor-paired project1
```

which outputs:

```
Template configuration file created at: ./proj1/config/proj1-template.yaml
Edit to finalize custom options, then prepare full sample config with:
bcbio_nextgen.py -w template ./proj1/config/proj1-template.yaml proj1 sample1.bam sample2.fq
```

The different options in the configuration file that can be used are:

* `fc_date`, `fc_name`: combined to form a prefix of intermediate files.
* `upload`: which directory to put the output files.
* `details`: list of sections, where each section describes a sample to process.
* `description`: name of sample - used in final output.
* `analysis`: which pipeline to run.
* `algorithm`: tune the analysis pipeline.
* `quality_format`: FASTQ quality.

### Variant calling
See [bcbio doc](https://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#variant-calling).

The `variantcaller` option can take a list with multiple callers
(`false`, freebayes, gatk-haplotype, haplotyper, platypus, mutect, mutect2,
scalpel, tnhaplotyper, tnscope, vardict, varscan, samtools, gatk).

Code
-----
See [bcbio doc](https://bcbio-nextgen.readthedocs.io/en/latest/contents/code.html).

The following two directories contain code of most interest for getting started:

* `bcbio-nextgen/bcbio/pipeline`, starting from `main.py`
* `bcbio-nextgen/bcbio/variation`


Miscellaneous
-------------
* Installation directory on Spartan: `/data/projects/punim0010/local/share/bcbio`
