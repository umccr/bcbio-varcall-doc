The commands run by bcbio when using `FreeBayes` as the variant caller can be
found in the `bcbio-nextgen-commands.log` log file. In summary:

### freebayes

```
freebayes \
  -f GRCh37.fa \ 
  --genotype-qualities \ # 
  --strict-vcf \
  --ploidy 2 \
  --targets work/freebayes/21/batch1-21_0_19327180-regions.bed \
  --min-repeat-entropy 1 \
  --no-partial-observations \
  --min-alternate-fraction 0.1 \
  --pooled-discrete \
  --pooled-continuous \
  --report-genotype-likelihood-max \
  --allele-balance-priors-off \
  work/prealign/tumor_downsample/tumor_100x-ready_chr21-dedup.bam \
  work/prealign/control_downsample/control_100x-ready_chr21-dedup.bam |
```

* Options:
  * `-f GRCh37.fa`: reference FASTA
  * `--genotype-qualities`: calculate marginal prob of genotypes and report as
    GQ for each sample
  * `--strict-vcf`: generate strict VCF format (FORMAT/GQ will be an int)
  * `--ploidy 2`: ploidy...
  * `--targets work/freebayes/21/batch1-21_0_19327180-regions.bed`: limit
    analysis to targets listed in this `bed` file
  * `--min-repeat-entropy 1`: To detect interrupted repeats, build across
    sequence until it has entropy > 1 bits per bp
  * `--no-partial-observations`: exclude observations which do not fully span the
    dynamically-determined detection window
  * `--min-alternate-fraction 0.1`: require at least this fraction of
    observations supporting an alternate allele within a single individual in
    order to evaluate the position
  * `--pooled-discrete`: Assume that samples result from pooled sequencing.
    Model pooled samples using discrete genotypes across pools. When using this
    flag, set --ploidy to the number of alleles in each sample.
  * `--pooled-continuous`: Output all alleles which pass input filters,
    regardles of genotyping outcome or model.
  * `--report-genotype-likelihood-max`: Report genotypes using the
    maximum-likelihood estimate provided from genotype likelihoods.
  * `--allele-balance-priors-off`: Disable use of aggregate probability of
    observation balance between alleles as a component of the priors.
  * `tumor.bam control.bam`: `BAM` files to be analysed

### bcftools filter

```
bcftools filter -i 'ALT="<*>" || QUAL > 5' |
```

* Include only sites where `ALT` has something or `QUAL` > 5.

### freebayes call-somatic

```
python -c 'from bcbio.variation import freebayes; freebayes.call_somatic("tumor_downsample", "control_downsample")' |
```

* `call_somatic` function found in
  [bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/variation/freebayes.py)
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

### awk

```
awk -F$'\t' -v OFS='\t' \
  '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $4) } {print}' |
```

* Just clean up the REF field in the VCF file

### bcftools view

```
bcftools view -s tumor_downsample,control_downsample -a - |
```

* `-s x,y`: Include both samples x and y
* `-a`: trim alternate alleles not seen in subset

### freebayes remove-missingalt

```
py -x 'bcbio.variation.freebayes.remove_missingalt(x)' |
```

* `remove_missingalt` function found in
  [bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/variation/freebayes.py)
* "Removes lines that are missing an alternative allele. During cleanup of extra alleles,
  bcftools has an issue in complicated cases with duplicate alleles and will end up
  stripping all alternative alleles. This removes those lines to avoid issues
  downstream."

### vcflib and vt

```
vcfallelicprimitives -t DECOMPOSED --keep-geno |
```

* Part of [vcflib](https://github.com/vcflib/vcflib)
* If multiple allelic primitives (gaps or mismatches) are specified in a single
  VCF record, split the record into multiple lines, but drop all INFO fields.
  "Pure" MNPs are split into multiple SNPs unless the `-m` flag is provided.
  Genotypes are phased where complex alleles have been decomposed, provided
  genotypes in the input.
    * `-t DECOMPOSED`: Tag records which are split apart of a complex allele with
      the flag `DECOMPOSED`.
    * `--keep-geno`: Maintain genotype-level annotations when decomposing.

```
vcffixup - |
```

* Part of [vcflib](https://github.com/vcflib/vcflib)
* Count the allele frequencies across alleles present in each record in the VCF
  file. (Similar to vcftools --freq.) Uses genotypes from the VCF file to
  correct AC (alternate allele count), AF (alternate allele frequency),
  NS (number of called), in the VCF records.

```
vcfstreamsort |
```

* Part of [vcflib](https://github.com/vcflib/vcflib)
* Reads VCF on stdin and guarantees that the positional order is correct
  provided out-of-order variants are no more than 100 positions in the VCF file
  apart.

```
vt normalize -n -r GRCh37.fa -q - |
```

* Part of [vt](https://github.com/atks/vt)
* Normalise variants (see <http://genome.sph.umich.edu/wiki/Vt#Normalization>)
    * `-n`: no fail on reference inconsistency.
      There is an underlying assumption that the REF field is consistent with the
      reference sequence use, vt will check for this and will fail if reference
      inconsistency is encountered; this may be relaexd with the -n option.
    * `-r`: reference FASTA file
    * `-q`: quiet

```
vcfuniqalleles
```

* Part of [vcflib](https://github.com/vcflib/vcflib)
* For each record, remove any duplicate alternate alleles that may have resulted
  from merging separate VCF files

```
vt uniq - 2> /dev/null |
```

* Part of [vt](https://github.com/atks/vt)
* Drops duplicate variants that appear later in the file.

```
bgzip -c > work/bcbiotx/tmp68Zr8g/batch1-21_0_19327180-raw.vcf.gz
```

* output final `bgzip` file
