The commands run by bcbio when using `VarScan` as the variant caller can be
found in the `bcbio-nextgen-commands.log` log file. In summary:

## varscan
```
varscan \
  -Xms681m -Xmx3181m \
  -Duser.language=en \
  -Duser.country=US \
  -XX:+UseSerialGC \
-Djava.io.tmpdir=work/bcbiotx/tmpuV45Kf \
somatic
```

* Run the `somatic` module from VarScan, which takes two pileups as arguments:

### pileup1

```
<(
samtools mpileup \
  -f GRCh37.fa \
  -d 1000 \
  -L 1000 \
  -l work/varscan/21/batch1-21_19329108_42402896-regions.bed \
  work/prealign/control_downsample/control_100x-ready_chr21-dedup.bam \
  | { ifne grep -v -P '\t0\t\t$' || true; }
  )
```

* `-f GRCh37.fa`: reference sequence file.
* `-d 1000`: max depth.
* `-L 1000`: max depth for indel calling.
* `-l regions.bed`: skip positions not in this `bed` file
* `control.bam`: input normal `BAM` file
* if the input stream is not empty, output lines which don't end with 0
  tab-tab...


### pileup2

```
<(
samtools mpileup \
  -f GRCh37.fa \
  -d 1000 \
  -L 1000 \
  -l work/varscan/21/batch1-21_19329108_42402896-regions.bed \
  work/prealign/tumor_downsample/tumor_100x-ready_chr21-dedup.bam \
  | { ifne grep -v -P '\t0\t\t$' || true; }
  ) \
```

* Same as previous pileup, only this one is for the tumor sample.

## varscan additional args

```
  --output-snp work/bcbiotx/tmpYw7los/batch1-21_19329108_42402896-snp.vcf \
  --output-indel work/bcbiotx/tmpYw7los/batch1-21_19329108_42402896-indel.vcf \
  --output-vcf \
  --min-coverage 5 \
  --p-value 0.98 \
  --strand-filter 1 \
  --min-var-freq 0.1
```

* `--output-snp snp.vcf`: output file for SNP calls
* `--output-indel indels.vcf`: output file for indel calls
* `--min-coverage 5`: min cov in normal and tumor to call variant
* `--p-value 0.98`: p-value threshold to call a heterozygote
* `--strand-filter 1`: if set to 1, removes variants with `>90%` strand bias
* `--min-var-freq 0.1`: min variant frequency to call a heterozygote

## varscan fix snps

```
cat work/bcbiotx/tmpWTXUUy/batch1-21_19329108_42402896-snp.vcf \
  | /data/projects/punim0010/local/share/bcbio/anaconda/bin/py -x \
  'bcbio.variation.varscan.fix_varscan_output(x, "control_downsample", "tumor_downsample")' \
  | awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $4) } {print}' \
  | awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $5) } {print}' \
  | ifne vcfuniqalleles \
  | /data/projects/punim0010/local/share/bcbio/anaconda/bin/py -x \
  'bcbio.variation.vcfutils.add_contig_to_header(x,
  "/data/projects/punim0010/local/share/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa")' \
  | bcftools filter \
    -m + \
    -s REJECT \
    -e "SS != '.' && SS != '2'" 2> /dev/null \
  | /data/projects/punim0010/local/share/bcbio/anaconda/bin/py -x \
  'bcbio.variation.varscan.spv_freq_filter(x, 1)' \
    | bgzip -c > work/bcbiotx/tmpqwbclp/batch1-21_19329108_42402896-snp-fix.vcf.gz
```

* Pipe `snp.vcf` into the first argument of the `fix_varscan_output` function
  which can be found in `bcbio-nextgen/bcbio/variation/varscan.py`. It takes a
  further two arguments `"control_downsample"` and `"tumor_downsample"`. This
  function:

  > "Fixes the ALT column and also fixes floating point values output as
    strings: FREQ, SSC".

* The `awk` commands simply clean the REF and ALT columns.
* `vcfuniqalleles` is part of `vcflib`: remove any duplicate alternate alleles
  that may have resulted from merging separate VCF files.
* Pipe into the first argument of the `add_contig_to_header` function
  which can be found in `bcbio-nextgen/bcbio/variation/vcfutils.py`. Second
  argument is the ref seq FASTA file.
* Pipe into `bcftools filter`:
    * `-m +`: append new FILTER strings (i.e. don't replace them)
    * `-s REJECT`: annotate FILTER column with 'REJECT'
    * `-e EXPRESSION`: excludes sites for which EXPRESSION is true. So filter
      out sites where SS is not equal to '.' or '2'.
* Pipe into the first argument of the `spv_freq_filter` function which can be
  found in `bcbio-nextgen/bcbio/variation/varscan.py`. This function:

  > "Filters VarScan calls based on the SPV value and frequency. Removes calls
  with SPV < 0.05 and a tumor FREQ > 0.35. False positives dominate these higher freq, low
  SPV calls. They appear to be primarily non-somatic/germline variants not
  removed by other filter".

* bgzip output vcf file

## varscan fix indels

* Same code as above, just for the `indels.vcf`

## tabix

```
tabix -f -p vcf work/bcbiotx/tmpb4NtZh/batch1-21_19329108_42402896-snp-fix.vcf.gz
tabix -f -p vcf work/bcbiotx/tmp0SPSvI/batch1-21_19329108_42402896-indel-fix.vcf.gz
```

* `tabix` both vcfs.
