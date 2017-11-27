The commands run by bcbio when using `VarDict` as the variant caller can be
found in the `bcbio-nextgen-commands.log` log file. In summary:

<!-- vim-markdown-toc GFM -->

* [vardict-java](#vardict-java)
* [testsomatic](#testsomatic)
* [var2vcf](#var2vcf)
* [add-contig](#add-contig)
* [bcftools filter1](#bcftools-filter1)
* [depth-freq-filter](#depth-freq-filter)
* [bcftools filter2](#bcftools-filter2)
* [seds](#seds)
* [freebayes-call-somatic](#freebayes-call-somatic)
* [awks](#awks)
* [vcfstreamsort](#vcfstreamsort)
* [bgzip](#bgzip)

<!-- vim-markdown-toc -->

### vardict-java

```
export VAR_DICT_OPTS='-Xms750m -Xmx3500m -XX:+UseSerialGC -Djava.io.tmpdir=/data/projects/punim0010/projects/Diakumis_bcbio_varcall_doc/git/scripts/bcbio/2017-08-28_cancer-normal/work/bcbiotx/tmpyVERSl' &&
vardict-java
-G /data/projects/punim0010/local/share/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
-f 0.1
-N tumor_downsample
-b "/data/projects/punim0010/projects/Diakumis_bcbio_varcall_doc/git/scripts/bcbio/2017-08-28_cancer-normal/work/prealign/tumor_downsample/tumor_100x-ready_chr21-dedup.bam|/data/projects/punim0010/projects/Diakumis_bcbio_varcall_doc/git/scripts/bcbio/2017-08-28_cancer-normal/work/prealign/control_downsample/control_100x-ready_chr21-dedup.bam"
-c 1
-S 2
-E 3
-g 4
-Q 10
-F 0x700
/data/projects/punim0010/projects/Diakumis_bcbio_varcall_doc/git/scripts/bcbio/2017-08-28_cancer-normal/work/vardict/21/batch1-21_19329108_42402896-regions-regionlimit.bed
```

* Specify min/max mem and tmp dir in `work/bcbiotx/`
* Run vardict-java with below options:
    * `-G GRCh37.fa`: input `FASTA` file. Needs to be indexed (`.fai`).
    * `-f 0.1`: allele frequency threshold.
    * `-N tumor_downsample`: sample name (not sure how it is picked up when two
      samples).
    * `-b "tumor.bam|control.bam"`: input `BAM` files. This way `VarDict` is run
      in paired mode. Need to be indexed (`.bai`).
    * `-c 1`: the column for chromosome.
    * `-S 2`: the column for region start.
    * `-E 3`: the column for region end.
    * `-g 4`: the column for gene name.
    * `-Q 10`: filter out reads with mapping quality less than 10.
    * `-F 0x700`: hex flag to filter out reads. `0x700` means:
        * not primary alignment
        * read fails platform/vendor qc
        * read is PCR or optical duplicate
    * `target_regions.bed`: target regions in `BED` format. bcbio splits
      chr into chunks so e.g. for chr21 there were three separate commands run
      with _nearly_ contiguous `BED` coordinates (based on file names).
    * Output:
        * Total of 51 columns
        * Cols 1-7: sample, gene, chr, start, end, ref, alt
        * Cols 8-13: For Sample 1: depth-tot, depth-totalt, depth-fwdref, depth-revref,
          depth-fwdalt, depth-revalt
        * Col 14-25: For Sample 1: genotype, AF, bias, pmean, pstd, qual, qstd, mapq, sn, hiaf, adjaf, nm
        * Col 26-31: For Sample 2: depth-tot, depth-totalt, depth-fwdref, depth-revref,
          depth-fwdalt, depth-revalt
        * Col 32-43: For Sample 2: genotype, AF, bias, pmean, pstd, qual, qstd, mapq, sn, hiaf, adjaf, nm
        * Col 44-51: shift3, msi, msilen, left-flank-seq, right-flank-seq, segment, status, type

* Pipe the output into...

### testsomatic

```
| testsomatic.R
```

* "Performs a statistical test for strand bias"
* Runs a Fisher's exact test four times for each input row:
    1. Use columns 10-13 for counts. Keep pvalues1 and oddratio1.
    2. Use columns 28-31 for counts. Keep pvalues2 and oddratio2.
    3. Use columns 9, tref, 27 and rref for counts (greater)
    4. Use columns 9, tref, 27 and rref for counts (less)
* Keep the smaller pvalue from the greater/less test in pvalues, and the odds
  ratio in oddratio.
* Insert the final columns into the input data frame:
    * d[, 1:25], pvalues1, oddratio1, d[, 26:43], pvalues2, oddratio2, d[, 44:dim(d)[2]], pvalues, oddratio
* Pipe the output into...

### var2vcf

```
| var2vcf_paired.pl
-P 0.9
-m 4.25
-f 0.1
-M
-N "tumor_downsample|control_downsample"
```

* Run `var2vcf_paired.pl` with following options:
    * `-P 0.9`:  maximum p-value
    * `-m 4.25`: max mean mismatches allowed
    * `-f 0.1`: minimum AF
    * `-M`: output only candidate somatic
    * `-N "tumor_downsample|control_downsample"`: sample names.
* Create  `%hash` where key is the chromosome, value is a hash where key is the
  start position, value is an array of arrays, each containing a single row from
  the previous dataframe. Example:

```
%hash = {
  '21' => {
    '12345' => [
                 [
                  'tumor_downsample', '21', '21', '12345', '12345',
                  'A', 'C', '2', '2', '0', '0', '1', '1', 'C/C', '1', ...
                  ]
                ],
    '12567' => [
                 [
                  'tumor_downsample', '21', '21', '12567', '12567',
                  'G', 'T', '4', '4', '0', '0', '1', '3', 'T/T', '1', ...
                  ]
                ],
```

* `@chrs` contains all the keys of `%hash`.
* Print a VCF header, followed by column names.
* Foreach chromosome in the hash keys:
    * Foreach position in the position keys:
        * tmp is the array/row of each position
            * loop over each position. Description of fields:

| Field        | Description |
|--------------|-------------|
| `$sample`    | Sample name |
| `$gene`      | Gene |
| `$chrt`      | Chromosome |
| `$start`     | Start Pos |
| `$end`       | End Pos |
| `$ref`       | Ref base |
| `$alt`       | Alt base |
| `$dpX`       | Total depth of coverage (for sample X) |
| `$vdX`       | Variant depth of coverage |
| `$rfwdX`     | Depth fwd ref |
| `$rrevX`     | Depth rev ref |
| `$vfwdX`     | Depth fwd alt |
| `$vrevX`     | Depth rev alt |
| `$gtX`       | Genotype |
| `$afX`       | Allele frequency |
| `$biasX`     | Strand bias |
| `$pmeanX`    | Mean position in reads |
| `$pstdX`     | Position STD in reads |
| `$qualX`     | Mean of the base (Phred) qualities of the bases for the variant |
| `$qstdX`     | Standard deviation as above |
| `$mapqX`     | Mapping quality (mean of mapping qualities of reads supporting variants |
| `$snX`       | Signal to noise |
| `$hiafX`     | AF using only high-quality bases |
| `$adjafX`    | Adj AF for indels due to local realignment |
| `$nmX`       | Mean mismatches in reads |
| `$sbfX`      | Strand bias Fisher p-value |
| `$oddratioX` | Strand bias Odds Ratio |
| `$shift3`    | No. of bases to be shifted to 3' for dels due to alternative alignment |
| `$msi`       | Microsatellite. If greater than 1 indicates Microsatellite Instability (MSI) |
| `$msilen`    | MSI unit repeat length in bp |
| `$lseq`      | 5' flanking sequence |
| `$rseq`      | 3' flanking sequence |
| `$seg`       | Segment |
| `$status`    | Somatic or germline status |
| `$type`      | Variant type (SNV, Insertion, Deletion, Complex) |
| `$pvalue`    | p-value |
| `$oddratio`  | Odds Ratio |

* Create two arrays `@filters` and `@filters2`, pushing into them info based on
  the value of the above fields for sample 1 and 2.
* Info checked:
    * mindepth1, minvardepth1 unless the status is strongsomatic etc.
    * mindepth2, minvardepth2
    * af1, pmean1, pstd1, qual1, mapq1, sn1, nm1, bias1
    * af2, pmean2, pstd2, qual2, mapq2, sn2, nm2, bias2
    * msi, bias
    * if the `-M` option was used, push pvalue etc.
* If there have been any filters pushed into `@filters`, join them into a
  `;`-separated string. Else, use `PASS`. Assign result to `$filter`.
* If the `-M` option was used, `PASS` the variant if it's good in the germline
  sample (`@filters2` = 0).
* Assign the genotype for both samples after comparing the allele frequencies:
    * the default `$GTFREQ` is 0.2
    * the default `$FREQ` (AF) is 0.02 (bcbio uses 0.1)
    * if 1-afX < GTFREQ => `1/1`, else
    * if afX >= 0.5 => `1/0`, else
    * if afX >= FREQ => `0/1`, else `0/0`
* Print pinfo1, pfilter and pinfo2
* Pipe the output vcf into...

### add-contig

```
| /data/projects/punim0010/local/share/bcbio/anaconda/bin/py -x
'bcbio.variation.vcfutils.add_contig_to_header(x, "/data/projects/punim0010/local/share/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa")'
```

* Runs the `add_contig_to_header` function in `bcbio-nextgen/bcbio/variation/vcfutils.py`
* The `file_contigs` function that is called is in `bcbio-nextgen/bcbio/bam/ref.py`
* Pipe output into...

### bcftools filter1

```
| bcftools filter -m '+' -s 'REJECT' -e 'STATUS !~ ".*Somatic"' 2> /dev/null
```

* Run `bcftools filter` with following options:
    * `-m '+'`: append new FILTER strings (i.e. don't replace them)
    * `-s 'REJECT'`: annotate FILTER column with 'REJECT'
    * `-e EXPRESSION`: exclude sites for which EXPRESSION is true. So filter out
      sites where STATUS does not contain the 'Somatic' string.
* Pipe output into...

### depth-freq-filter

```
| /data/projects/punim0010/local/share/bcbio/anaconda/bin/py -x
'bcbio.variation.vardict.depth_freq_filter(x, 0, "False")'
```

Runs the `depth_freq_filter` function in
`bcbio-nextgen/bcbio/variation/vardict.py`:

* First add two FILTER headers to the vcf file (above the column header line)
* Then if the variant doesn't pass several filters, append the appropriate flag
  to the FILTER column, or replace if already was '.' or 'PASS' (read code for
  what filters are specified).
* Pipe output into...

### bcftools filter2

```
| bcftools filter -i 'QUAL >= 0'
```

* Include only sites with non-negative QUAL

### seds

```
| sed 's/\\.*Somatic\\/Somatic/'
```

* Get a clean 'somatic' status

```
| sed 's/REJECT,Description=".*">/REJECT,Description="Not Somatic via VarDict">/'
```

* Clean up FILTER row for REJECT ID

### freebayes-call-somatic

```
| /data/projects/punim0010/local/share/bcbio/anaconda/bin/python -c
'from bcbio.variation import freebayes; freebayes.call_somatic("tumor_downsample", "control_downsample")' 
```

* For VarDict I believe it simply checks the AF of the tumor compared to the
  normal based on a threshold. If there is supporting evidence, adds 'SOMATIC'
  to the INFO column. Else, adds 'REJECT' to the FILTER column.
* Pipe output into...

### awks

```
| awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $4) } {print}'
```

* Cleans up Ref column

```
| awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $5) } {print}'
```

* Cleans up Alt column

```
| awk -F$'\t' -v OFS='\t' '$1!~/^#/ && $4 == $5 {next} {print}'
```

* I *think* above filters out rows where Ref and Alt are the same.

* Pipe output into...

### vcfstreamsort

```
| /data/projects/punim0010/local/share/bcbio/anaconda/bin/vcfstreamsort
```

* Part of `vcflib`. Reads VCF on stdin and guarantees that the positional order
  is correct provided out-of-order variants are no more than 100 positions in
  the VCF file apart.
* Pipe output into...

### bgzip

```
| bgzip -c >
/data/projects/punim0010/projects/Diakumis_bcbio_varcall_doc/git/scripts/bcbio/2017-08-28_cancer-normal/work/bcbiotx/tmpyVERSl/batch1-21_19329108_42402896-raw.vcf.gz
```

* bgzips final vcf file

