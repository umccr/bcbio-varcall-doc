The commands run by bcbio when using `VarDict` as the variant caller can be
found in the `bcbio-nextgen-commands.log` log file:

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

__chunk1__

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
    * output has one set of columns for each of the two samples.
* Pipe the output into...

```
| testsomatic.R
```

__chunk2__
* `testsomatic.R` "performs a statistical test for strand bias"


* Pipe the output into...
```
| var2vcf_paired.pl
-P 0.9
-m 4.25
-f 0.1
-M
-N "tumor_downsample|control_downsample"
```

__chunk3__

* chunk4
```
| /data/projects/punim0010/local/share/bcbio/anaconda/bin/py -x
'bcbio.variation.vcfutils.add_contig_to_header(x, "/data/projects/punim0010/local/share/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa")'
```

* chunk5
```
| bcftools filter -m '+' -s 'REJECT' -e 'STATUS !~ ".*Somatic"' 2> /dev/null
```

* chunk6
```
| /data/projects/punim0010/local/share/bcbio/anaconda/bin/py -x
'bcbio.variation.vardict.depth_freq_filter(x, 0, "False")'
```
```
| bcftools filter -i 'QUAL >= 0'
```
```
| sed 's/\\.*Somatic\\/Somatic/'
```
```
| sed 's/REJECT,Description=".*">/REJECT,Description="Not Somatic via VarDict">/'
```
```
| /data/projects/punim0010/local/share/bcbio/anaconda/bin/python -c
'from bcbio.variation import freebayes; freebayes.call_somatic("tumor_downsample", "control_downsample")' 
```
```
| awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $4) } {print}'
```
```
| awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $5) } {print}'
```
```
| awk -F$'\t' -v OFS='\t' '$1!~/^#/ && $4 == $5 {next} {print}'
```
```
| /data/projects/punim0010/local/share/bcbio/anaconda/bin/vcfstreamsort
```
```
| bgzip -c >
/data/projects/punim0010/projects/Diakumis_bcbio_varcall_doc/git/scripts/bcbio/2017-08-28_cancer-normal/work/bcbiotx/tmpyVERSl/batch1-21_19329108_42402896-raw.vcf.gz
```

