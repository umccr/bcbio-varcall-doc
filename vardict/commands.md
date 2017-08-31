The commands run by bcbio when using `VarDict` as the variant caller are:

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

* chunk1:
    * specify min/max mem and tmp dir.
    *


* chunk2
```
| testsomatic.R
```

* chunk3
```
| var2vcf_paired.pl -P 0.9 -m 4.25 -f 0.1  -M  -N "tumor_downsample|control_downsample"
```

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

