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

* argument is the read output from...

## samtools1

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

* argument is the read output from...

## samtools2

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

## varscan fix indels

```
cat work/bcbiotx/tmpWTXUUy/batch1-21_19329108_42402896-indel.vcf \
  | /data/projects/punim0010/local/share/bcbio/anaconda/bin/py -x \
  'bcbio.variation.varscan.fix_varscan_output(x, "control_downsample", "tumor_downsample")' \
  | awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $4) } {print}' \
  | awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $5) } {print}' \
  | ifne vcfuniqalleles \
  | /data/projects/punim0010/local/share/bcbio/anaconda/bin/py -x \
  'bcbio.variation.vcfutils.add_contig_to_header(x,
  "/data/projects/punim0010/local/share/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa")' \
  | bcftools filter \
  -m + -s REJECT \
  -e "SS != '.' && SS != '2'" 2> /dev/null \
  | /data/projects/punim0010/local/share/bcbio/anaconda/bin/py -x \
  'bcbio.variation.varscan.spv_freq_filter(x, 1)' \
  | bgzip -c > work/bcbiotx/tmpIJ1Wc3/batch1-21_19329108_42402896-indel-fix.vcf.gz
```

## tabix

```
tabix -f -p vcf work/bcbiotx/tmpb4NtZh/batch1-21_19329108_42402896-snp-fix.vcf.gz
tabix -f -p vcf work/bcbiotx/tmp0SPSvI/batch1-21_19329108_42402896-indel-fix.vcf.gz
```
