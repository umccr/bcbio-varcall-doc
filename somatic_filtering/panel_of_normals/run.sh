

for f in 1VD/na12878_2.vcf.gz 2VD/na12878_3.vcf.gz 3VD/na12878_3VD.vcf.gz 4KC/na12878_4KC.vcf.gz PG/na12878_1.vcf.gz ; do
   (echo "${f/.vcf.gz/.pon.vcf}" ;
    vcfanno -lua normals/code.lua normals/vcfanno.n10.toml $f | bgzip -c > ${f/.vcf.gz/.pon.vcf.gz} ;
    tabix ${f/.vcf.gz/.pon.vcf.gz} ;
    bcftools filter -i "INFO/PoN_CNT>=1" ${f/.vcf.gz/.pon.vcf.gz} -Oz -o ${f/.vcf.gz/.pon.n1.vcf.gz} ;
    bcftools filter -i "INFO/PoN_CNT>=2" ${f/.vcf.gz/.pon.vcf.gz} -Oz -o ${f/.vcf.gz/.pon.n2.vcf.gz} ;
    /home/vlad/validation/rtgeval.kit/run-eval \
     	-s /home/vlad/bcbio/genomes/Hsapiens/GRCh37/rtg/GRCh37.sdf \
        -b /home/vlad/bcbio/genomes/Hsapiens/GRCh37/validation/giab-NA12878-NA24385-somatic/truth_regions.bed \
        /home/vlad/bcbio/genomes/Hsapiens/GRCh37/validation/giab-NA12878-NA24385-somatic/truth_small_variants.vcf.gz \ 
        ${f/.vcf.gz/.pon.n1.vcf.gz}
    /home/vlad/validation/rtgeval.kit/run-eval \
     	-s /home/vlad/bcbio/genomes/Hsapiens/GRCh37/rtg/GRCh37.sdf \
        -b /home/vlad/bcbio/genomes/Hsapiens/GRCh37/validation/giab-NA12878-NA24385-somatic/truth_regions.bed \
        /home/vlad/bcbio/genomes/Hsapiens/GRCh37/validation/giab-NA12878-NA24385-somatic/truth_small_variants.vcf.gz \ 
        ${f/.vcf.gz/.pon.n2.vcf.gz}
    ) & 
done

bgzip na12878_1vd_vs_2vd-vardict-annotated-pon.vcf
tabix -p vcf na12878_1vd_vs_2vd-vardict-annotated-pon.vcf.gz

### Filtering

To filter out variants met at least twice in the panel:
```bash
bcftools filter -i "INFO/PoN_CNT>1" na12878_1vd_vs_2vd-vardict-annotated-pon.vcf.gz -Oz -o na12878_1vd_vs_2vd-vardict-annotated-pon-n1.vcf.gz
```

### Evaluation

Using Heng Li's [rtgeval](https://github.com/lh3/rtgeval) tool to compare against GiaB variants truth set. The tool normalizes variants, compares to truth set in specified regions, and reports numbers of FP/FN/TN.

Subsample to 1 main sample:
```bash
bcftools view -s 1VD na12878_1vd_vs_2vd-vardict-annotated-pon-n1.vcf.gz -Oz -o na12878_1vd_vs_2vd-vardict-annotated-pon-n1-sample.vcf.gz
```

```bash
cd /home/vlad/validation
./rtgeval.kit/run-eval -s /home/vlad/bcbio/genomes/Hsapiens/GRCh37/rtg/GRCh37.sdf \
    -b $HOME/bcbio/genomes/Hsapiens/GRCh37/validation/giab-NA12878-NA24385-somatic/truth_regions.bed \
    $HOME/bcbio/genomes/Hsapiens/GRCh37/validation/giab-NA12878-NA24385-somatic/truth_small_variants.vcf.gz \ 
    na12878_1vd_vs_2vd-vardict-annotated-pon-n1-sample.vcf.gz
```
