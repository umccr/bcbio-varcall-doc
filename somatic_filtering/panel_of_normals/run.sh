NORMALS_DIR=/home/vlad/validation/normals

for f in mixture/24385-ensemble-annotated.vcf.gz ; do
   (echo "${f/.vcf.gz/.pon.vcf}" ;
    vcfanno -lua $NORMALS_DIR/code.lua $NORMALS_DIR/vcfanno.n10.toml $f | bgzip -c > ${f/.vcf.gz/.pon.vcf.gz} ;
    bcftools filter -e "INFO/PoN_CNT>=1" ${f/.vcf.gz/.pon.vcf.gz} -Oz -o ${f/.vcf.gz/.pon.n1.vcf.gz} ;
    bcftools filter -e "INFO/PoN_CNT>=2" ${f/.vcf.gz/.pon.vcf.gz} -Oz -o ${f/.vcf.gz/.pon.n2.vcf.gz} ;
    tabix -p vcf ${f/.vcf.gz/.pon.n1.vcf.gz} ;
    tabix -p vcf ${f/.vcf.gz/.pon.n2.vcf.gz} ;
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
