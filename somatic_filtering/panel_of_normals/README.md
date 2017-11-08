Panel of normals
================

Using a sample from a normal tissue as a baseline for cancer variant calling allows to get rid of pre-existing "germline" calls in the tissue. 
However, this method might still miss to filter out false positive somatic variants due to the following reasons:
- sequencing or alignment artefacts that failed to be filtered in "tumor"
- real germline variants, not called in "normal" due to insufficient coverage in normal samples in those regions (e.g. due to low GC, unbalanced structural variants).

Applying a set unrelated normals to filter false positives might help to reduce the FP number by partly addressing the second factor, becuase the poorly covered regions from the related normal might be sufficiently covered in unrelated normals.

Matt Eldridge used an SNV panel from 149 blood normals from the oesophageal ICGC project to flag variants that are observed in:
* at least 5 unrelated normals
* in each, minimum variant allele fraction of 0.05
* in each, minimum allele count of 3.
MuTect2 from DREAM challenge synthetic 4 dataset results:
	precision 0.91 -> 0.97 while the
	sensitivity is only marginally reduced, remaining at 0.75
ICGC benchmark datasets, MuTect2 and Strelka, but not in quite such spectacular fashion, possibly reflecting the relative similarity of the sequencing carried out for these datasets compared with that run on our oesophageal samples.

## Samples

We are building a panel from normal samples used in UMCCR. Currently, the following:
```
MH17B001P004
MH17B001P013
VPH52_Blood
VPH54_Blood
VPH56_Blood
VPH58_Blood
VPH59_Blood
VPH61_Blood
WPT-013
```

## Annotating

The following code can be used to annotate a VCF sample against the set of normals. It will add an `INFO` key `PoN_CNT`, showing the number of hits in the panel:
```bash
python make_vcfanno_cnf.py > vcfanno.toml
cd ~/validation/normals
vcfanno -lua code.lua vcfanno.toml MH17B001P013-germline-ensemble-annotated.vcf.gz > MH17B001P013-germline-ensemble-annotated-pon.vcf.gz
tabix -p vcf MH17B001P013-germline-ensemble-annotated-pon.vcf.gz
```

## Evaluation

Using Heng Li's [rtgeval](https://github.com/lh3/rtgeval) tool to compare against GiaB variants truth set. The tool normalizes variants, compares to truth set in specified regions, and reports numbers of FP/FN/TN.

```bash
cd /home/vlad/validation
./rtgeval.kit/run-eval -s /home/vlad/bcbio/genomes/Hsapiens/GRCh37/rtg/GRCh37.sdf -b $HOME/bcbio/genomes/Hsapiens/GRCh37/validation/giab-NA12878-NA24385-somatic/truth_regions.bed $HOME/bcbio/genomes/Hsapiens/GRCh37/validation/giab-NA12878-NA24385-somatic/truth_small_variants.vcf.gz MH17B001P013-germline-ensemble-annotated.vcf.gz
```

## NA12878/NA24385 somatic-like mixture


## Using 2 differrent NA12878 as a T/N pair

Two NA12878 samples running as a T/N pair ideally should yeild zero somatic variants. However, in real world the output will be non-empty. Since all real variants should cancel out, the remainings are assumed to be false positives. We count how many of them are removed by applying the panel of normals, making sure that it doesn't reduce the true positive count significantly.

```
Sample                   TP         FP         FN         Precision  Sensitivity
PG                       1,100,526  2,167,729  4,193      33.67%     99.62%
1VD                      1,103,426  2,170,848  1,293      33.70%     99.88%
1VD /PoN                 377,173    189,382    727,546    66.57%     34.14%
1VD/2VD                  491        6,095      1,104,228  7.46%      0.04%        0000.n10.vcf.gz
1VD/2VD /PoN.10 - 1hit   442        5,189      1,104,277  7.85%      0.04%        0000.n10.hit1.vcf.gz
1VD/2VD /PoN.10 - 2hits  506        5,341      1,104,213  8.65%      0.05%        0000.n10.hit2.vcf.gz
2VD                      1,103,564  2,171,013  1,155      33.70%     99.90%
3VD                      1,103,733  2,171,143  986        33.70%     99.91%
4KC                      1,103,368  2,171,035  1,351      33.70%     99.88%
Truth                    1,104,541  171        171        99.98%     99.98%
```







