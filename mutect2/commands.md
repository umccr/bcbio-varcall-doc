The commands run by bcbio when using `MuTect2` as the variant caller can be
found in the `bcbio-nextgen-commands.log` log file. In summary (from
[docs](https://github.com/broadinstitute/gatk/tree/master/docs/mutect)):

  > Mutect2 emits candidate variants with a set of annotations. After that
    FilterMutectCalls produces filtered calls by subjecting these variants to a
    series of hard filters that reject sites if some annotation is out of an
    allowable range.

More info
[here](https://software.broadinstitute.org/gatk/gatkdocs/4.beta.1/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php).

## MuTect2
```
gatk-launch \
--javaOptions '-Xms681m -Xmx3181m -XX:+UseSerialGC -Djava.io.tmpdir=tmpR4gahi' \
Mutect2 \
-R GRCh37.fa \
--annotation ClippingRankSumTest \
--annotation DepthPerSampleHC \
--annotation FisherStrand \
--annotation MappingQualityRankSumTest \
--annotation MappingQualityZero \
--annotation QualByDepth \
--annotation ReadPosRankSumTest \
--annotation RMSMappingQuality \
--annotation MappingQuality \
--annotation DepthPerAlleleBySample \
--annotation Coverage \
--readValidationStringency LENIENT \
-I tumor_100x-ready_chr21-dedup.bam \
--tumorSampleName tumor_downsample \
-I control_100x-ready_chr21-dedup.bam \
--normalSampleName control_downsample \
-L work/mutect2/21/batch1-21_0_19327180-regions.bed \
--interval_set_rule INTERSECTION \
-ploidy 2 \
-O work/bcbiotx/tmpR4gahi/batch1-21_0_19327180-raw.vcf.gz &&
```

Source code found in the GATK
[GitHub repo](https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect)


* Inputs are specified that way because
([code](https://github.com/broadinstitute/gatk/blob/a482f09909b36f772e3dd4d1e5c030f15ee7ecc5/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/M2ArgumentCollection.java#L17)):

```
//TODO: HACK ALERT HACK ALERT HACK ALERT
//TODO: GATK4 does not yet have a way to tag inputs, eg -I:tumor tumor.bam
//-I:normal normal.bam,
//TODO: so for now we require the user to specify bams *both* as inputs,
//with -I tumor.bam -I normal.bam
//TODO: *and* as sample names e.g. -tumor tumorSampleName -normal
//normalSampleName
```

### Annotation Modules

* [ClippingRankSumTest](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_ClippingRankSumTest.php) -
Rank Sum Test for hard-clipped bases on REF vs ALT reads (INFO field):
    - Are there more or less hard clips for reference alleles compared to
      alternate alleles. Ideally want close to 0. Negative/Positive = more
      alt/ref allele reads have hard-clipped bases.

* [DepthPerSampleHC](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_DepthPerSampleHC.php) -
Depth of informative coverage for each sample as considedered by HaplotypeCaller (HC) (FORMAT field):
    - Count of reads considered informative by HaplotypeCaller (i.e. allele it
      carries can be easily distinguished, as opposed to a read partially
      overlapping an STR).

* [FisherStrand](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_FisherStrand.php) -
Strand bias estimated using Fisher's Exact Test (INFO field):
    - You can sometimes get a different number of reads from a particular strand
      (forward/reverse) supporting one allele vs. the other. The Fisher's Exact
      Test determines if this difference is statistically significant. Output is
      a Phred-scaled p-value. Higher value means more likely to be bias. More
      bias indicates false positive calls.

* [MappingQualityRankSumTest](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityRankSumTest.php) -
Rank Sum Test for mapping qualities of REF versus ALT reads (INFO field):
    - Compare mapping qualities of REF versus ALT reads. Ideally want close to
      0. Negative/Positive = ALT/REF reads have lower mapping qualities than
      REF/ALT reads. In practice, we only filter out low negative values because
      we want to filter out ALT sites supported by low quality reads.

* [MappingQualityZero](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityZero.php) -
Count of all reads with MAPQ = 0 across all samples (INFO field):
    - Count of how many reads have MAPQ = 0 across all samples. High counts
      indicate regions where it is difficult to make confident calls.

* [QualByDepth](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_QualByDepth.php) -
Variant call confidence normalized by depth of sample reads supporting a variant (INFO field):
    - Since the QUAL score is highly influenced by depth of coverage, you can
      normalise based on the depth and get a more objective picture of how well
      supported the call is.

* [ReadPosRankSumTest](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_ReadPosRankSumTest.php) -
Rank Sum Test for relative positioning of REF versus ALT alleles within reads (INFO field):
    - Check for bias in the read position of REF versus ALT sites.


* [RMSMappingQuality](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_RMSMappingQuality.php) -
Root Mean Square of the mapping quality of reads across all samples (INFO field):
    - Estimate of overall mapping quality of ALT reads.

* MappingQuality (not sure)

* [DepthPerAlleleBySample](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_DepthPerAlleleBySample.php) -
Depth of coverage of each allele per sample (FORMAT field):
    - Unfiltered count of reads that support REF/ALT for an individual sample.
      Order of values is REF, ALT1, ALT2 etc.

* [Coverage](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_Coverage.php) -
Total depth of coverage per sample and over all samples (INFO, FORMAT fields):
    - At the variant site level (INFO), the DP value is the unfiltered depth of
      coverage over all samples.
    - At the sample level (FORMAT), the DP value is the count of reads that
      passed the internal qc metrics of the caller (e.g. MAPQ > 17).


## FilterMutectCalls

```
unset JAVA_HOME &&
export PATH=/data/projects/punim0010/local/share/bcbio/anaconda/bin:$PATH &&
gatk-launch \
--javaOptions '-Xms681m -Xmx3181m -XX:+UseSerialGC -Djava.io.tmpdir=tmpR4gahi' \
FilterMutectCalls \
--variant work/bcbiotx/tmpR4gahi/batch1-21_0_19327180-raw.vcf.gz \
--output work/bcbiotx/tmpR4gahi/batch1-21_0_19327180.vcf.gz
```

* GATK [documentation](https://software.broadinstitute.org/gatk/gatkdocs/4.beta.2/org_broadinstitute_hellbender_tools_walkers_mutect_FilterMutectCalls.php)
* Source code in GATK
[GitHub repo](https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect)

### Filters

For source code, check out the `Mutect2FilteringEngine.java` in the above repo.

* `tumor_lod`: minimum likelihood of an allele as determined by the somatic
  likelihoods model required to pass. `TUMOR_LOD_THRESHOLD = 5.3`.
    * LOD threshold for calling tumor variant. Only variants with tumor LODs
      exceeding the threshold can pass filtering
* `maxEventsInHaplotype`: the maximum allowable number of called variants
  co-occurring in a single assembly region. If the number of called variants
  exceeds this they will all be filtered. Note that this filter is misnamed
  because it counts the total number of events over all haplotypes in an
  assembly region. `maxEventsInHaplotype = 2`.
    * Variants coming from a haplotype with more than this many events are filtered
* `uniqueAltReadCount`: the minimum number of unique (start position, fragment
  length) pairs required to make a call. This count is a proxy for the number
  of unique molecules (as opposed to PCR duplicates) supporting an allele.
  Normally PCR duplicates are marked and filtered by the GATK engine, but in
  UMI-aware calling this may not be the case, hence the need for this filter.
  `uniqueAltReadCount = 0`.
    * Filter a variant if a site contains fewer than this many unique
      (i.e. deduplicated) reads supporting the alternate allele
* `maxAltAllelesThreshold`: the maximum allowable number of alt alleles at a
  site. By default only biallelic variants pass the filter. `numAltAllelesThreshold = 1`
    * Filter variants with too many alt alleles.
* `max_germline_posterior`: the maximum posterior probability, as determined by
  the above germline probability model, that a variant is a germline event.
  `maxGermlinePosterior = 0.025`.
    * Maximum posterior probability that an allele is a germline variant.
      The default value of this argument was chosen to achieve roughly one false
      positive per covered megabase according to a back-of-the-envelope
      calculation assuming the neutral model of allele selection, with
      parameters estimated from ExAC.
* `normal_artifact_lod`: the maximum acceptable likelihood of an allele in the
  normal _by the somatic likelihoods model_. This is different from the normal
  likelihood that goes into the germline model, which makes a diploid assumption.
  Here we compute the normal likelihood as if it were a tumor in order to detect
  artifacts. `NORMAL_ARTIFACT_LOD_THRESHOLD = 0.0`.
    * LOD threshold for calling normal artifacts. Measure of the minimum
      evidence to support that a variant observed in the tumor is not
      also present in the normal as an artifact i.e. not as a germline event.
* `strandArtifactPosteriorProbability`: the posterior probability of a strand
  artifact, as determined by the model described above, required to apply the
  strand artifact filter. This is necessary but not sufficient - we also
  require the estimated max a posteriori allele fraction to be less than
  `strandArtifactAlleleFraction`. The second condition prevents filtering
  real variants that also have significant strand bias, i.e. a true
  variant that _also_ has some artifactual reads. `STRAND_ARTIFACT_POSTERIOR_PROB_THRESHOLD = 0.99`.
    * Filter a variant if the probability of strand artifact exceeds this
      number.
* `strandArtifactAlleleFraction`: `STRAND_ARTIFACT_ALLELE_FRACTION_THRESHOLD = 0.01`.
    * Filter a variant if the MAP estimate of allele fraction given artifact is
      below this number.
* `minMedianBaseQuality`: the minimum median base quality of bases supporting a SNV.
  `minMedianBaseQuality = 20`.
    * Filter variants for which the median base quality of alt reads is too
      low.
* `minMedianMappingQuality`: the minimum median mapping quality of reads supporting an allele.
  `minMedianMappingQuality = 30`.
    * Filter variants for which the median mapping quality of alt reads is too
      low.
* `maxMedianFragmentLengthDifference`: the maximum difference between the median fragment
  lengths reads supporting alt and reference alleles. Note that fragment
  length is based on where paired reads are mapped, not the actual physical fragment length.
  `maxMedianFragmentLengthDifference = 10000`.
    * Filter variants for which the median fragment length for alt reads is very
      different from the median for ref reads.
* `minMedianReadPosition`: the minimum median length of bases supporting an
  allele from the closest end of the read. Indels positions are measured by
  the end farthest from the end of the read. `minMedianReadPosition = 5`.
    * Filter variants for which the median position of alt alleles within reads
      is too near the end of reads.
* `contaminationTable`: if `FilterMutectCalls` is passed a `contaminationTable`
  from `CalculateContamination` it will filter alleles with allele fraction less
  than the whole-bam contamination in the table. `contaminationTable = null`.
    * Table containing contamination information.
