configfile: 'config.yaml'
shell.prefix("set -euo pipefail; ")

rule all:
    input:
        expand('{bam_dir}/{sample}_chr{chrom}_2pc.bam',
               bam_dir = config['out_dir'],
               sample = config['samples'],
               chrom = config['subset_chrom'])

rule subset_bam:
    input:
        config['bam_dir'] + '{sample}.bam'
    output:
        config['out_dir'] + '{sample}_chr{subset_chrom}.bam'
    log:
        config['log_dir'] + '{sample}_chr{subset_chrom}.bam.log'
    threads:
        1
    shell:
        "samtools view {input} -b {wildcards.subset_chrom} > {output} 2> {log}"

rule index_bam:
    input:
        config['out_dir'] + '{sample}_chr{subset_chrom}.bam'
    output:
        config['out_dir'] + '{sample}_chr{subset_chrom}.bam.bai'
    log:
        config['log_dir'] + '{sample}_chr{subset_chrom}.bam.bai.log'
    shell:
        "samtools index {input} {output} 2> {log}"

rule sample_sam:
    input:
        config['out_dir'] + '{sample}_chr{subset_chrom}.bam'
    output:
        config['out_dir'] + '{sample}_chr{subset_chrom}_2pc.bam'
    log:
        config['log_dir'] + '{sample}_chr{subset_chrom}_2pc.bam.log'
    threads:
        1
    params:
        seed = 123.02
    shell:
        "samtools view -b -s {params.seed} {input} > {output} 2> {log}"
