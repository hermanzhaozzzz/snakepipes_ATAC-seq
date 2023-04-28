# --------------------------------------------------------------->>>>>>>
# pipeline
# --------------------------------------------------------------->>>>>>>
# 1. cutadapte and bowtie2 alignment BAM file 
# 2. BAM file sort and remove duplication 
# 3. Keep only MAPQ30

SAMPLES = [
    'Liu-Untreated-rep1',
    'Liu-Untreated-rep2',
    'Liu-Untreated-rep3',
    'OurATACSeq_1',
    'OurATACSeq_2',
    'OurATACSeq_4',
    'OurATACSeq_5',
    'N6-1310-1',
    'N6-1310-2',
    'N6-DI-1',
    'N6-DI-2',
    'N6-UNES-1',
    'N6-UNES-2',
    'N6-WT-1',
    'N6-WT-2',
]


THREADS = '24'



# --------------------------------------------------------------->>>>>>>
# index and files
# --------------------------------------------------------------->>>>>>>
HG38_BWA_INDEX = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.bwa_index"
HG38_BT2_INDEX = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.bowtie2_index"
HG38_BWA_ChrM_INDEX = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/chrM.fa.bwa_index"
HG38_BT2_ChrM_INDEX = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/chrM.fa.bowtie2_index"
HG38_FA = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa"
HG38_ChrM_FA = "/home/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/chrM.fa"

# --------------------------------------------------------------->>>>>>>
# software
# --------------------------------------------------------------->>>>>>>
PYTHON = "/home/zhaohuanan/zhaohn_HD/miniconda3/envs/snakepipes_target-seq/bin/python" # ok
PICARD = "/home/zhaohuanan/zhaohn_HD/1.apps/picard/picard.jar"
GATK38 = "/home/zhaohuanan/zhaohn_HD/1.apps/GenomeAnalysisTK/GenomeAnalysisTK.jar"
VARSCAN = "/home/zhaohuanan/zhaohn_HD/miniconda3/bin/../share/varscan-2.4.4-0/VarScan.jar"
with os.popen("which cutadapt") as path:
    CUTADAPT = path.read().strip()
    print('PATH cutadapt:', CUTADAPT)
with os.popen("which bowtie2") as path:
    BOWTIE2 = path.read().strip()
    print('PATH bowtie2:', BOWTIE2)
with os.popen("which samtools") as path:
    SAMTOOLS = path.read().strip()
    print('PATH samtools:', SAMTOOLS)
with os.popen("which java") as path:
    JAVA = path.read().strip()
    print('PATH java:', JAVA)
with os.popen("which samtools") as path:
    SAMTOOLS = path.read().strip()



# --------------------------------------------------------------->>>>>>>
# snakepipes
# --------------------------------------------------------------->>>>>>>
rule all:
    input:
#         expand("../fastq/{sample}_combined_R1.fastq.gz", sample=SAMPLES),
#         expand("../fastq/{sample}_combined_R2.fastq.gz", sample=SAMPLES),
        expand("../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort.bam", sample=SAMPLES),
        expand("../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort.bam", sample=SAMPLES),
        expand("../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.bam", sample=SAMPLES),
        expand("../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam", sample=SAMPLES),
        expand("../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort.mapping_stats", sample=SAMPLES),
        expand("../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort.mapping_stats", sample=SAMPLES),
        expand("../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.mapping_stats", sample=SAMPLES),
        expand("../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.mapping_stats", sample=SAMPLES),
        expand("../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.bam.bai", sample=SAMPLES),
        expand("../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam.bai", sample=SAMPLES),
        expand("../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort.bam.bai", sample=SAMPLES),
        expand("../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort.bam.bai", sample=SAMPLES),
        expand("../mpileup/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.mpileup", sample=SAMPLES),
        expand("../mpileup/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.mpileup",sample=SAMPLES),
        expand("../bmat/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.bmat",sample=SAMPLES),
        expand("../bmat/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.bmat",sample=SAMPLES),
        expand("../vcf/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.vcf", sample=SAMPLES),
        expand("../vcf/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.vcf", sample=SAMPLES)
        
rule bowtie2_mapping:
    input:
        "../fastq/{sample}_combined_R1.fastq.gz",
        "../fastq/{sample}_combined_R2.fastq.gz"
    output:
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.sam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.sam"
    log:
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.sam.log",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.sam.log",
    params:
        tag_RG = "{sample}",
        tag_ID = "ID:{sample}",
        tag_SM = "SM:{sample}",
        tag_PL = "PL:ILLUMINA"
    shell:
        """
        {BOWTIE2} -x {HG38_BT2_INDEX} --rg-id {params.tag_RG} --rg {params.tag_ID} --rg {params.tag_SM} --rg {params.tag_PL} -p {THREADS} -1 {input[0]} -2 {input[1]} -S {output[0]} > {log[0]} 2>&1
        {BOWTIE2} -x {HG38_BT2_ChrM_INDEX} --rg-id {params.tag_RG} --rg {params.tag_ID} --rg {params.tag_SM} --rg {params.tag_PL} -p {THREADS} -1 {input[0]} -2 {input[1]} -S {output[1]} > {log[1]} 2>&1
        """    

rule filter_MAPQ30:
    input:
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.sam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.sam"
    output:
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30.bam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30.bam"
    shell:
        """
        {SAMTOOLS} view -hb -q 30 -F 4 -F 8 -o {output[0]} {input[0]}
        {SAMTOOLS} view -hb -q 30 -F 4 -F 8 -o {output[1]} {input[1]}
        """

rule samtools_sort_by_position:
    input:
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30.bam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30.bam"
    output:
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort.bam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort.bam"
    shell:
        """
        {SAMTOOLS} sort -O BAM -o {output[0]} -T {output[0]}.temp -@ {THREADS} -m 4G {input[0]}
        {SAMTOOLS} sort -O BAM -o {output[1]} -T {output[1]}.temp -@ {THREADS} -m 4G {input[1]}
        """

rule mark_duplicate:
    input:
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort.bam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort.bam"
    output:
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.bam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.matrix",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.matrix"
    log:
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.log",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.log"

    shell:
        """
        {JAVA} -Xms80g -Xmx80g -XX:ParallelGCThreads={THREADS} -jar {PICARD} MarkDuplicates I={input[0]} O={output[0]} M={output[1]} ASO=coordinate REMOVE_DUPLICATES=true 2>{log[0]}
        {JAVA} -Xms80g -Xmx80g -XX:ParallelGCThreads={THREADS} -jar {PICARD} MarkDuplicates I={input[1]} O={output[2]} M={output[3]} ASO=coordinate REMOVE_DUPLICATES=true 2>{log[1]}
        """


rule samtools_index:
    input:
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.bam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort.bam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort.bam"
    output:
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.bam.bai",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam.bai",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort.bam.bai",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort.bam.bai"

    shell:
        """
        {SAMTOOLS} index -@ {THREADS} {input[0]} {output[0]} \\
        
        {SAMTOOLS} index -@ {THREADS} {input[1]} {output[1]} \\
        
        {SAMTOOLS} index -@ {THREADS} {input[2]} {output[2]} \\
        
        {SAMTOOLS} index -@ {THREADS} {input[3]} {output[3]}
        """        
rule get_mapping_stats_rmdup:
    input:
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort.bam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort.bam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.bam",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam"
    output:
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort.mapping_stats",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort.mapping_stats",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.mapping_stats",
        "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.mapping_stats"
    shell:
        """
        {SAMTOOLS} stats -@ {THREADS} --remove-overlaps --reference {HG38_FA} {input[0]} > {output[0]}
        {SAMTOOLS} stats -@ {THREADS} --remove-overlaps --reference {HG38_ChrM_FA} {input[1]} > {output[1]}
        {SAMTOOLS} stats -@ {THREADS} --remove-overlaps --reference {HG38_FA} {input[2]} > {output[2]}
        {SAMTOOLS} stats -@ {THREADS} --remove-overlaps --reference {HG38_ChrM_FA} {input[3]} > {output[3]}
        """
rule bam_mpileup_genome:
    input:
        bam = "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.bam",
        bam_index = "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.bam.bai"
    output:
        "../mpileup/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.mpileup"
    shell:
        "samtools mpileup --reference {HG38_FA} -q 30 -Q 30 -o {output} {input.bam}"
rule bam_mpileup_chrM:
    input:
        bam = "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam",
        bam_index = "../bam/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam.bai"
    output:
        "../mpileup/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.mpileup"
    shell:
        "{SAMTOOLS} mpileup --reference {HG38_ChrM_FA} -q 30 -Q 30 -o {output} {input.bam}"

rule mpileup_to_bmat_genome:
    input:
        "../mpileup/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.mpileup"
    output:
        "../bmat/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.bmat"
    log:
        "../bmat/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.bmat.log"
    shell:
        "{PYTHON} program/parse-mpileup-V04.py -i {input} -o {output} -p 24 -n 0 > {log} 2>&1"
rule mpileup_to_bmat_chrM:
    input:
        "../mpileup/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.mpileup"
    output:
        "../bmat/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.bmat"
    log:
        "../bmat/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.bmat.log"
    shell:
        "{PYTHON} program/parse-mpileup-V04.py -i {input} -o {output} -p 24 -n 0 > {log} 2>&1"
rule varscan_genome:
    input:
        "../mpileup/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.mpileup"
    output:
        "../vcf/293T-DdCBE-ATACSeq-{sample}_bt2_hg38.MAPQ30_sort_rmdup.vcf"
    shell:
        """
        {JAVA} -Xms80g -Xmx80g -XX:ParallelGCThreads={THREADS} -jar {VARSCAN} mpileup2snp \
        {input} \
        --min-coverage 8 \
        --min-reads2 2 \
        --min-avg-qual 30 \
        --min-var-freq 0.001 \
        --min-freq-for-hom 0.75 \
        --p-value 0.95 \
        --strand-filter 1 \
        --output-vcf 1 > {output}
        """
rule varscan_chrM:
    input:
        "../mpileup/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.mpileup"
    output:
        "../vcf/293T-DdCBE-ATACSeq-{sample}_bt2_hg38_chrM.MAPQ30_sort_rmdup.vcf"
    shell:
        """
        {JAVA} -Xms80g -Xmx80g -XX:ParallelGCThreads={THREADS} -jar {VARSCAN} mpileup2snp \
        {input} \
        --min-coverage 8 \
        --min-reads2 2 \
        --min-avg-qual 30 \
        --min-var-freq 0.001 \
        --min-freq-for-hom 0.75 \
        --p-value 0.95 \
        --strand-filter 1 \
        --output-vcf 1 > {output}
        """