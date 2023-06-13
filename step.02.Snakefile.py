# ——————————————————>>>>>>>>>>
# Project Name: ATAC-seq standard protocol
# Author: Hua-nan ZHAO
# E-mail: hermanzhaozzzz@gmail.com
# Update log:
#     2023-04-28: start project
# ——————————————————>>>>>>>>>>
import os
import json
# ------------------------------------------------------------------->>>>>>>>>>
# FUNCTIONS
# ------------------------------------------------------------------->>>>>>>>>>
def print_head(SAMPLES, MODE):
    print('----------\nSAMPLES:')
    [print('\t' + i) for i in SAMPLES]
    print('----------\nMODE:')
    print('\t' + MODE)
    print('----------\n\n')

def check_cmd(x):
    return any(
        os.access(os.path.join(path, x), os.X_OK) 
        for path in os.environ["PATH"].split(os.pathsep)
    )

def check_read(x):
    if x == "PE":
        read = ['R1', 'R2']
    elif x == "SE":
        read = ['SE']
    else:
        raise ValueError()
    return read
# ------------------------------------------------------------------->>>>>>>>>>
# SAMPLE INFO
# ------------------------------------------------------------------->>>>>>>>>>
with open('./samples.json') as f:
    dt = json.loads(f.read())

SAMPLES = dt['samples']
MODE = dt['seq_mode']
THREAD = dt['thread']
R1_AD = dt["r1_ad"]
R2_AD = dt["r2_ad"]
MACS2_GSIZE = dt["macs2_gsize"]


READ = check_read(MODE)

print_head(SAMPLES, MODE)
print(READ)
# ------------------------------------------------------------------->>>>>>>>>>
# DATABASE INFO
# ------------------------------------------------------------------->>>>>>>>>>
GENOME = dt['genome']
BOWTIE2_INDEX_GENOME = dt['bowtie2_index_genome']
BOWTIE2_INDEX_chrM = dt['bowtie2_index_chrM']
BOWTIE2_INDEX_PLASMID = dt['bowtie2_index_plasmid']
# ------------------------------------------------------------------->>>>>>>>>>
# SOFTWARE INFO
# ------------------------------------------------------------------->>>>>>>>>>
# check if cmd exists
assert check_cmd("fastp")
assert check_cmd("bowtie2")
assert check_cmd("bedtools")
assert check_cmd("samtools")
assert check_cmd("samblaster")
assert check_cmd("gatk")
assert check_cmd("macs2")
assert check_cmd("bamCoverage")
assert check_cmd("computeMatrix")
assert check_cmd("plotHeatmap")
FASTP = "fastp"
BOWTIE2 = 'bowtie2'
BEDTOOLS = 'bedtools'
SAMTOOLS = 'samtools'
SAMBLASTER = 'samblaster'
MACS2 = 'macs2'
GATK4 = 'gatk'
BAM_COVERAGE = 'bamCoverage'
COMPUTE_MATRIX = 'computeMatrix'
PLOT_HEATMAP = 'plotHeatmap'
# ------------------------------------------------------------------->>>>>>>>>>
# rule all
# ------------------------------------------------------------------->>>>>>>>>>
rule all:
    input:
        expand("../bam/{sample}.sortp_chrM.bam.flagstat", sample=SAMPLES),
        expand("../bam/{sample}.sortn_genome_rmchrM.bam", sample=SAMPLES),
        expand("../bam/{sample}.sortn_genome_rmchrM.bam.flagstat", sample=SAMPLES),
        expand("../bam/{sample}.sortp_genome_rmchrM_rmdup.bam", sample=SAMPLES),
        expand("../bam/{sample}.sortp_genome_rmchrM_rmdup.bam.flagstat", sample=SAMPLES),
        expand("../plot/atac-seq_insert-size_{sample}.pdf", sample=SAMPLES),
        expand("../bam/{sample}.sortn_final_ummapped_try_plasmid.bam.flagstat", sample=SAMPLES),
        expand("../count_table/{sample}.TSS_2Kbp.count_table", sample=SAMPLES),
        expand("../bam/{sample}.sortp_genome_rmchrM_rmdup.reheader.bam", sample=SAMPLES),
        expand("../bigwig/{sample}.RPKM.bw", sample=SAMPLES),

# ------------------------------------------------------------------->>>>>>>>>>
# trim adaptor
# ------------------------------------------------------------------->>>>>>>>>>
rule fastp_trim_adaptor:
    input:
        fwd="../fastq/{sample}_R1.fastq.gz",
        rev="../fastq/{sample}_R2.fastq.gz"
    output:
        fwd=temp("../fastq_trim_adaptor/{sample}_R1_cutadapt.fq.gz"),
        rev=temp("../fastq_trim_adaptor/{sample}_R2_cutadapt.fq.gz"),
        html="../fastq_trim_adaptor/{sample}.html",
        json="../fastq_trim_adaptor/{sample}.json"
    log:
        "../fastq_trim_adaptor/{sample}.log"
    message:
        "trim_adaptor {input}"
    shell:
        """
        {FASTP} -w {THREAD} -h {output[html]} -j  {output[json]} \
            -i {input[fwd]} -I {input[rev]} -o {output[fwd]} \
            -O {output[rev]} 2> {log}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# bowtie2 mapping
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# BOWTIE2_INDEX_chrM
rule bowtie2_align_chrM:
    input:  # TODO改参数1000 2000  先mapping线粒体
        "../fastq_trim_adaptor/{sample}_R1_cutadapt.fq.gz",
        "../fastq_trim_adaptor/{sample}_R2_cutadapt.fq.gz"
    output:
        temp("../bam/{sample}.sortn_chrM.bam")
    log:
        "../bam/{sample}.sortn_chrM.bam.bowtie2.log"
    shell:
        """
        {BOWTIE2} --threads {THREAD} --very-sensitive -x {BOWTIE2_INDEX_chrM} -X 2000 \
            -1 {input[0]} -2 {input[1]} 2> {log} \
            | {SAMTOOLS} view -@ {THREAD} -Sb - > {output[0]}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# samtools sort by position
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule BAM_sort_by_position:
    input:
        "../bam/{sample}.sortn_chrM.bam"
    output:
        temp("../bam/{sample}.sortp_chrM.bam")
    shell:
        """
        {SAMTOOLS} sort -O BAM -o {output} -@ {THREAD} {input}
        """
rule flagstat_first_map_to_chrM:
    input:
        "../bam/{sample}.sortp_chrM.bam"
    output:
        "../bam/{sample}.sortp_chrM.bam.flagstat"
    shell:
        """
        {SAMTOOLS} flagstat {input} -@ {THREAD} > {output}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# pick_unmapped_reads step1 unmapped bam
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule pick_unmapped_reads1:
    input:
        "../bam/{sample}.sortp_chrM.bam"
    output:
        temp("../bam/{sample}.sortp_chrM.unmapped.bam")
    shell:
        """
        {SAMTOOLS} view -@ {THREAD} -b -f 4 {input} > {output}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# pick_unmapped_reads step2 unmapped bam sort by name
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule pick_unmapped_reads2:
    input:
        "../bam/{sample}.sortp_chrM.unmapped.bam"
    output:
        temp("../bam/{sample}.sortp_chrM.unmapped_sortn.bam")
    shell:
        """
        {SAMTOOLS} sort -O BAM -o {output} -n -@ {THREAD} {input}
        """     
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# pick_unmapped_reads step3 bam to fastq
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# bamtofastq  bedtools  [OPTIONS] -i <BAM -fq <FASTQ>
# 参数:
# Opti      Description      on
# fq2      默认-fq输出的是一个fastq文件，添加-fq2参数可以将成对的fastq文件分别输出到两个文件中。但是输)      文件需要先对reads按名字进行排序(samtools sort -n aln.bam aln.qsort)
# bedtools bamtofastq -i x.bam -fq /dev/stdout -fq2
rule pick_unmapped_reads3:
    input:
        "../bam/{sample}.sortp_chrM.unmapped_sortn.bam"
    output:
        fq1 = temp("../fastq_unmapped/{sample}_R1.fastq"),
        fq2 = temp("../fastq_unmapped/{sample}_R2.fastq")
    log:
        "../fastq_unmapped/{sample}.log"
    shell:
        """
        {BEDTOOLS} bamtofastq -i {input} -fq {output.fq1} -fq2 {output.fq2}  > {log} 2>&1
        """
rule pick_unmapped_reads4:
    input:
        fq1 = "../fastq_unmapped/{sample}_R1.fastq",
        fq2 = "../fastq_unmapped/{sample}_R2.fastq"
    output:
        fq1 = "../fastq_unmapped/{sample}_R1.fastq.gz",
        fq2 = "../fastq_unmapped/{sample}_R2.fastq.gz"
    shell:
        """
        pigz -p {THREAD} {input.fq1}
        pigz -p {THREAD} {input.fq2}
        """      
## the later step will remove chrM from the bam file and coordinate sort the bam
## so I did not cooridnate sort the bam at this step to save some time.
# -X/--maxins <int>maximum fragment length (500)
# 考虑到基本的ATAC-seq包含了最长2000左右的reads，所以bowtie2的X参数需要修改
rule bowtie2_align_genome:
    input:  # -X 2000  先mapping线粒体
        "../fastq_unmapped/{sample}_R1.fastq.gz",
        "../fastq_unmapped/{sample}_R2.fastq.gz"
    output:
        "../bam/{sample}.sortn_genome_rmchrM.bam"
    log:
        bowtie2 = "../bam/{sample}.sortn_genome_rmchrM.bam.bowtie2.log",
        markdup = "../bam/{sample}.sortn_genome_rmchrM.bam.markdup.log"
    shell:  # samblaster mark duplicates for read id grouped reads. I do not coordinate sort the bam
        """
        {BOWTIE2} --threads {THREAD} --very-sensitive -x {BOWTIE2_INDEX_GENOME} -X 2000 \
            -1 {input[0]} -2 {input[1]} 2> {log.bowtie2} \
            | {SAMBLASTER} 2> {log.markdup} \
            | {SAMTOOLS} view -@ {THREAD} -Sb - > {output[0]}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# check number of reads mapped by samtools flagstat
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule flagstat_second_map_to_genome:
    input:
        "../bam/{sample}.sortn_genome_rmchrM.bam"
    output:
        "../bam/{sample}.sortn_genome_rmchrM.bam.flagstat"
    shell:
        """
        {SAMTOOLS} flagstat {input} -@ {THREAD} > {output}
        """
## shifting the reads are only critical for TF footprint, for peak calling and making bigwigs, it should be fine using the bams without shifting
# https://sites.google.com/site/atacseqpublic/atac-seq-analysis-methods/offsetmethods
rule remove_duplicates:
    input:
        "../bam/{sample}.sortn_genome_rmchrM.bam"
    output:
        "../bam/{sample}.sortp_genome_rmchrM_rmdup.bam",
        "../bam/{sample}.sortp_genome_rmchrM_rmdup.bam.bai"
    shell:  # remove duplicates and reads on chrM, coordinate sort the bam; samblaster expects name sorted bamq
        """
        {SAMTOOLS} view -h {input} \
            | {SAMBLASTER} --removeDups \
            | {SAMTOOLS} view -Sb -F 4 - \
            | {SAMTOOLS} sort -@ {THREAD} -T ../temp_file/{input}.tmp -o {output[0]}
        
        {SAMTOOLS} index {output[0]}
        """
rule flagstat_second_map_to_genome_rmdup:
    input:
        "../bam/{sample}.sortp_genome_rmchrM_rmdup.bam"
    output:
        "../bam/{sample}.sortp_genome_rmchrM_rmdup.bam.flagstat"
    shell:
        """
        {SAMTOOLS} flagstat {input} -@ {THREAD} > {output}
        """
rule plot_insert_size_distribution:
    input:
        "../bam/{sample}.sortp_genome_rmchrM_rmdup.bam"
    output:
        txt = "../plot/atac-seq_insert-size_{sample}.txt",
        pdf = "../plot/atac-seq_insert-size_{sample}.pdf",
        reheader_bam = "../bam/{sample}.sortp_genome_rmchrM_rmdup.reheader.bam"
    shell:
        """
        {SAMTOOLS} reheader -c 'grep -v ^@PG' {input} > {output.reheader_bam}
        {SAMTOOLS} index {output.reheader_bam}
        {GATK4} CollectInsertSizeMetrics \
            --INPUT {output.reheader_bam} \
            --OUTPUT {output.txt} \
            --Histogram_FILE {output.pdf} \
            --METRIC_ACCUMULATION_LEVEL ALL_READS
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# pick_unmapped_reads step1 unmapped bam
# ------------------------------------------------------------------------------------------>>>>>>>>>>   
rule pick_final_unmapped_reads1:
    input:
        "../bam/{sample}.sortn_genome_rmchrM.bam"
    output:
        temp("../bam/{sample}.sortn_final_unmapped.bam")
    shell:
        """
        {SAMTOOLS} view -@ {THREAD} -b -f 4 {input} > {output}
        """ 
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# pick_unmapped_reads step2 bam to fastq
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# bamtofastq  bedtools  [OPTIONS] -i <BAM -fq <FASTQ>
# 参数:
# Opti      Description      on
# fq2      默认-fq输出的是一个fastq文件，添加-fq2参数可以将成对的fastq文件分别输出到两个文件中。但是输)      文件需要先对reads按名字进行排序(samtools sort -n aln.bam aln.qsort)
# bedtools bamtofastq -i x.bam -fq /dev/stdout -fq2
rule pick_final_unmapped_reads2:
    input:
        "../bam/{sample}.sortn_final_unmapped.bam"
    output:
        fq1 = temp("../fastq_unmapped2/{sample}_final_R1.fastq"),
        fq2 = temp("../fastq_unmapped2/{sample}_final_R2.fastq")
    log:
        "../fastq_unmapped2/{sample}_final.log"
    shell:
        """
        {BEDTOOLS} bamtofastq -i {input} -fq {output.fq1} -fq2 {output.fq2}  > {log} 2>&1
        """
rule pick_final_unmapped_reads3:
    input:
        fq1 = "../fastq_unmapped2/{sample}_final_R1.fastq",
        fq2 = "../fastq_unmapped2/{sample}_final_R2.fastq"
    output:
        fq1 = "../fastq_unmapped2/{sample}_final_R1.fastq.gz",
        fq2 = "../fastq_unmapped2/{sample}_final_R2.fastq.gz"
    shell:
        """
        pigz -p {THREAD} {input.fq1}
        pigz -p {THREAD} {input.fq2}
        """      
rule bowtie2_align_plasmid:
    input:  # -X 2000  先mapping线粒体
        "../fastq_unmapped2/{sample}_final_R1.fastq.gz",
        "../fastq_unmapped2/{sample}_final_R2.fastq.gz"
    output:
        "../bam/{sample}.sortn_final_ummapped_try_plasmid.bam"
    log:
        "../bam/{sample}.sortn_genome_rmchrM.bam.bowtie2.log"
    shell:  # samblaster mark duplicates for read id grouped reads. I do not coordinate sort the bam
        """
        {BOWTIE2} --threads {THREAD} --very-sensitive -x {BOWTIE2_INDEX_PLASMID} -X 2000 \
            -1 {input[0]} -2 {input[1]} 2> {log} \
            | {SAMTOOLS} view -@ {THREAD} -Sb - > {output[0]}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# check number of reads mapped by samtools flagstat
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule flagstat_third_map_to_plasmid:
    input:
        "../bam/{sample}.sortn_final_ummapped_try_plasmid.bam"
    output:
        "../bam/{sample}.sortn_final_ummapped_try_plasmid.bam.flagstat"
    shell:
        """
        {SAMTOOLS} flagstat {input} -@ {THREAD} > {output}
        """
# ------------------------------------------------------------------------------------------>>>>>>>>>>
# tss score
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule bedtools_coverage:
    input:
        bam = "../bam/{sample}.sortp_genome_rmchrM_rmdup.bam",
        bed = "ref_data/hg38_ucsc_refseq_NM_TSS_extend2kb.sorted.bed"
    output:
        "../count_table/{sample}.TSS_2Kbp.count_table"
    shell:
        """
        {BEDTOOLS} coverage \
            -a {input.bed} \
            -b {input.bam} \
            -sorted \
            -g ref_data/genome_ucsc_hg38.fa.fai \
            -counts > {output}
        """
rule bam2bigwig_RPKM:
    input:
        "../bam/{sample}.sortp_genome_rmchrM_rmdup.reheader.bam"
    output:
        "../bigwig/{sample}.RPKM.bw"
    shell:
        """
        {BAM_COVERAGE} --bam {input} -o {output} -of bigwig --scaleFactor 1 --binSize 10 -p {THREAD} --normalizeUsing RPKM
        """

# # https://github.com/taoliu/MACS/issues/145
# rule call_peaks_macs2:
#     input: 
#         "../bam/{sample}.sortp_rmchrM.bam"
#     output: 
#         bed = "../macs2_callpeak/{sample}_macs2_peaks.broadPeak"
#     params:
#         name = lambda wildcards: wildcards['sample']
#     log: 
#         "../macs2_callpeak/{sample}_macs2_"
#     shell:  # for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
#         """
#         {MACS2} callpeak -t {input[0]} \
#             --keep-dup all -f BAMPE -g {MACS2_GSIZE} -B \
#             --outdir 06peak_macs2 -n {params.name} -p {config[macs2_pvalue]} \
#             --broad --broad-cutoff {config[macs2_pvalue_broad]} &> {log}
#         """
