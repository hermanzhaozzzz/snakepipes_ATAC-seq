BAMs = ['GFP-NLS_REP-1', 'GFP_REP-1', 'ND6-DddAwt_REP-1', 'ND6-DddAwt_REP-2', 'SIRT6-DddA11_REP-1', 'SIRT6-DddA11_REP-2']
BEDs = [
    # 'CTCF_medium', 'CTCF_strong', 'CTCF_weak', 'CTCF_within-offs', 'offs-ND6-DEP', 'offs-ND6-IND', 'random', 
    'tss_has_ctcf', 'tss_no_ctcf'
]



rule all:
    input:
        expand('../bed_for_count/AnalysisChanges_BAM_{bam}_BED_{bed}.count_table', bam=BAMs, bed=BEDs)
rule count_table:
    input:
        bed='../bed_for_count/2023-05-12_{bed}_for_count.bed',
        bam='../bam/ATACSeq_{bam}.sortp_genome_rmchrM_rmdup.bam'
    output:
        '../bed_for_count/AnalysisChanges_BAM_{bam}_BED_{bed}.count_table'
    shell:
        """
        bedtools coverage \
            -a {input.bed} -b {input.bam} > {output}
        """

