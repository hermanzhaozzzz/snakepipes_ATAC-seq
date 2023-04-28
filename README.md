# snakepipes_atac-seq

## step01
run step01.ATAC-seq_mapping_to_hg38_and_chrM.py

## step02
```
cd ../bam

ls *hg38.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-1310-1_bt2_hg38.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-1310-2_bt2_hg38.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-DI-1_bt2_hg38.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-DI-2_bt2_hg38.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-UNES-1_bt2_hg38.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-UNES-2_bt2_hg38.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-WT-1_bt2_hg38.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-WT-2_bt2_hg38.MAPQ30_sort.bam


ls *hg38.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-1310-1_bt2_hg38.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-1310-2_bt2_hg38.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-DI-1_bt2_hg38.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-DI-2_bt2_hg38.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-UNES-1_bt2_hg38.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-UNES-2_bt2_hg38.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-WT-1_bt2_hg38.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-WT-2_bt2_hg38.MAPQ30_sort_rmdup.bam

ls *hg38_chrM.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-1310-1_bt2_hg38_chrM.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-1310-2_bt2_hg38_chrM.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-DI-1_bt2_hg38_chrM.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-DI-2_bt2_hg38_chrM.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-UNES-1_bt2_hg38_chrM.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-UNES-2_bt2_hg38_chrM.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-WT-1_bt2_hg38_chrM.MAPQ30_sort.bam
293T-DdCBE-ATACSeq-N6-WT-2_bt2_hg38_chrM.MAPQ30_sort.bam

ls *hg38_chrM.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-1310-1_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-1310-2_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-DI-1_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-DI-2_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-UNES-1_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-UNES-2_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-WT-1_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam
293T-DdCBE-ATACSeq-N6-WT-2_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam
```

```
samtools idxstats 293T-DdCBE-ATACSeq-N6-1310-1_bt2_hg38.MAPQ30_sort.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-1310-2_bt2_hg38.MAPQ30_sort.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-DI-1_bt2_hg38.MAPQ30_sort.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-DI-2_bt2_hg38.MAPQ30_sort.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-UNES-1_bt2_hg38.MAPQ30_sort.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-UNES-2_bt2_hg38.MAPQ30_sort.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-WT-1_bt2_hg38.MAPQ30_sort.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-WT-2_bt2_hg38.MAPQ30_sort.bam  | awk -F "\t" '{print $3}'


samtools idxstats 293T-DdCBE-ATACSeq-N6-1310-1_bt2_hg38.MAPQ30_sort_rmdup.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-1310-2_bt2_hg38.MAPQ30_sort_rmdup.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-DI-1_bt2_hg38.MAPQ30_sort_rmdup.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-DI-2_bt2_hg38.MAPQ30_sort_rmdup.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-UNES-1_bt2_hg38.MAPQ30_sort_rmdup.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-UNES-2_bt2_hg38.MAPQ30_sort_rmdup.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-WT-1_bt2_hg38.MAPQ30_sort_rmdup.bam  | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-WT-2_bt2_hg38.MAPQ30_sort_rmdup.bam  | awk -F "\t" '{print $3}'


samtools idxstats 293T-DdCBE-ATACSeq-N6-1310-1_bt2_hg38_chrM.MAPQ30_sort.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-1310-2_bt2_hg38_chrM.MAPQ30_sort.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-DI-1_bt2_hg38_chrM.MAPQ30_sort.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-DI-2_bt2_hg38_chrM.MAPQ30_sort.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-UNES-1_bt2_hg38_chrM.MAPQ30_sort.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-UNES-2_bt2_hg38_chrM.MAPQ30_sort.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-WT-1_bt2_hg38_chrM.MAPQ30_sort.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-WT-2_bt2_hg38_chrM.MAPQ30_sort.bam | awk -F "\t" '{print $3}'

samtools idxstats 293T-DdCBE-ATACSeq-N6-1310-1_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-1310-2_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-DI-1_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-DI-2_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-UNES-1_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-UNES-2_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-WT-1_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam | awk -F "\t" '{print $3}'
samtools idxstats 293T-DdCBE-ATACSeq-N6-WT-2_bt2_hg38_chrM.MAPQ30_sort_rmdup.bam | awk -F "\t" '{print $3}'
```
