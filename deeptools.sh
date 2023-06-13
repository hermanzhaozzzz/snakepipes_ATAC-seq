computeMatrix reference-point \
-S \
../bigwig/ATACSeq_GFP-NLS_REP-1.RPKM.bw \
../bigwig/ATACSeq_GFP_REP-1.RPKM.bw \
../bigwig/ATACSeq_ND6-DddAwt_REP-1.RPKM.bw \
../bigwig/ATACSeq_ND6-DddAwt_REP-2.RPKM.bw \
../bigwig/ATACSeq_SIRT6-DddA11_REP-1.RPKM.bw \
../bigwig/ATACSeq_SIRT6-DddA11_REP-2.RPKM.bw \
-R \
../bed_for_count/2023-05-12_CTCF_medium_for_count.bed \
../bed_for_count/2023-05-12_CTCF_strong_for_count.bed \
../bed_for_count/2023-05-12_CTCF_weak_for_count.bed \
../bed_for_count/2023-05-12_CTCF_within-offs_for_count.bed \
../bed_for_count/2023-05-12_random_for_count.bed \
../bed_for_count/2023-05-12_tss_has_ctcf_for_count.bed \
../bed_for_count/2023-05-12_tss_no_ctcf_for_count.bed \
../bed_for_count/2023-05-12_offs-ND6-DEP_for_count.bed \
../bed_for_count/2023-05-12_offs-ND6-IND_for_count.bed \
-o ../ATAC-seq.RPKM.out.mat.gz \
--referencePoint center \
--beforeRegionStartLength 2000 \
--afterRegionStartLength 2000 \
--skipZeros \
--binSize 10 \
--samplesLabel \
GFP-NLS_REP-1 \
GFP_REP-1 \
ND6-DddAwt_REP-1 \
ND6-DddAwt_REP-2 \
SIRT6-DddA11_REP-1 \
SIRT6-DddA11_REP-2 \
--numberOfProcessors 22









plotHeatmap -m ../ATAC-seq.RPKM.out.mat.gz \
--perGroup \
-out coverage.heatmap.base-on-region.pdf \
--dpi 200 --heatmapHeight 40 --heatmapWidth 3.5 \
--plotTitle "Coverage distribution of ATAC-seq"  --legendLocation upper-left --zMax 45 \
--legendLocation lower-left \
--regionsLabel \
CTCF_medium \
CTCF_strong \
CTCF_weak \
CTCF_within-offs \
random \
tss_has_ctcf \
tss_no_ctcf \
offs-ND6-DEP \
offs-ND6-IND \
--colorList white,lightblue,purple
