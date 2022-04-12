source ~/.bashrc_miniconda3
conda activate abc_enhancer_gene
#only need this once for formatting HiC data and calculate powerlaw fit
module load juicer
cd /home/FCAM/jcotney/ANALYSIS/ABC_Enhancer_Gene_Prediciton
#see directories for reference files and enhancer annotations along with scripts for generation
#gene expression outputs from Rail-RNA via Tara
#calculate average expression per timepoint/tissue as necessary
#use chromhmm calls or DHS peaks overlapping active regions
#filter gene expression based on the gene lists we are providing in gencode.v19.annotation.genes.known.sorted.bed
#might also need to further filter for expression level
#including low GINI index scoring genes from our analysis as ubiquitously_expressed_genes
#filtering of DHS peaks for enhancer states

grep -f enhancer_states.txt REDONE_NCC_R1_25state_dense.sorted.bed > REDONE_NCC_R1_25state_dense.sorted.enhancer_states.bed
grep -f enhancer_states.txt REDONE_NCC_R2_25state_dense.sorted.bed > REDONE_NCC_R2_25state_dense.sorted.enhancer_states.bed
multiIntersectBed -i REDONE_NCC_R*_25state_dense.sorted.enhancer_states.bed | awk '{ if ($4 > 1) print $1"\t"$2"\t"$3}' > reproducible_CNCC_enhancer_states.bed
multiIntersectBed -i REDONE_NCC_R*_25state_dense.sorted.bed | awk '{ if ($4 > 1) print $1"\t"$2"\t"$3}' | bedtools merge -i stdin | bedtools intersect -a stdin -b reproducible_CNCC_enhancer_states.bed > CNCC_DNAse_overlap_enhancers.bed

cat CNCC_tpm.txt | awk '{if ($2 >=0.5) print $0}' > CNCC_tpm.filtered.txt


python ~/TOOLS/ABC-Enhancer-Gene-Prediction/src/run.neighborhoods.py \
--candidate_enhancer_regions enhancer_lists/CNCC_DNAse_overlap_enhancers.bed \
--genes gencode.v19.annotation.genes.known.sorted.bed \
--H3K27ac /home/FCAM/jvanoudenhove/ANALYSIS/CHIPSEQ/NCC_CHIPSEQ/Human/ChromImpute/CHROMIMPUTE_APPLY/impute_NCC_R1_H3K27ac.pval.signal.bigWig,/home/FCAM/jvanoudenhove/ANALYSIS/CHIPSEQ/NCC_CHIPSEQ/Human/ChromImpute/CHROMIMPUTE_APPLY/impute_NCC_R2_H3K27ac.pval.signal.bigWig \
--DHS /home/FCAM/jvanoudenhove/ANALYSIS/CHIPSEQ/NCC_CHIPSEQ/Human/ChromImpute/CHROMIMPUTE_APPLY/impute_NCC_R1_DNase.pval.signal.bigWig,/home/FCAM/jvanoudenhove/ANALYSIS/CHIPSEQ/NCC_CHIPSEQ/Human/ChromImpute/CHROMIMPUTE_APPLY/impute_NCC_R2_DNase.pval.signal.bigWig \
--expression_table CNCC_tpm.filtered.txt \
--ubiquitously_expressed_genes gencode_gini_ubiquitous.txt \
--chrom_sizes ~/GENOME/hg19/dna/hg19_nh.no_chrY.no_chrM.chrom.sizes \
--cellType CNCC \
--outdir ABC_output_CNCC/DHS_Filtered_Neighborhoods/


python ~/TOOLS/ABC-Enhancer-Gene-Prediction/src/predict.py \
--enhancers ABC_output_CNCC/DHS_Filtered_Neighborhoods/EnhancerList.txt \
--genes ABC_output_CNCC/DHS_Filtered_Neighborhoods/GeneList.txt \
--HiCdir HiC/juicer/ \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold .03 \
--cellType CNCC \
--outdir ABC_output_CNCC/DHS_Filtred_Predictions_wHiC/ \
--make_all_putative

cat ABC_output_CNCC/DHS_Filtred_Predictions_wHiC/EnhancerPredictions.bedpe | awk '{if ($2 < $6) {print $1"\t"$2"\t"$6"\t"$7"\t"int($8*1000)"\t"$8"\tCS13-12383\t#7A67EE\t"$1"\t"$2"\t"$3"\t"$7"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$10} else {print $4"\t"$6"\t"$2"\t"$7"\t"int($8*1000)"\t"$8"\tCS13-12383\t#7A67EE\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t"$1"\t"$2"\t"$3"\t"$7"\t"$10}}' > ABC_output_CNCC/DHS_Filtred_Predictions_wHiC/EnhancerPredictions.interact
cat ABC_output_CNCC/DHS_Filtred_Predictions_wHiC/EnhancerPredictions.interact | awk '{if ($6 > .05) print $0}' > ABC_output_CNCC/DHS_Filtred_Predictions_wHiC/EnhancerPredictions.05.interact