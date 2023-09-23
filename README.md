# NLRP3

Pathway PRS and Rare variant analysis described in detail in relevant branches. 

SMR analysis was performed using Yang Lab software with standart settings (https://yanglab.westlake.edu.cn/)


Pathway PRS pipeline

1. Perform standard QC on your base (GWAS) and target (genotyping) data (reference: https://choishingwan.github.io/PRS-Tutorial/) - script not shown here
2. Create a bed file for genes in the pathways (please make sure their coordinates are consistent with the reference genome of your base and target data) - The_NLRP3_inflammasome.bed
3. Calculate PRS using these genes only with PRSice - step1.Pathway_PRS_main.sh
4. Perform logistic regressions between the pathway PRS and the phenotype for each cohort - step2.combined_logistic_regression.r
5. Meta analyze the results - step3.meta_analysis.r