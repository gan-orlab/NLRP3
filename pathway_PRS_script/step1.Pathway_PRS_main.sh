#!/bin/bash
# step 1
## this script file simply runs the pathway_prs_by_bed.sh with different arguments, generating pathway PRS
## for all cohorts

name="The_NLRP3_inflammasome.bed:The_NLRP3_inflammasome"
output_folder=/home/liulang/runs/lang/pathway_prs_konstantin/results/;

###PPMI
covariate=PPMI_covar.cov;
pheno=PPMI_covar.pheno;
target=PPMI_chr#_rsid_reformat.QC;
cohort=PPMI;
##main
bash pathway_prs_by_bed.sh $covariate $pheno $target $name $output_folder $cohort;

###APDGC
covariate=APDGC_covar.cov;
pheno=APDGC_covar.pheno;
target=APDGC_chr#_rsid_reformat.QC;
cohort=APDGC
##main
bash pathway_prs_by_bed.sh $covariate $pheno $target $name $output_folder $cohort;

###IPDGC
covariate=IPDGC_covar.cov;
pheno=IPDGC_covar.pheno;
target=IPDGC_chr#_rsid_reformat.QC;
cohort=IPDGC;
##main
bash pathway_prs_by_bed.sh $covariate $pheno $target $name $output_folder $cohort;

###NIND
covariate=NIND_covar.cov;
pheno=NIND_covar.pheno;
target=NIND_chr#_rsid_reformat.QC;
cohort=NIND;
##main
bash pathway_prs_by_bed.sh $covariate $pheno $target $name $output_folder $cohort;

###NGRC
covariate=NGRC_covar.cov;
pheno=NGRC_covar.pheno;
target=NGRC_chr#_rsid_reformat.QC;
cohort=NGRC;
##main
bash pathway_prs_by_bed.sh $covariate $pheno $target $name $output_folder $cohort;

###omnix
covariate=omnix_PD_covar_nomissing.cov;
pheno=omnix_PD_covar_nomissing.pheno;
target=all_omnix_chr#_rsid_reformat.QC;
cohort=omnix;
##main
bash pathway_prs_by_bed.sh $covariate $pheno $target $name $output_folder $cohort;

###UKBB
covariate=UKBB_pheno_covariate.cov
pheno=UKBB_pheno_covariate.pheno
target=ukb22828_c#_b0_v3.QC
cohort=UKBB
bash pathway_prs_by_bed.sh $covariate $pheno $target $name $output_folder $cohort

