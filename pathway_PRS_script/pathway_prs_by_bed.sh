#!/bin/bash
##r version: 4.0.2
##PRSice version 2.3.5
##substitute /path/to/PD/GWAS to PD GWAS file
##substitute /path/to/PRSice/folder to your PRSice directory
##covariate file is a tsv file with ID, age, sex and first 10 PCs
##pheno file is a tsv file with ID and phenotype (the phenotype column is named "phenotype")
##target file is a QCed plink file
##name is in the format of "bedfile:the pathway name"
##output_folder is the directory to the folder
covariate=$1
pheno=$2
target=$3
name=$4
output_folder=$5
cohort=$6
output_name=${output_folder}/${cohort}
Rscript /path/to/PRSice/folder/PRSice.R \
    --prsice /path/to/PRSice/folder/PRSice_linux \
    --a1 A1 \
    --a2 A2 \
    --base /path/to/PD/GWAS \
    --beta \
    --bar-levels 1 \
    --fastscore \
    --ignore-fid \
    --binary-target T \
    --clump-kb 250kb \
    --print-snp  \
    --clump-p 0.05 \
    --clump-r2 0.100000 \
    --bed ${name} \
    --thread 10 \
    --perm 10000 \
    --prevalence 0.005 \
    --cov $covariate \
    --pheno $pheno \
    --pheno-col phenotype \
    --out $output_name \
    --pvalue p \
    --score avg \
    --snp SNP \
    --stat b \
    --target $target \
