#!/usr/bin/env Rscript
#r version 4.1.0
## this script was designed to perform logistic regression for pathway PRS and convert
## the coefficients into odd ratio for multiple pathways and 
## multiple cohorts, putting them all together in preparation of a meta analysis

# scroll down for the main program

# Load the required packages
library(data.table)
library(dplyr)
library(magrittr)
library(broom)
# logistic regression
lr <- function(prs,data,pathway){
    prs <- data.frame(fread(prs))
    data <- data.frame(fread(data))
    data.col = colnames(data)
    # Extract the covariates
    cov_name = data.col[!(data.col %in% c("sample.id","phenotype"))]
    # Extract the phenotype and covariate
    cov_phenotype = data[,c("sample.id",cov_name,"phenotype")]
    #merge data 
    prs = prs[,c("IID",pathway)]#only IID and last column are selected
    ##convert to zscore!!!!!!
    prs[,dim(prs)[2]] = (unlist(prs[,dim(prs)[2]]) - mean(unlist(prs[,dim(prs)[2]])))/sd(unlist(prs[,dim(prs)[2]]))
    
    df = merge(prs,cov_phenotype,by.x = "IID",by.y = "sample.id")
    print(dim(df))
    prs_column_name = colnames(df)[2]
    variables_name = c(prs_column_name,cov_name)
    pheno_name="phenotype"
    #change phenotype if they are either 1 or 2
    if (2 %in% df[,pheno_name]){
    df[,pheno_name] = ifelse(df[,pheno_name] == 2, 1, ifelse(df[,pheno_name] == 1, 0, NA))
    }

    # Fit a logistic regression model with the specified covariates
    formula <- as.formula(paste(paste("phenotype","~"), paste(variables_name, collapse="+")))
    print(formula)
    model <- glm(formula, data = df, family = "binomial")
    return(model)
}

# summarize results for one pathway
main_log_or <- function(df_dir,pathway,output_folder){
    prs.result <- NULL
    for (i in 1:dim(df_dir)[1]){
        cohort = df_dir[i,"cohort"]
        data = df_dir[i,"data"]
        prs = df_dir[i,"prs"]
        # logistic regression
        model = lr(prs,data,pathway)
        model.df = tidy(model) # summarize model statistics into a dataframe
        coef = model.df[model.df$term == pathway,]
        or = coef["estimate"]
        se = coef["std.error"]
        pval = coef["p.value"]
        prs.result %<>% rbind(.,
        data.frame(cohort=cohort,
                    P=as.numeric(pval), 
                    logOR=as.numeric(or),
                    SE=as.numeric(se)))
    }
    write.csv(prs.result,paste(paste(output_folder,pathway,sep=""),"_correlation_analysis(OR).csv",sep = ""),row.names=FALSE)
}

# perform logistic regression on multiple pathways
main_with_multiple_pathway_logOR <- function(df_dir,pathways_lst,output_folder){
    for (pathway in pathways_lst){
        main_log_or(df_dir,pathway,output_folder)
    }
}


# main program
## .best files are pathway PRS calculated from last step
## data file here contains both phenotype and the covariates
omnix_data="omnix_PD_covar_nomissing.csv"
omnix_prs="omnix.best"
###PPMI
PPMI_data="PPMI_covar.txt"
PPMI_prs="PPMI.best"
###APDGC
APDGC_data="APDGC_covar.txt"
APDGC_prs="APDGC.best"
###IPDGC
IPDGC_data="IPDGC_covar.txt"
IPDGC_prs="IPDGC.best"
###NIND
NIND_data="NIND_covar.txt"
NIND_prs="NIND.best"
###NGRC
NGRC_data="NGRC_covar.txt"
NGRC_prs="NGRC.best"
###UKBB
UKBB_data="UKBB_pheno_covariate_downsampled.csv"
UKBB_prs="UKBB.best"
cohort <- c("mcgill","PPMI","APDGC","IPDGC","NIND","NGRC","UKBB")
data <- c(omnix_data,PPMI_data,APDGC_data,IPDGC_data,NIND_data,NGRC_data,UKBB_data)
prs <- c(omnix_prs,PPMI_prs,APDGC_prs,IPDGC_prs,NIND_prs,NGRC_prs,UKBB_prs)
## create a table and iterate through it
df_dir <- data.frame(cohort = cohort, data = data, prs = prs)
## pathway list (only 1 in this case)
name=c("The_NLRP3_inflammasome")
output_folder = "/path/to/output"
#cohort sample size
#3651, 581,923,10709, 1686, 3940, 6621
main_with_multiple_pathway_logOR(df_dir,name,output_folder)
