# last step: perform meta analysis for all pathways from multiple cohorts

#remotes::install_github("MathiasHarrer/dmetar")

library(dmetar)
library(data.table)
library(glue)
library(esc)
library(meta)
suppressMessages(library(tidyverse))
options(tidyverse.quiet = TRUE)
library(magrittr)
library(metafor)


# meta analysis for one pathway
main_pathway = function(pathway){
    df = read_csv(glue("{pathway}_correlation_analysis(logOR).csv"))
    df[df$cohort == "mcgill","cohort"] = "McGill"
    m.gen <- metagen(TE = logOR,
                    seTE = SE,
                    studlab = cohort,
                    data = df,
                    method.random.ci = "classic",
                    fixed = TRUE,
                    pval= P,
                    random = TRUE,
                    method.tau = "REML",
                    hakn = TRUE,
                    title = glue("{pathway}_pathway_PRS"))
    #convert to OR
    m.gen$TE.common = exp(m.gen$TE.common)
    m.gen$seTE.common = m.gen$TE.common * m.gen$seTE.common
    m.gen$lower.common = m.gen$TE.common - 1.96 * m.gen$seTE.common
    m.gen$upper.common = m.gen$TE.common + 1.96 * m.gen$seTE.common
    m.gen$TE.random = exp(m.gen$TE.random)
    m.gen$seTE.random = m.gen$TE.random * m.gen$seTE.random
    m.gen$lower.random = m.gen$TE.random - 1.96 * m.gen$seTE.random
    m.gen$upper.random = m.gen$TE.random + 1.96 * m.gen$seTE.random
    #convert cohort-level data to OR
    m.gen$TE = exp(m.gen$TE)
    m.gen$seTE = m.gen$TE * m.gen$seTE
    m.gen$lower = m.gen$TE - 1.96 * m.gen$seTE
    m.gen$upper = m.gen$TE + 1.96 * m.gen$seTE
    print(summary(m.gen))
    #change the plotwidth to 7cm if you only wanna include either random effect or fixed effect
    pdf(file = glue("{pathway}_forestplot.pdf"), width = 10, height = 5)
    forest.meta(m.gen,
            leftcols = c("cohort"),
            leftlabs = c("Cohort"),
            rightlabs = c("OR","95%-CI", "Weight(Fixed)","Weight(Random)"),
            text.common = "Fixed effect model",
            text.random = "Random effect model  ",
            xlim = c(0,2),
            plotwidth = "5cm",
            ref = 1,
            fontsize = 14,
            comb.fixed = TRUE,
            comb.random = TRUE,
            print.I2 = TRUE,
            print.tau2 = FALSE,
            lty.random = NULL,
            fs.hetstat = 10)
    dev.off()
    return(m.gen)
}

# main program
main = function(){
  name=c("The_NLRP3_inflammasome"
  )
  result = NULL
  for (pathway in name){
    m.gen = main_pathway(pathway)
    result %<>% rbind(.,
        data.frame(pathway=pathway,
                    Fixed = as.numeric(m.gen$TE.common), 
                    Fixed_se = as.numeric(m.gen$seTE.common),
                    Fixed_lower = as.numeric(m.gen$lower.common),
                    Fixed_upper = as.numeric(m.gen$upper.common),
                    Fixed_pval = as.numeric(m.gen$pval.common),
                    Random = as.numeric(m.gen$TE.random), 
                    Random_se = as.numeric(m.gen$seTE.random),
                    Random_lower = as.numeric(m.gen$lower.random),
                    Random_upper = as.numeric(m.gen$upper.random),
                    Random_pval = as.numeric(m.gen$pval.random),
                    Heterogeneity_I2 = as.numeric(m.gen$I2),
                    Heterogeneity_I2_lower = as.numeric(m.gen$lower.I2),
                    Heterogeneity_I2_upper = as.numeric(m.gen$upper.I2),
                    Heterogeneity_I2_pval = as.numeric(m.gen$pval.Q)
                   )
        )
  }
  write.csv(result,"meta_analysis_pathway_prs.csv",row.names=FALSE)
  return(result)
}

results = main()



