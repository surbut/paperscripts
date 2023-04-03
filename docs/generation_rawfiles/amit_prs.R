#rm(list=ls())
library('data.table')
library("CVrisk")
library("dplyr")
#######
prs=fread("~/Dropbox/phenotypes/CAD_prs.txt")
dim(na.omit(prs)) 
# 487409      3
cad_pheno=fread("~/Dropbox (Personal)/phenotypes/UKBB_HardCAD_202109.tsv")[,c(1:7,10)]
#cad_pheno=fread("~/Dropbox (Personal)/phenotypes/Coronary_Artery_Disease_HARD.tab.tsv")[,c(1:8,11)]
cad_pheno[cad_pheno$prevalent_disease==1,"incident_disease"]=0
cad_pheno$phenos.enrollment=cad_pheno$enroll_age
cad_pheno$phenos.has_CAD=cad_pheno$has_disease
dim(cad_pheno)
# 502485      7 Starting with 

### 
d2=fread("~/Dropbox/output/Diabetes_Type_2.tab.tsv.gz")
d2[d2$prevalent_disease==1,"incident_disease"]=0
d3=data.frame(sample_id=d2$sample_id,prev_disease_dm=d2$prevalent_disease,inc_disease_dm=d2$incident_disease)
cad_pheno=merge(cad_pheno,d3,by="sample_id")
dim(cad_pheno)
# 502485      9

m=merge(prs,cad_pheno,by.x = "SAMPLE",by.y = "sample_id")

dim(m)
# 487276      9 that have wc'd prs

#subeject with pheno and no PRS
nrow(cad_pheno)-nrow(m)

### lipids info
lipids=fread("~/Dropbox/ukbb-lipids-meds.txt")
sum(rowSums(is.na(lipids[,c("eid","statin0","ldladj","hdladj","choladj","anylipidmed0")]))>0)
# exclude 73414
lip=na.omit(lipids[,c("eid","statin0","ldladj","hdladj","choladj","anylipidmed0")])
dim(lip)
# 429122     6

# baselin info
ukb=fread("~/Dropbox/big_ukb_file.txt",header=T,sep="\t")
cov=data.frame("iid"=ukb$id,"sex"=ukb$Sex,"sbp"=ukb$SBP_adjMeds,"smoke"=ukb$SmokingStatus,"bp_med"=ukb$BP_Meds,"Race"=ukb$in_white_British_ancestry_subset,"bmi"=ukb$BMI)
cov$Race=ifelse(cov$Race==1,"white","other")
##
sum(rowSums(is.na(cov))>0)
## exclude 48973 without baseling info

dim(cov)
cov=na.omit(cov)
dim(na.omit(cov))
# 453656      8

length(setdiff(prs$SAMPLE,cad_pheno$sample_id))
# 133 individuals had PRS and not phenotypes,
m2=merge(m,cov,by.x = "SAMPLE",by.y = "iid")

dim(na.omit(m2))
# 452641     18

## difference in baseline covariates
nrow(m)-nrow(m2)

df=merge(lip,m2,by.x="eid",by.y="SAMPLE")
dim(df)
#393572     23   21 intersection between 429K with lipids and folks with prs and phenotypes
## missing lipid info
nrow(m2)-nrow(df)



n=nrow(df)

## now remove statin folks
dim(df[which(df$statin0==1),])
df=df[-which(df$statin0==1),]
dim(df)
# 329740     23
## on statin
n-nrow(df)

## compute ascvd scores

df$sex=ifelse(as.factor(df$sex)=="Male","male","female")
df$old_smoke=df$smoke
df$smoke=ifelse(df$smoke==2,1,0)
df$SCORESUM=df$prs=df$SCORE
#df$prs_quant=(df$SCORESUM-mean(df$SCORESUM))/sd(df$SCORESUM)
df$prs_quant=(scale(df$prs))
df$prs.r=(rank(df$prs_quant)/length(df$prs_quant))
df$phenos.CAD_censor_age=df$censor_age

df$age.r=(rank(df$phenos.CAD_censor_age)/length(df$prs_quant))
df$ldladj.r=(rank(df$ldladj)/length(df$ldladj))
rownames(df)=as.character(df$eid)

df2=compute_CVrisk(df,scores = c("ascvd_10y_accaha"),
                   age = "phenos.enrollment", race = "Race", gender = "sex", bmi = "bmi", sbp = "sbp",hdl = "hdladj", totchol = "choladj", bp_med = "bp_med", smoker = "smoke",
                   diabetes = "prev_disease_dm", lipid_med = "anylipidmed0",
                   fh_heartattack = NULL, cac = NULL)

df3=compute_CVrisk2(df,scores = c("as2"),
                   age = "phenos.enrollment", race = "Race", gender = "sex", bmi = "bmi", sbp = "sbp",hdl = "hdladj", totchol = "choladj", bp_med = "bp_med", smoker = "smoke",
                   diabetes = "prev_disease_dm", lipid_med = "anylipidmed0",
                   fh_heartattack = NULL, cac = NULL)
df$ascvd_10y_accaha=df2$ascvd_10y_accaha
df$ascvd_10y_accaha_all=df3$as2
df$ascvd.r=rank(df$ascvd_10y_accaha)/length(df$ascvd_10y_accaha)
df$ascvdcat_all <- cut(df$ascvd_10y_accaha_all, breaks=c(0, 7.5, 20,Inf), labels=c("low", "intermediate","high"))
df$ascvdcat <- cut(df$ascvd_10y_accaha, breaks=c(0, 7.5, 20,Inf), labels=c("low", "intermediate","high"))
df$prscat <- cut(df$prs.r, breaks=c(0, 0.20,0.80,0.975,1), labels=c("low", "intermediate","int-high","high"))
df$prscat <- cut(df$prs.r, breaks=c(0, 0.20,0.80,1), labels=c("low", "intermediate","high"))
df$age_to_censor=df$phenos.CAD_censor_age-df$phenos.enrollment

sum(df$age_to_censor<0)
dfc=df[-which(df$age_to_censor<0),]
# 
# 
# dim(dfc)
# dim(dfc)



saveRDS(dfc,"~/Dropbox/paper_scripts/output/amit_df.rds")

#####


library(khsmisc)
library(khsverse)

design <- tibble::tribble(
  ~left,               ~n_left, ~right,              ~n_right,
  "Study base",        502485,    "Lack quality-controlled genotypes for risk score generation",       15209 ,
  "Polygenic score and Phenotypic information",  487276  ,     "Participants with\nmissing baseline covariate information (smoking, anti-hypertensive use, race, and sex)", 34635,
  "With baseline information", 452641,     "Missing Lipid information",                  59069,
  "With Lipid information",  393572,     "On Statin therapy",                 63832,
  "Polygenic score and Phenotypic information off statin therapy",  329740,     "Remains with coronary artery disease",                1903,
  "Complete-case set", 327837,     "",                  NA_integer_)
  

# Plot

pdf("Figs/exclusion_flow.pdf")
e=exclusion_flowchart(design, width = 2)
dev.off()




table1(~ factor(sex) + ldladj + sbp + factor(bp_med) +factor(smoke) +
         factor(prev_disease_dm) +factor(incident_disease) + factor(prscat) + 
         phenos.enrollment +prs.r|ascvdcat,data=df)




# 
# prs=fread("~/Dropbox/phenotypes/CAD_prs.txt")
# lipids=read.table("~/Dropbox/ukbb-lipids-meds.txt",sep = "\t",header = T)
# phenos=readRDS("~/Dropbox/phenotypes/merged_phenos.rds")
# phenos$incident_disease.x[phenos$prevalent_disease.x==1]=0
# phenos$incident_disease.y[phenos$prevalent_disease.y==1]=0
# phenos$incident_disease.x.1[phenos$prevalent_disease.x.1==1]=0
# 
# phenos$true_censor=as.numeric(as.character(apply(phenos,1,function(x){min(x["censor_age.x"],x["death_censor_age.x"])})))
# 
# phenos$true_censor.y=as.numeric(as.character(apply(phenos,1,function(x){min(x["censor_age.y"],x["death_censor_age.y"])})))
# 
# phenos$true_censor_x.1=as.numeric(as.character(apply(phenos,1,function(x){min(x["censor_age.x.1"],x["death_censor_age.x.1"])})))
# 
# 
# phenos=data.frame("iid"=phenos$sample_id,"enrollment"=phenos$enroll_age.x,"has_CAD"=phenos$has_disease.x,"CAD_censor_age"=phenos$true_censor,
#                           "DM"=phenos$disease.y,"has_DM"=phenos$prevalent_disease.y,"Inc_DM_censor_age"=phenos$true_censor.y,"HTN"=phenos$disease.x.1,
#                   "has_HTN"=phenos$has_disease.x.1,"HTN_censor_age"=phenos$true_censor_x.1,"Death_Age"=phenos$death_censor_age.x,"has_Died"=phenos$has_died.x)
# 
# 
# ukb=read.table("~/Dropbox/big_ukb_file.txt",header=T,sep="\t")
# cov=data.frame("iid"=ukb$id,"sex"=ukb$Sex,"sbp"=ukb$SBP_adjMeds,"DM"=ukb$Prev_Diabetes_All,"smoke"=ukb$SmokingStatus,"bp_med"=ukb$BP_Meds,"Race"=ukb$in_white_British_ancestry_subset,"bmi"=ukb$BMI)
# cov$Race=ifelse(cov$Race==1,"white","other")
# #cov=readRDS("~/Dropbox/cov.rds")
# dim(na.omit(cov))
# #447143      8
# ### ensure cases of HTN/DM before CAD
# # 
# # g=grep(pattern = "_age",x =
# #          names(phenos))
# # 
# # temp=data.frame(apply(phenos[,g[-1]],2,function(x){
# #   x<phenos[,"CAD_censor_age"]|(x==phenos[,"CAD_censor_age"]&phenos[,"has_CAD"]==1)}))
# # colnames(temp)=c("DM_prevCAD","HTN_prevCAD")
# # sum((phenos$has_HTN==0&temp$HTN_prevCAD==TRUE))
# # sum(phenos$has_DM==0&temp$DM_prevCAD==TRUE)
# # w=which(phenos$has_DM==0&temp$DM_prevCAD==TRUE)
# # temp[w,"DM_prevCAD"]=FALSE
# # w=which(phenos$has_HTN==0&temp$HTN_prevCAD==TRUE)
# # temp[w,"HTN_prevCAD"]=FALSE
# # 
# # sum((phenos$has_HTN==0&temp$HTN_prevCAD==TRUE))## check to make sure no
# # #if trait wasn't there not counted as TREU
# # sum(phenos$has_DM==0&temp$DM_prevCAD==TRUE)
# # 
# # temp=data.frame(apply(temp,2,function(x){ifelse(x==TRUE,1,0)}))
# 
# #prs=prs[,c("IND_ID","RAW_EUR_CAD")];colnames(prs)=c("iid","prs")
# m=merge(prs,phenos,by.x = "SAMPLE",by.y = "iid")
# #rm(p2)
# dim(na.omit(lipids[,c("eid","statin0","ldladj","hdladj","choladj","anylipidmed0")]))
# df=na.omit(merge(lipids[,c("eid","statin0","ldladj","hdladj","choladj","anylipidmed0")],m,by.x="eid",by.y="SAMPLE"))
# dim(df)
# # rm(m)
# # rm(prs)
# df=merge(df,cov,by.x = "eid",by.y="iid")
# df=na.omit(df)
# # 387932 
# df=df[-which(df$statin0==1),]
# dim(df)
# # rm(cov)
# #df$SCORESUM=-df$SCORESUM
# 
# df$sex=ifelse(as.factor(df$sex)=="Male","male","female")
# 
# df$smoke=ifelse(df$smoke==2,1,0)
# df$prs=df$SCORE
# #df$SCORESUM=df$RAW_EUR_CAD
# #df$prs_quant=(df$SCORESUM-mean(df$SCORESUM))/sd(df$SCORESUM)
# df$prs_quant=(scale(df$prs))
# df$phenos.enrollment=df$enrollment
# df$phenos.CAD_censor_age=df$CAD_censor_age
# df$prs.r=(rank(df$prs_quant)/length(df$prs_quant))
# df$age.r=(rank(df$CAD_censor_age)/length(df$prs_quant))
# df$ldladj.r=(rank(df$ldladj)/length(df$ldladj))
# rownames(df)=as.character(df$eid)
# #df$race=ifelse(ukb[rownames(df),"in_white_British_ancestry_subset"],"white","other")
# #df$htn.r=(rank(df$HTN_prevCAD)/length(df$ldladj))
# #df$ldladj.r=(rank(df$ldladj)/length(df$ldladj))
# suppressPackageStartupMessages(library("CVrisk"))
# df2=compute_CVrisk(df,scores = c("ascvd_10y_accaha", "ascvd_10y_frs", "ascvd_10y_frs_simple"),
#                    age = "phenos.enrollment", race = "Race", gender = "sex", bmi = "bmi", sbp = "sbp",hdl = "hdladj", totchol = "choladj", bp_med = "bp_med", smoker = "smoke",
#                    diabetes = "has_DM", lipid_med = "anylipidmed0",
#                    fh_heartattack = NULL, cac = NULL)
# 
# df$ascvd.r=rank(df$ascvd_10y_accaha)/length(df$ascvd_10y_accaha)
# 
# df$ascvd_10y_accaha=df2$ascvd_10y_accaha
# df$ascvdcat <- cut(df2$ascvd_10y_accaha, breaks=c(0, 7.5, 20,Inf), labels=c("low", "intermediate","high"))
# df$prscat <- cut(df$prs.r, breaks=c(0, 0.20,0.80,0.975,1), labels=c("low", "intermediate","int-high","high"))
# df$prscat <- cut(df$prs.r, breaks=c(0, 0.20,0.80,1), labels=c("low", "intermediate","high"))
# df$age_to_censor=df$phenos.CAD_censor_age-df$phenos.enrollment
# dfc=df[-which(df$age_to_censor<0),]
# 
# 
# dim(dfc)
# 
# ################
# 
# 


dfukb_baseline=readRDS("output/dfukbbaseline.rds")
dat2=merge(x = dad,y = dfukb_baseline,by.x = "eid",by.y = "identifier",all.x = TRUE)
dat2$bp_med2=ifelse(dat2$f.6153.0.0==2|dat2$f.6153.0.1==2|dat2$f.6153.0.2==2|dat2$f.6177.0.0==2|dat2$f.6177.0.1==2|dat2$f.6177.0.2==2,1,0)
dat2$bp_med2[which(is.na(dat2$bp_med2))]=0

dad=dat
dad$bp_med2=dat2$bp_med2
saveRDS(da,"~/paperscripts/output/amit_df_bp.rds")

compute_CVrisk2=function (df, scores = c("ascvd_10y_accaha", "as2", 
                         "ascvd_10y_frs_simple", "chd_10y_mesa", "chd_10y_mesa_cac"), 
          age, gender, race, sbp = NULL, bmi = NULL, hdl = NULL, totchol = NULL, 
          bp_med = NULL, smoker = NULL, diabetes = NULL, lipid_med = NULL, 
          fh_heartattack = NULL, cac = NULL) 
{
  all_args <- as.list(environment())
  valid_pred <- c("age", "gender", "race", "sbp", "bmi", "hdl", 
                  "totchol", "bp_med", "smoker", "diabetes", "lipid_med", 
                  "fh_heartattack", "cac")
  pred_args <- list()
  for (var in valid_pred) {
    if (!is.null(all_args[[var]])) 
      pred_args[[var]] <- df[[all_args[[var]]]]
  }
  results <- sapply(scores, function(x) do.call(x, pred_args))
  row.names(results) <- NULL
  return(cbind(df, results))
}

as2=function (race = "white", gender = c("male", "female"), age, 
          totchol, hdl, sbp, bp_med, smoker, diabetes, ...) 
{
  if (!all(gender %in% c("male", "female")) | missing(gender)) {
    stop("gender must be either 'male' or 'female'")
  }
  ascvd_pooled_coef <- NULL
  utils::data(ascvd_pooled_coef, envir = environment())
  race <- ifelse(race %in% c("white", "aa"), race, "white")
  race_sex <- data.frame(race, gender)
  race_sex$id <- as.numeric(row.names(race_sex))
  pooled_coef <- merge(race_sex, ascvd_pooled_coef)
  pooled_coef <- pooled_coef[order(pooled_coef$id), ]
  sbp_treated <- ifelse(bp_med == 1, sbp, 1)
  sbp_untreated <- ifelse(bp_med == 0, sbp, 1)
  indv_sum <- log(age) * pooled_coef$ln_age + log(age)^2 * 
    pooled_coef$ln_age_squared + log(totchol) * pooled_coef$ln_totchol + 
    log(age) * log(totchol) * pooled_coef$ln_age_totchol + 
    log(hdl) * pooled_coef$ln_hdl + log(age) * log(hdl) * 
    pooled_coef$ln_age_hdl + log(sbp_treated) * pooled_coef$ln_treated_sbp + 
    log(sbp_treated) * log(age) * pooled_coef$ln_age_treated_sbp + 
    log(sbp_untreated) * pooled_coef$ln_untreated_sbp + log(sbp_untreated) * 
    log(age) * pooled_coef$ln_age_untreated_sbp + smoker * 
    pooled_coef$smoker + smoker * log(age) * pooled_coef$ln_age_smoker + 
    diabetes * pooled_coef$diabetes
  risk_score <- round((1 - (pooled_coef$baseline_survival^exp(indv_sum - 
                                                                pooled_coef$group_mean))) * 100, 2)
  ifelse(risk_score < 1, 1, risk_score)
}
#  