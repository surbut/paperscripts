library('reshape') 
library('dplyr') 
library("survival") 
library("ggplot2") 
library("data.table") 

fwk=read.table("~/Desktop//Framingham/85598/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v32.p13.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v32.pht006027.v3.p13.c1.vr_wkthru_ex09_1_1001s.HMB-IRB-MDS.txt",header = T,skip = 1,sep="\t") 
sum(na.omit(fwk$DMRX1)) ## only 37 folks with DM
fwk=fwk[,c("dbGaP_Subject_ID","SEX","AGE1","CURRSMK1","HDL1","TC1","LIPRX1","CALC_LDL1","BMI1","SBP1","DMRX1","HRX1")]

length(which(rowSums(is.na(fwk))>0))
length(which(rowSums(is.na(fwk))==0))

d=fwk[which(rowSums(is.na(fwk))==0),]

## details here https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/dataset.cgi?study_id=phs000007.v33.p14&pht=3316

outcomes=read.table("~/Desktop//Framingham/85598/PhenoGenotypeFiles/RootStudyConsentSet_phs000007.Framingham.v32.p13.c1.HMB-IRB-MDS/PhenotypeFiles/phs000007.v32.pht003316.v9.p13.c1.vr_survcvd_2018_a_1267s.HMB-IRB-MDS.txt",sep="\t",header = T)
df=merge(fwk, outcomes[,c("dbGaP_Subject_ID","shareid","chddate", "chd")])
dim(df)
df$chd.age=df$chddate/365.25+df$AGE1

d2=merge(d, outcomes[,c("dbGaP_Subject_ID","shareid","chddate", "chd")])
dim(d2)
dim(na.omit(d2))

## 3660   15
sum(na.omit(d2$LIPRX1)==1)

d3=d2[-which(d2$LIPRX1==1),]
dim(d2)[1]-dim(d3)[1]
## 21
d4=d3[-which(d3$chddate<=0&d3$chd==1),]
dim(d4)[1]-dim(d3)[1]
##  51
dim(d4)
#3588
saveRDS(df,'~/Dropbox/fh_full.rds')
prs=fread("~/Dropbox (Personal)/Fram_allchr_CAD_c1.profile")

df2=merge(df, prs, by.x = "shareid", by.y = "V1")
d5=merge(d4, prs, by.x = "shareid", by.y = "V1")
dim(d4)[1]-dim(d5)[1]

nrow(df)-nrow((df2))
saveRDS(df2,'~/Dropbox/fh_prs.rds')
sum(na.omit(fwk$LIPRX1)==1)
sum(na.omit(df))

#####



library(khsmisc)
library(khsverse)

design <- tibble::tribble(
  ~left,               ~n_left, ~right,              ~n_right,
  "Study base with Outcome Data",        3821 ,   "Lack Total Cholesterol, HDL, Smoking, Systolic or Blood Pressure Information at baseline ",       154 ,
  "Contain All Covariates",  3660,     "On Statin at Baseline", 21,
  "Contain All Covariates and Statin-free",  3639,     "With CAD at baseline", 51,
  "With baseline covariate and outcome data", 3588,   "Missing PRS Info",                 959,
  "PRS, Pheno, Covariate info", 2629,  "",                  NA_integer_)


##

e=exclusion_flowchart(design, width = 2)
e
###


