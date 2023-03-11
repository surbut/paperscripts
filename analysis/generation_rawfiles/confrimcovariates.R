
library("data.table")
library("ukbpheno")
library("ggplot2")
library("survival")
library("dplyr")
library("tidyr")

# the directory with datafiles
pheno_dir="~/Dropbox (Personal)/pheno_dir/"


fukbtab <- paste(pheno_dir,"ukb47823.tab",sep="")
# meta data file
fhtml <- paste(pheno_dir,"ukb47823.html",sep="")

# hospital inpatient data
fhesin <- paste(pheno_dir,"ukb_showcase_9.01/hesin.txt",sep="")
fhesin_diag <- paste(pheno_dir,"ukb_showcase_9.01/hesin_diag.txt",sep="")
fhesin_oper <- paste(pheno_dir,"ukb_showcase_9.01/hesin_oper.txt",sep="")


# lst.harmonized.data=readRDS("~/Dropbox/pheno_dir/lst.harmonized.data.rds")
# 
# # # f.53.0.0 contains baseline visit date 
# # df_reference_dt_v0<-lst.harmonized.data$dfukb[,c("identifier","f.53.0.0")]
# 
# # # f.52.0.0 month of birth
# 
# d4<-lst.harmonized.data$dfukb[,c("identifier","f.52.0.0","f.34.0.0")] 
# # # f.34.0.0 contains year of birth, f.52 is month of birth, create birthday on 15th of month
# 
# lst.harmonized.data$dfukb$Birthdate<-as.Date(with(d4,paste(f.34.0.0,f.52.0.0,15,sep="-")),"%Y-%m-%d")
# 
# df_reference_dt_v0<-lst.harmonized.data$dfukb[,c("identifier","Birthdate")]

## create the new variable for reference date

dfhtml <- read_ukb_metadata(fhtml)
# 
# # Rename the identifier column in the metadata
# 
dfhtml[which(dfhtml$field.tab=="f.eid"),]$field.tab<-"identifier"
# 
# # Age at assessment center visit, sex, BMI, HbA1c, glucose,insulin within 1 year of diagnosis,UK Biobank assessment center location, date of visit, LDL cholesterol, HDL, SBP, 
# 
#baseline_fields<-c(21003,31,34,21001,53,54)
# #21001,54,53) #,30780,30760,4080)
# 

#baseline_fields<-c(21003,31,34,21001,53,54)
# #baseline_fields<-c(33) #,30780,30760,4080)
# 
# # Extract these variables from main dataset

baseline_fields<-c(21003,31,6177,6153,20116,20160,4080)

# 
dfukb_baseline <- read_ukb_tabdata(fukbtab,dfhtml,fields_to_keep = baseline_fields)
# 
## smoking status (never, past, current) 20116
## ever smokied (yes no) 20160
## BP 6153 6177 female ales,2	Blood pressure medication 
## 4080 is BP

bigu=fread("~/Dropbox/big_ukb_file.txt")
e=bigu[,c("id","Sex","SmokingStatus","SmokingStatusv2","ever_smoked","SBP","BP_Meds","SBP_adjMeds")]
ub=dfukb_baseline[,c("identifier","f.31.0.0","f.4080.0.0","f.6153.0.0","f.6153.0.1","f.6153.0.2","f.6177.0.0","f.6177.0.1","f.6177.0.2","f.20116.0.0","f.20160.0.0","f.21003.0.0")]

m=merge(e,ub,by.x="id",by.y="identifier")
h=m[which(m$BP_Meds==1),]

## confrimed that has at least one instance
