### Numbers for paper

#Follow up ranges for FOS

cad=fread("~/Dropbox/UKBB_HardCAD_202109.tsv")
summary(cad$censor_age[cad$has_disease==0]-cad$enroll_age[cad$has_disease==0])

fh=readRDS("~/Dropbox/fh_full.rds")

summary(fh$chddate[fh$chd==0])/365


## Phenotype definiteions in /medpop/esp2/pradeep/UKBiobank/r201910/bigquery/disease/cvdi

## Diabetes_Type_1.tab or Diabetes_Type_2.tab 

## Hazard rates UKB

df=readRDS("output/amit_df.rds")
dat=data.frame(df %>%group_by(round(phenos.enrollment,0)) %>%summarise(length(phenos.enrollment)))[c(3:32),]

hazards=data.frame(readRDS("output/hazards_rate_UKB_amit.rds"))
hrmat=hazards[,grep("HR",colnames(hazards))];
hrmat$age=dat$round.phenos.enrollment..0.
summary(hrmat)

upperbound=data.frame(readRDS("output/hazards_rate_upperbound_amit.rds"))
lowerbound=data.frame(readRDS("output/hazards_rate_lowerbound_amit.rds"))

lowerbound$age=upperbound$age=dat$round.AGE1..0.

## Hazard rates FHS
df=readRDS("~/Dropbox/fh_full.rds")
dat=data.frame(df %>%group_by(round(AGE1,0)) %>%summarise(length(AGE1)))
dat=dat[c(11:50),]

hazards=readRDS("output/hazards_fhs.rds")
hrmat=hazards[,grep("HR",colnames(hazards))];
hrmat$age=dat$round.AGE1..0.
summary(hrmat)

(hrmat[,c("age","prs.HR")])

low=data.frame(readRDS("output/hazards_rate_lowerbound_fh.rds"))
upper=data.frame(readRDS("output/hazards_rate_upperbound_fh.rds"))
low$age=upper$age=dat$round.AGE1..0.



### For example, the hazard ratio of lifetime CAD for systolic blood pressure decreased from 4.61 (95% CI XXX-XXX) at age 18y to 
#1.68 (95% CI XXX-XXX) at age 70y. The trends were overall similar for low-density lipoprotein (LDL) cholesterol , smoking, and diabetes (only computed in UKB given that scarcity of diabetes in FOS) (Figure 1A and 1B).  

## At younger ages (40-43 years), genomic risk predicted over 4-fold (95% CI XXX-XXX) more events compared to nongenomic risk 

df=readRDS("~/Dropbox/fh_withprs_pce.rds")
df=df[-which(is.na(df$ascvd_10y_accaha)),]
dff=df

df=na.omit(data.frame(readRDS("output/amit_df.rds")))
df$chd=df$phenos.has_CAD


df$AGE1=df$phenos.enrollment

m=rbind(df[,c("AGE1","prs.r","ascvd_10y_accaha","chd")],dff[,c("AGE1","prs.r","ascvd_10y_accaha","chd")])

df=m
df=df[df$AGE1>=40&df$AGE1<=70,]
i=df %>% mutate(ints=cut(AGE1,breaks=c(39.99,45,50,55,60,65,70,75)))
fulla=i%>%group_by(ints)%>%summarise("PRS"=sum(chd==1&ascvd_10y_accaha<7.5&prs.r>0.80)/sum(chd==1),"PCE"=sum(chd==1&ascvd_10y_accaha>7.5&prs.r<0.80)/sum(chd==1),"Both"=sum(chd==1&ascvd_10y_accaha>7.5&prs.r>0.80)/sum(chd==1),"n"=length(ints))

## variance of a proportion is sqrt(pq/n)

1.96*sqrt((fulla$PRS*(1-fulla$PRS))/fulla$n)+fulla$PRS
1.96*sqrt((fulla$PCE*(1-fulla$PCE))/fulla$n)+fulla$PCE
1.96*sqrt((fulla$Both*(1-fulla$Both))/fulla$n)+fulla$PCE


-1.96*sqrt((fulla$PRS*(1-fulla$PRS))/fulla$n)+fulla$PRS
-1.96*sqrt((fulla$PCE*(1-fulla$PCE))/fulla$n)+fulla$PCE
-1.96*sqrt((fulla$Both*(1-fulla$Both))/fulla$n)+fulla$PCE

### hazard tables

df=readRDS("output/amit_df.rds")
dat=data.frame(df %>%group_by(round(phenos.enrollment,0)) %>%summarise(length(phenos.enrollment)))[c(3:32),]


hazards=data.frame(readRDS("output/hazards_rate_UKB_amit.rds"))
upperbound=data.frame(readRDS("output/hazards_rate_upperbound_amit.rds"))
lowerbound=data.frame(readRDS("output/hazards_rate_lowerbound_amit.rds"))

hrmat=hazards[,grep("HR",colnames(hazards))];
hrmat$age=dat$round.phenos.enrollment..0.

ormat=hazards[,grep("OR",colnames(hazards))];ormat$age=dat$round.phenos.enrollment..0.
r2mat=hazards[,grep("r2",colnames(hazards))];r2mat$age=dat$round.phenos.enrollment..0.

colnames(hrmat)=colnames(ormat)=colnames(r2mat)=c("Polygenic Score","Pooled cohort equations","LDL Cholesterol","Smoking","Age","Systolic blood pressure","Diabetes mellitus","Body Mass Index","Male Sex","High-density Lipoprotein","Total cholesterol","age")
r2mat=round(r2mat,3)
rownames(r2mat)=r2mat$age

upperbound$age=dat$round.phenos.enrollment..0.
u=melt(upperbound,id.vars="age")

lowerbound$age=dat$round.phenos.enrollment..0.
l=melt(lowerbound,id.vars="age")


round(u$value,2)
round(l$value,2)
round(haz$value,2)
## convert to longform

haz = melt(hrmat, id.vars = "age")

u$value=round(u$value,2)
l$value=round(l$value,2)
haz$value=round(haz$value,2)

haz$low.se = l$value
haz$upper.se = u$value

h=cbind(haz[1:nrow(dat),],haz[31:60,],haz[61:90,],
        haz[91:120,],haz[151:180,],haz[181:210,],
        haz[211:240,],haz[241:270,],haz[271:300,],
        haz[301:330,])

write.table(h,"output/hrindex_UKB.csv",quote = F,sep=",",row.names = F)
write.table(r2mat,"output/PVE_UKB.csv",quote = F,sep=",",row.names = F)
rm(list=ls())
## for FOS

df1=readRDS("~/Dropbox/fh_full.rds")
df2=readRDS("~/Dropbox/fh_prs.rds")
dat=data.frame(df1 %>%group_by(AGE1)%>%summarise(l=length(AGE1)))
dat=dat[c(11:50),]

hazards=data.frame(readRDS("output/hazards_fhs.rds"))
upperbound=data.frame(readRDS("output/hazards_rate_upperbound_fh.rds"))
lowerbound=data.frame(readRDS("output/hazards_rate_lowerbound_fh.rds"))


hrmat=hazards[,grep("HR",colnames(hazards))];
hrmat$age=dat$AGE1

ormat=hazards[,grep("OR",colnames(hazards))];ormat$age=dat$AGE1
r2mat=hazards[,grep("r2",colnames(hazards))];r2mat$age=dat$AGE1

colnames(hrmat)=colnames(ormat)=colnames(r2mat)=c("Polygenic Score","LDL Cholesterol","Smoking","Age","Systolic blood pressure","Total cholesterol","Body Mass Index","Male Sex","High-density Lipoprotein","age")

r2mat=round(r2mat,3)
rownames(r2mat)=r2mat$age

colnames(hrmat)=c("Polygenic Score","LDL Cholesterol","Smoking","Age","Systolic blood pressure","Total cholesterol","Body Mass Index","Male Sex","High-density Lipoprotein","age")

upperbound$age=dat$AGE1
u=melt(upperbound,id.vars="age")

lowerbound$age=dat$AGE1
l=melt(lowerbound,id.vars="age")



## convert to longform
haz = melt(hrmat, id.vars = "age")

u$value=round(u$value,2)
l$value=round(l$value,2)
haz$value=round(haz$value,2)

haz$low.se = l$value
haz$upper.se = u$value

h=cbind(haz[1:nrow(dat),],haz[41:80,],haz[81:120,],
        haz[121:160,],haz[161:200,],haz[201:240,],
        haz[241:280,],haz[281:320,],haz[321:360,])

write.table(h,"output/hrindex_FHs.csv",quote = F,sep=",",row.names = F)
write.table(r2mat,"output/PVE_fhs.csv",quote = F,sep=",",row.names = F)