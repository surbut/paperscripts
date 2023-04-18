library('data.table')
library("CVrisk")
library("dplyr")

prs=fread("~/Dropbox/phenotypes/CAD_prs.txt")
dim(na.omit(prs)) 
# 487409      3
cad_pheno=fread("~/Dropbox (Personal)/phenotypes/UKBB_HardCAD_202109.tsv")[,c(1:7,10)]
df=readRDS("output/amit_df.rds")


## here scale all together
mall=merge(prs,cad_pheno[,c("sample_id","enroll_age")],by.x="SAMPLE",by.y="sample_id")
mall$SCORE=scale(mall$SCORE)
m2=mall[SAMPLE%in%df$eid,]

all=aggregate(SCORE ~ round(enroll_age,0), data = mall, median)
sub=aggregate(SCORE ~ round(enroll_age,0), data = m2, median)
colnames(sub)=colnames(all)

meld=merge(all,sub,by="round(enroll_age, 0)")
colnames(meld)=c("age","all_median","subset_median")
mt=melt(meld[3:33,],id.vars="age")

ggplot(mt,aes(age,y = value,col=variable))+geom_point()+ylim(c(-0.2,0.2))


## now do it with scaling by df
all=aggregate(SCORE ~ round(enroll_age,0), data = mall, median)
sub=aggregate(prs_quant ~ round(phenos.enrollment,0), data = df, median)
colnames(sub)=colnames(all)

meld=merge(all,sub,by="round(enroll_age, 0)")
colnames(meld)=c("age","all_median","subset_median")
mt=melt(meld[3:33,],id.vars="age")

ggplot(mt,aes(age,y = value,col=variable))+geom_point()+ylim(c(-0.2,0.2))

###


hist(df$prs_quant[df$phenos.enrollment<55],freq=FALSE)
lines(density(df$prs_quant[df$phenos.enrollment>55&df$phenos.enrollment<65]))
ines(density(df$prs_quant[df$phenos.enrollment>65))

plot(density(mall$SCORE[mall$enroll_age<55]),col="green",main="Density of PRS Distribution by Age of Enrollment")
lines(density(mall$SCORE[mall$enroll_age>55&mall$enroll_age<65]),col="red")
lines(density(mall$SCORE[mall$enroll_age>65]),col="blue")
legend("topright",legend = c("<55","55-65",">65"),col=c("green","red","blue"),pch=19)

