# --------------------------------------------------------------------------------
#Title: Dietary quality scores, metabolic signature, and cardiometabolic diseases
#Purpose:Estimate contribution of overall microbial community to MS
#Study:MLVS
#Programmer:Huan Yun
#Date:20231012
#Note:1) we did not transform MS; 2) for microbiome, we use TSS and scale them; 3) we did not fix any covariate in the model; 4) we only use filtered species
# --------------------------------------------------------------------------------

#load packages
library(data.table)
library(readxl)
library(dplyr)
library(openxlsx)
library(cvTools)
library(glmnet)
library(CMAverse)

#set working dir
setwd("/udd/nhhyu/DP_T2D/FinalOutput")

#--------------------------------------------------------------------------------------------
#
#           Step1: run elastic net model to get the contribution of overall microbiome to the variation of metabolic signatures
#
#--------------------------------------------------------------------------------------------
#load mlvs data
load("/udd/nhhyu/DP_T2D/data/Microbiome/ToBeUsed.Taxon.RData") #dataset was created by "FigureS13.R" script

#load filtered results
taxa_c <- read_excel("/udd/nhhyu/DP_T2D/results/Feature-wide association results.xlsx", sheet = "sheet_a") %>% as.data.frame()
microvar <- subset(taxa_c, type == "species")$microName
length(microvar) #151

#standarise microbiome data (TSS transformation and then scale)
Taxon_MLVS_Std <- ToBeUsed.Taxon
Taxon_MLVS_Std[microvar] <- apply(Taxon_MLVS_Std[microvar], 2, scale)

#get name of dietary score
diet <- c("amed_av","ahei_av","dash_av","pdi_av","hpdi_av","updi_av","edip_av","edih_av")
ms <- c("amed_c","ahei_c","dash_c","pdi_c","hpdi_c","updi_c","edip_c","edih_c")

#define dataset
df <- as.data.frame(Taxon_MLVS_Std)

#define interaction term
for (i in diet){
  for (j in microvar){
    df$temp=df[,i]*df[,j]
    names(df)[which(names(df)=="temp")]= paste("int_",i,"_",j,sep='') 
  }
}

#define outcome
outvar <- ms

#define microbiome
microvar

#define covariates
covar <- c("age_fec","bmi_bld","totMETs_paq","smoke_bld","probio_2m_fec","antibio_12m_fec","colsc_2m_fec","acid_2m_fec",
           "stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6")

#check data including MS (raw value) and microbiome (TSS+scale)
rbind(summary(Taxon_MLVS_Std$amed_c),summary(Taxon_MLVS_Std$ahei_c),summary(Taxon_MLVS_Std$dash_c),summary(Taxon_MLVS_Std$pdi_c),
      summary(Taxon_MLVS_Std$hpdi_c),summary(Taxon_MLVS_Std$updi_c),summary(Taxon_MLVS_Std$edip_c),summary(Taxon_MLVS_Std$edih_c))

#split data
nrFolds <- 10
folds <- sample(rep_len(1:nrFolds, nrow(df)),replace=F) 

#run elastic net model
Predict = data.frame(amed=NA,ahei=NA,dash=NA,pdi=NA,hpdi=NA,updi=NA,edip=NA,edih=NA)
for (y_bi in outvar){
  for(k in 1:10){ 
    fold <- which(folds == k) 
    data.train <- df[-fold,]
    data.test <- df[fold,]
    x.train <- as.matrix(data.train[,microvar])
    y.train <- as.matrix(data.train[,y_bi])
    x.test <- as.matrix(data.test[,microvar])
    cv <- cv.glmnet(x.train, y.train, alpha = 0.5, nfold=10) 
    Predict[fold,y_bi] = data.frame(predict(cv, x.test ,type="response", s="lambda.min"))
  }
}

#define function to estimate variability
rsq <- function(x, y) summary(lm(y~x))$r.squared

#calculate R and R square
performance <- cbind(c(rsq(df$amed_c,Predict$amed_c),rsq(df$ahei_c,Predict$ahei_c),rsq(df$dash_c,Predict$dash_c),rsq(df$pdi_c,Predict$pdi_c),
                        rsq(df$hpdi_c,Predict$hpdi_c),rsq(df$updi_c,Predict$updi_c),rsq(df$edip_c,Predict$edip_c),rsq(df$edih_c,Predict$edih_c))) %>% as.data.frame()

rownames(performance) <- outvar
colnames(performance) <- "R2"

#save results
write.csv(performance,file = "Figure7.1_MLVS.csv")

#--------------------------------------------------------------------------------------------
#
#           Step2: test the mediated role of overall microbiome in diet-metabolic signature relationship
#
#--------------------------------------------------------------------------------------------
#combine predicted score with mediation analysis
exp <- c("amed_av","ahei_av","dash_av","pdi_av","hpdi_av","updi_av","edip_av","edih_av")
med <- c("amed_p","ahei_p","dash_p","pdi_p","hpdi_p","updi_p","edip_p","edih_p")

Predict <- Predict[,-c(1:8)]
names(Predict) <- med

df2 <- cbind(df,Predict[med])

#inverse-normal transformation
inormal <- function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

df2[,c(exp,med)] <- apply(df2[,c(exp,med)],2,inormal)

out <- c("amed_c","ahei_c","dash_c","pdi_c","hpdi_c","updi_c","edip_c","edih_c")
df2[,out] <- apply(df2[,out],2,scale)

#get controlled median value of mediators
metab_dataset <- df2 %>% dplyr::select(med)

controlled_med<- metab_dataset%>%
  dplyr::select(med)%>%data.table::melt()%>%
  group_by(variable)%>%
  summarize(median=median(value,na.rm = T))

controlled_med <- as.vector(controlled_med$median)

#check covariates
ba_cov <- apply(df2[covar],2,summary)

#check exposures-missing
ba_exp <- apply(df2[exp],2,summary)

#check mediators-no missing values
ba_med <- apply(df2[med],2,summary)

#check outcomes-no misssing values
ba_out <- apply(df2[out],2,summary)

#run mediation model
cmest1 <- cmest(data=subset(df2,amed_av>-10),model="rb",outcome="amed_c",event=NULL,exposure="amed_av",mediator="amed_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=100,yvar=1)
cmest2 <- cmest(data=subset(df2,ahei_av>-10),model="rb",outcome="ahei_c",event=NULL,exposure="ahei_av",mediator="ahei_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=100,yvar=1,multimp=T) 
cmest3 <- cmest(data=subset(df2,dash_av>-10),model="rb",outcome="dash_c",event=NULL,exposure="dash_av",mediator="dash_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=100,yvar=1,multimp=T) 
cmest4 <- cmest(data=subset(df2,pdi_av>-10),model="rb",outcome="pdi_c",event=NULL,exposure="pdi_av",mediator="pdi_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=100,yvar=1,multimp=T) 
cmest5 <- cmest(data=subset(df2,hpdi_av>-10),model="rb",outcome="hpdi_c",event=NULL,exposure="hpdi_av",mediator="hpdi_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=100,yvar=1,multimp=T) 
cmest6 <- cmest(data=subset(df2,updi_av>-10),model="rb",outcome="updi_c",event=NULL,exposure="updi_av",mediator="updi_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=100,yvar=1,multimp=T) 
cmest7 <- cmest(data=subset(df2,edip_av>-10),model="rb",outcome="edip_c",event=NULL,exposure="edip_av",mediator="edip_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=100,yvar=1,multimp=T) 
cmest8 <- cmest(data=subset(df2,edih_av>-10),model="rb",outcome="edih_c",event=NULL,exposure="edih_av",mediator="edih_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=100,yvar=1,multimp=T)

#get mediation results
res1 <- as.data.frame(as.matrix(summary(cmest1)$summarydf)); res1$mediator <- "amed"
res2 <- as.data.frame(as.matrix(summary(cmest2)$summarydf)); res2$mediator <- "ahei"
res3 <- as.data.frame(as.matrix(summary(cmest3)$summarydf)); res3$mediator <- "dash"
res4 <- as.data.frame(as.matrix(summary(cmest4)$summarydf)); res4$mediator <- "pdi"
res5 <- as.data.frame(as.matrix(summary(cmest5)$summarydf)); res5$mediator <- "hpdi"
res6 <- as.data.frame(as.matrix(summary(cmest6)$summarydf)); res6$mediator <- "updi"
res7 <- as.data.frame(as.matrix(summary(cmest7)$summarydf)); res7$mediator <- "edip"
res8 <- as.data.frame(as.matrix(summary(cmest8)$summarydf)); res8$mediator <- "edih"

res_med <- rbind(res1,res2,res3,res4,res5,res6,res7,res8)

#save mediation results
write.csv(res_med, file = "Figure7.2_MLVS.csv")

#--------------------------------------------------------------------------------------------
#
#           Step3: test the hypothesis for doing mediation analysis
#
#--------------------------------------------------------------------------------------------
#define exposure, mediator, and outcome
cov <- covar
out <- c("amed_c","ahei_c","dash_c","pdi_c","hpdi_c","updi_c","edip_c","edih_c")

df2[,out] <- apply(df2[,out],2,scale)

#association between exposure and outcome
b <- matrix(NA,8,8)
se <- matrix(NA,8,8)
p <- matrix(NA,8,8)

rownames(b) <- rownames(se) <- rownames(p) <- exp
colnames(b) <- colnames(se) <- colnames(p) <- out

for (i in exp){
  for (j in out){
    df2$Score <- df2[,i]
    df2$y <- df2[,j]
    fit <- lm(y~Score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6,data=df2)
    b[i,j] = coef(summary(fit))["Score",1]
    se[i,j] = coef(summary(fit))["Score",2]
    p[i,j] = coef(summary(fit))["Score",4]
  }
}

med1 <- cbind(b[,1],se[,1],p[,1],b[,2],se[,2],p[,2],b[,3],se[,3],p[,3],b[,4],se[,4],p[,4],b[,5],se[,5],p[,5],b[,6],se[,6],p[,6],b[,7],se[,7],p[,7],b[,8],se[,8],p[,8]) 

#association between exposure and mediator
b <- matrix(NA,8,8)
se <- matrix(NA,8,8)
p <- matrix(NA,8,8)

rownames(b) <- rownames(se) <- rownames(p) <- exp
colnames(b) <- colnames(se) <- colnames(p) <- med

for (i in exp){
  for (j in med){
    df2$Score <- df2[,i]
    df2$y <- df2[,j]
    fit <- lm(y~Score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6,data=df2)
    b[i,j] = coef(summary(fit))["Score",1]
    se[i,j] = coef(summary(fit))["Score",2]
    p[i,j] = coef(summary(fit))["Score",4]
  }
}

med2 <- cbind(b[,1],se[,1],p[,1],b[,2],se[,2],p[,2],b[,3],se[,3],p[,3],b[,4],se[,4],p[,4],b[,5],se[,5],p[,5],b[,6],se[,6],p[,6],b[,7],se[,7],p[,7],b[,8],se[,8],p[,8]) 

#association between mediator and outcome
b <- matrix(NA,8,8)
se <- matrix(NA,8,8)
p <- matrix(NA,8,8)

rownames(b) <- rownames(se) <- rownames(p) <- med
colnames(b) <- colnames(se) <- colnames(p) <- out

for (i in med){
  for (j in out){
    df2$Score <- df2[,i]
    df2$y <- df2[,j]
    fit <- lm(y~Score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6,data=df2)
    b[i,j] = coef(summary(fit))["Score",1]
    se[i,j] = coef(summary(fit))["Score",2]
    p[i,j] = coef(summary(fit))["Score",4]
  }
}

med3 <- cbind(b[,1],se[,1],p[,1],b[,2],se[,2],p[,2],b[,3],se[,3],p[,3],b[,4],se[,4],p[,4],b[,5],se[,5],p[,5],b[,6],se[,6],p[,6],b[,7],se[,7],p[,7],b[,8],se[,8],p[,8]) 

#combine all results
med_tot <- rbind(c(med1[1,1:3],med2[1,1:3],med3[1,1:3]),
                        c(med1[2,4:6],med2[2,4:6],med3[2,4:6]),
                              c(med1[3,7:9],med2[3,7:9],med3[3,7:9]),
                                    c(med1[4,10:12],med2[4,10:12],med3[4,10:12]),
                                          c(med1[5,13:15],med2[5,13:15],med3[5,13:15]),
                                                c(med1[6,16:18],med2[6,16:18],med3[6,16:18]),
                                                      c(med1[7,19:21],med2[7,19:21],med3[7,19:21]),
                                                            c(med1[8,22:24],med2[8,22:24],med3[8,22:24])) %>% as.data.frame

colnames(med_tot) <- c("Beta1","SE1","P1","Beta2","SE2","P2","Beta3","SE3","P3")
rownames(med_tot) <- c("AMED","AHEI","DASH","PDI","hPDI","uPDI","EDIP","EDIH")

#save results 
write.csv(med_tot, file = "TableS15_MLVS.csv")