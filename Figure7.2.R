# --------------------------------------------------------------------------------
#Title:Dietary quality scores, metabolic signature, and cardiometabolic diseases
#Purpose:Estimate contribution of overall microbial community to MS
#Study:MBS
#Programmer:Huan Yun
#Date:20231012
#Note:1)we did not transform MS; 2)for microbiome, we use TSS and scale them; 3)we did not fix any covariate in the model; 4)we only use filtered species
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
#read data
load("/udd/nhhyu/DP_T2D/data/Microbiome/ToBeUsed.Taxon.MBS.RData")
Taxon_MBS <- as.data.frame(cov_avg)

#read diet data from nhs2 (here we use dietary data from NHS2)
nhs2 <- fread("/udd/nhhyu/DP_T2D/data/Phenotype/NHS2/Diet_NHS2_0711.csv")

#merge diet data
Taxon_MBS <- left_join(Taxon_MBS, nhs2[,c("id","amed_11","ahei2010_11","dashav","pdi_11","hpdi_11","updi_11","edip11","edih11")], by = "id")

#get filtered species name
anno_f <- fread("/udd/nhhyu/DP_T2D/data/Microbiome/Filtered_Anno_MBS.csv")
microvar <- anno_f$labname #153

#standarise microbiome data (TSS+scale)
Taxon_MBS_Std <- Taxon_MBS
Taxon_MBS_Std[microvar] <- apply(Taxon_MBS_Std[microvar], 2, scale)

#define function
rsq <- function(x, y) summary(lm(y~x))$r.squared

#define dataset
df <- as.data.frame(Taxon_MBS_Std)

#define outcome
outvar <- c("amed2","ahei2","dash2","pdi2","hpdi2","updi2","edip2","edih2")

#define microbiome
microvar

#define covar
covar <- c("ageyr","race","bmi","act","alco","smoke","calor","stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4",
           "stooltype_fec.5","stooltype_fec.6","stooltype_fec.7","stooltype_fec.8","antibio_12m_fec","colsc_2m_fec","probio_2m_fec","acid_2m_fec")

#check data including MS (raw value) and microbiome (TSS+scale)
rbind(summary(Taxon_MBS_Std$amed2),summary(Taxon_MBS_Std$ahei2),summary(Taxon_MBS_Std$dash2),summary(Taxon_MBS_Std$pdi2),
      summary(Taxon_MBS_Std$hpdi2),summary(Taxon_MBS_Std$updi2),summary(Taxon_MBS_Std$edip2),summary(Taxon_MBS_Std$edih2))

df <- subset(df,amed2>0)

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

#define function
rsq <- function(x, y) summary(lm(y~x))$r.squared

#calculate R and R square
performance <- cbind(c(rsq(df$amed2,Predict$amed2),rsq(df$ahei2,Predict$ahei2),rsq(df$dash2,Predict$dash2),rsq(df$pdi2,Predict$pdi2),
                      rsq(df$hpdi2,Predict$hpdi2),rsq(df$updi2,Predict$updi2),rsq(df$edip2,Predict$edip2),rsq(df$edih2,Predict$edih2))) %>% as.data.frame()

rownames(performance) <- outvar
colnames(performance) <- "R2"

#save results
write.csv(performance,file = "Figure7.1_MBS.csv")

#--------------------------------------------------------------------------------------------
#
#           Step2: test the mediated role of overall microbiome in diet-metabolic signature relationship
#
#--------------------------------------------------------------------------------------------
#combine predicted score with mediation analysis
exp <- c("amed_11","ahei2010_11","dashav","pdi_11","hpdi_11","updi_11","edip11","edih11")
med <- c("amed_p","ahei_p","dash_p","pdi_p","hpdi_p","updi_p","edip_p","edih_p")
out <- c("amed2","ahei2","dash2","pdi2","hpdi2","updi2","edip2","edih2")

Predict <- Predict[,-c(1:8)]
names(Predict) <- med

df2 <- cbind(df,Predict[med])

#inverse-normal transformation
inormal <- function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

df2[,c(exp,med)] <- apply(df2[,c(exp,med)],2,inormal)

#mediation analysis
metab_dataset <- df2 %>% dplyr::select(c(med))

controlled_med<- metab_dataset%>%
  dplyr::select(med)%>%data.table::melt()%>%
  group_by(variable)%>%
  summarise(median=median(value,na.rm = T))

controlled_med <- as.vector(controlled_med$median)

#check covariates
ba_cov <- apply(df2[covar],2,summary)

#check exposures-missing
ba_exp <- apply(df2[exp],2,summary)

#check mediators-no missing values
ba_med <- apply(df2[med],2,summary)

#check outcomes-no misssing values
ba_out <- apply(df2[out],2,summary)

#impute missing continuous covariate using median value
medfunc <- function(x){
  missplace <- which(is.na(x))
  x[missplace] <- median(x, na.rm = T)/2
  x
}

catfunc <- function(x){
  missplace <- which(is.na(x))
  x[missplace] <- 0
  x
}

df2[,c("ageyr","bmi","alco","act","calor")] <- apply(df2[,c("ageyr","bmi","alco","act","calor")],2,medfunc)
df2[,c("stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6","stooltype_fec.7","stooltype_fec.8","antibio_12m_fec","colsc_2m_fec","probio_2m_fec","acid_2m_fec")] <- apply(df2[,c("stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6","stooltype_fec.7","stooltype_fec.8","antibio_12m_fec","colsc_2m_fec","probio_2m_fec","acid_2m_fec")],2,catfunc)
df2[which(is.na(df2$smoke)),"smoke"] <- 1

#run mediation model
cmest1 <- cmest(data=subset(df2,amed_11>-10),model="rb",outcome="amed2",event=NULL,exposure="amed_11",mediator="amed_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=200,yvar=1)
cmest2 <- cmest(data=subset(df2,ahei2010_11>-10),model="rb",outcome="ahei2",event=NULL,exposure="ahei2010_11",mediator="ahei_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=200,yvar=1) 
cmest3 <- cmest(data=subset(df2,dashav>-10),model="rb",outcome="dash2",event=NULL,exposure="dashav",mediator="dash_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=200,yvar=1) 
cmest4 <- cmest(data=subset(df2,pdi_11>-10),model="rb",outcome="pdi2",event=NULL,exposure="pdi_11",mediator="pdi_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=200,yvar=1) 
cmest5 <- cmest(data=subset(df2,hpdi_11>-10),model="rb",outcome="hpdi2",event=NULL,exposure="hpdi_11",mediator="hpdi_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=200,yvar=1) 
cmest6 <- cmest(data=subset(df2,updi_11>-10),model="rb",outcome="updi2",event=NULL,exposure="updi_11",mediator="updi_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=200,yvar=1) 
cmest7 <- cmest(data=subset(df2,edip11>-10),model="rb",outcome="edip2",event=NULL,exposure="edip11",mediator="edip_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=200,yvar=1) 
cmest8 <- cmest(data=subset(df2,edih11>-10),model="rb",outcome="edih2",event=NULL,exposure="edih11",mediator="edih_p",basec=covar,EMint=FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=200,yvar=1)

#get mediation results
res1 <- as.data.frame(as.matrix(summary(cmest1)$summarydf)); res1$mediator <- "amed"
res2 <- as.data.frame(as.matrix(summary(cmest2)$summarydf)); res2$mediator <- "ahei"
res3 <- as.data.frame(as.matrix(summary(cmest3)$summarydf)); res3$mediator <- "dash"
res4 <- as.data.frame(as.matrix(summary(cmest4)$summarydf)); res4$mediator <- "pdi"
res5 <- as.data.frame(as.matrix(summary(cmest5)$summarydf)); res5$mediator <- "hpdi"
res6 <- as.data.frame(as.matrix(summary(cmest6)$summarydf)); res6$mediator <- "updi"
res7 <- as.data.frame(as.matrix(summary(cmest7)$summarydf)); res7$mediator <- "edip"
res8 <- as.data.frame(as.matrix(summary(cmest8)$summarydf)); res8$mediator <- "edih"

res_med <- rbind(res1["pm",],res2["pm",],res3["pm",],res4["pm",],res5["pm",],res6["pm",],res7["pm",],res8["pm",])

#save mediation results
write.csv(res_med, file = "Figure7.2_MBS.csv")

#--------------------------------------------------------------------------------------------
#
#           Step3: test the hypothesis for doing mediation analysis
#
#--------------------------------------------------------------------------------------------
#define exposure, mediator, and outcome
cov <- covar

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
    fit <- lm(y~Score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec,data=df2)
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
    fit <- lm(y~Score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec,data=df2)
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
    fit <- lm(y~Score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec,data=df2)
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
write.csv(med_tot, file = "TableS15_MBS.csv")