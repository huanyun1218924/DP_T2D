# --------------------------------------------------------------------------------
#Title:Dietary quality scores, metabolic signature, and cardiometabolic diseases
#Purpose: Explore the interaction effect between diet and microbiome on plasma metabolic signatures
#Note: in this part, we will analyse both overall and individual species. First, we will compare the contribution with and without interaction term to the variation of metabolic signatures; Second, we will perform LRT test to get p for interaction between diet and gut microbiome (stratified by continuous, median, carrier and non-carrier) on metabolic signatures
#Study: MLVS, MBS
#Programmer: Huan Yun
#Date: 20240119
# --------------------------------------------------------------------------------

#load packages
library(data.table)
library(dplyr)
library(readxl)
library(openxlsx)
library(cvTools)
library(glmnet)
library(ggplot2)
library(lmtest)
library(gridExtra)
library(cowplot)

#set working dir
setwd("/udd/nhhyu/DP_T2D/FinalOutput")

#--------------------------------------------------------------------------------------------
#
#           step1: prepare mlvs data
#
#--------------------------------------------------------------------------------------------
#load mlvs data
load("ToBeUsed.Taxon.RData")

#load filtered results
taxa_c <- read_excel("/udd/nhhyu/DP_T2D/results/Feature-wide association results.xlsx", sheet = "sheet_a") %>% as.data.frame()
microvar <- subset(taxa_c, type == "species")$microName
length(microvar) #151

Bug_MLVS <- microvar

#get annotation file for mlvs
Anno_MLVS <- taxa_c[,1:15]

#get name of dietary score and signature
diet <- c("amed_av","ahei_av","dash_av","pdi_av","hpdi_av","updi_av","edip_av","edih_av")
ms <- c("amed_c","ahei_c","dash_c","pdi_c","hpdi_c","updi_c","edip_c","edih_c")

#unify dietary score and metabolic signature
names(ToBeUsed.Taxon)[match(c(diet,ms),names(ToBeUsed.Taxon))] <- c("amed_score","ahei_score","dash_score","pdi_score","hpdi_score","updi_score","edip_score","edih_score","amed","ahei","dash","pdi","hpdi","updi","edip","edih")

#standarise microbiome data (TSS transformation and then scale)
Taxon_MLVS_Tss <- ToBeUsed.Taxon
Taxon_MLVS_Std <- ToBeUsed.Taxon
Taxon_MLVS_Std[microvar] <- apply(Taxon_MLVS_Std[microvar], 2, scale)

#--------------------------------------------------------------------------------------------
#
#           step2: prepare mbs data
#
#--------------------------------------------------------------------------------------------
#read data
load("/udd/nhhyu/DP_T2D/data/Microbiome/ToBeUsed.Taxon.MBS.RData")
Taxon_MBS <- as.data.frame(cov_avg)

#remove raw diet data
Taxon_MBS <- Taxon_MBS[,-c(796:800)]

#read cohort diet data from cohort
nhs2 <- fread("/udd/nhhyu/DP_T2D/data/Phenotype/NHS2/Diet_NHS2_0711.csv")

#merge cohort diet data with MBS
Taxon_MBS <- left_join(Taxon_MBS, nhs2[,c("id","amed_11","ahei2010_11","dashav","pdi_11","hpdi_11","updi_11","edip11","edih11")], by = "id")

#get filtered species name
anno_f <- fread("/udd/nhhyu/DP_T2D/data/Microbiome/Filtered_Anno_MBS.csv")
Bug_MBS <- anno_f$labname 
length(Bug_MBS) #153

#get annotation file for mbs
Anno_MBS <- anno_f

#get name of dietary score and signature
diet <- c("amed_11","ahei2010_11","dashav","pdi_11","hpdi_11","updi_11","edip11","edih11")
ms <- c("amed2","ahei2","dash2","pdi2","hpdi2","updi2","edip2","edih2")

#unify dietary score and metabolic signature
names(Taxon_MBS)[match(c(diet,ms),names(Taxon_MBS))] <- c("amed_score","ahei_score","dash_score","pdi_score","hpdi_score","updi_score","edip_score","edih_score","amed","ahei","dash","pdi","hpdi","updi","edip","edih")

#standarise microbiome data (TSS+scale)
Taxon_MBS_Tss <- Taxon_MBS
Taxon_MBS_Std <- Taxon_MBS
Taxon_MBS_Std[Bug_MBS] <- apply(Taxon_MBS_Std[Bug_MBS], 2, scale)

#--------------------------------------------------------------------------------------------
#
#           step3: run elastic net regression to get overall microbial score--tss+scale
#
#--------------------------------------------------------------------------------------------
#define outcome
outvar <- c("amed","ahei","dash","pdi","hpdi","updi","edip","edih")

###MLVS
#define dataset
df <- as.data.frame(Taxon_MLVS_Std)
df <- subset(df,amed>0)

#define microbiome
microvar

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

names(Predict) <- c("amed_p","ahei_p","dash_p","pdi_p","hpdi_p","updi_p","edip_p","edih_p")
Predict$id <- df$id

#combine predicted microbial score with tss dataset
Taxon_MLVS_Tss <- left_join(Taxon_MLVS_Tss,Predict,by="id")

###MBS
#define dataset
df <- as.data.frame(Taxon_MBS_Std)
df <- subset(df,amed>0)

#define microbiome
microvar <- Bug_MBS

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

names(Predict) <- c("amed_p","ahei_p","dash_p","pdi_p","hpdi_p","updi_p","edip_p","edih_p")
Predict$id <- df$id

#combine predicted microbial score with tss dataset
Taxon_MBS_Tss <- left_join(Taxon_MBS_Tss,Predict,by="id")

#--------------------------------------------------------------------------------------------
#
#           step4: get significant species associated with metabolic signatures
#
#--------------------------------------------------------------------------------------------
#read list in which we listed the taxa selected by elastic net regression
names(Anno_MLVS)[1] <- "Variable"
names(Anno_MLVS)[10] <- "Species"

names(Anno_MBS)[1] <- "Variable"
names(Anno_MBS)[8] <- "Species"

sig_list_M1 <- fread("/udd/nhhyu/DP_T2D/data/Microbiome/List_bug_elas_M1.csv")
names(sig_list_M1)[9] <- "Species"
sig_list_M1 <- rbind(left_join(sig_list_M1[which(sig_list_M1$Study=="MLVS"),],Anno_MLVS[,c("Variable","Species")],by="Species"),
                     left_join(sig_list_M1[which(sig_list_M1$Study=="MBS"),],Anno_MBS[,c("Variable","Species")],by="Species"))


#get significant species
name <- c("AMED","AHEI-2010","DASH","PDI","hPDI","uPDI","EDIP","EDIH")

EN_MLVS <- list()
EN_MBS <- list()

for (i in name){
  a0 <- sig_list_M1[which(sig_list_M1$Study=="MLVS"&sig_list_M1$DP==i&sig_list_M1$FDR<=0.05),]$Variable
  b0 <- sig_list_M1[which(sig_list_M1$Study=="MBS"&sig_list_M1$DP==i&sig_list_M1$FDR<=0.05),]$Variable
  EN_MLVS <- append(EN_MLVS,list(a0))
  EN_MBS <- append(EN_MBS,list(b0))
}

pre <- c("amed_p","ahei_p","dash_p","pdi_p","hpdi_p","updi_p","edip_p","edih_p")

#MLVS: for elas, we get ast, median, carrier variable; for microbial score, we get scale and median
for (i in c(EN_MLVS[[1]],EN_MLVS[[2]],EN_MLVS[[3]],EN_MLVS[[4]],EN_MLVS[[5]],EN_MLVS[[6]],EN_MLVS[[7]],EN_MLVS[[8]])){
  Taxon_MLVS_Tss$new = asin(sqrt(Taxon_MLVS_Tss[,i]))
  names(Taxon_MLVS_Tss)[which(names(Taxon_MLVS_Tss)=="new")] = paste(i,"_std",sep='')
  Taxon_MLVS_Tss$new2 = ifelse(Taxon_MLVS_Tss[,i]>median(Taxon_MLVS_Tss[,i]),"1","0")
  names(Taxon_MLVS_Tss)[which(names(Taxon_MLVS_Tss)=="new2")] = paste(i,"_med",sep='')
  Taxon_MLVS_Tss$new3 = ifelse(Taxon_MLVS_Tss[,i]>0,"1","0")
  names(Taxon_MLVS_Tss)[which(names(Taxon_MLVS_Tss)=="new3")] = paste(i,"_int",sep='')
}

for (i in pre){
  Taxon_MLVS_Tss$new = scale(Taxon_MLVS_Tss[,i])
  names(Taxon_MLVS_Tss)[which(names(Taxon_MLVS_Tss)=="new")] = paste(i,"_std",sep='')
  Taxon_MLVS_Tss$new2 = ifelse(Taxon_MLVS_Tss[,i]>median(Taxon_MLVS_Tss[,i]),"1","0")
  names(Taxon_MLVS_Tss)[which(names(Taxon_MLVS_Tss)=="new2")] = paste(i,"_med",sep='')
}

Int_MLVS <- list(c("amed_p_std","amed_p_med",paste(EN_MLVS[[1]],"_std",sep=''),paste(EN_MLVS[[1]],"_med",sep=''),paste(EN_MLVS[[1]],"_int",sep='')),
                 c("ahei_p_std","ahei_p_med",paste(EN_MLVS[[2]],"_std",sep=''),paste(EN_MLVS[[2]],"_med",sep=''),paste(EN_MLVS[[2]],"_int",sep='')),
                 c("dash_p_std","dash_p_med",paste(EN_MLVS[[3]],"_std",sep=''),paste(EN_MLVS[[3]],"_med",sep=''),paste(EN_MLVS[[3]],"_int",sep='')),
                 c("pdi_p_std","pdi_p_med",paste(EN_MLVS[[4]],"_std",sep=''),paste(EN_MLVS[[4]],"_med",sep=''),paste(EN_MLVS[[4]],"_int",sep='')),
                 c("hpdi_p_std","hpdi_p_med",paste(EN_MLVS[[5]],"_std",sep=''),paste(EN_MLVS[[5]],"_med",sep=''),paste(EN_MLVS[[5]],"_int",sep='')),
                 c("updi_p_std","updi_p_med",paste(EN_MLVS[[6]],"_std",sep=''),paste(EN_MLVS[[6]],"_med",sep=''),paste(EN_MLVS[[6]],"_int",sep='')),
                 c("edip_p_std","edip_p_med",paste(EN_MLVS[[7]],"_std",sep=''),paste(EN_MLVS[[7]],"_med",sep=''),paste(EN_MLVS[[7]],"_int",sep='')),
                 c("edih_p_std","edih_p_med",paste(EN_MLVS[[8]],"_std",sep=''),paste(EN_MLVS[[8]],"_med",sep=''),paste(EN_MLVS[[8]],"_int",sep='')))

Int_MLVS2 <- list(c("amed_p_med",paste(EN_MLVS[[1]],"_med",sep=''),paste(EN_MLVS[[1]],"_int",sep='')),
                  c("ahei_p_med",paste(EN_MLVS[[2]],"_med",sep=''),paste(EN_MLVS[[2]],"_int",sep='')),
                  c("dash_p_med",paste(EN_MLVS[[3]],"_med",sep=''),paste(EN_MLVS[[3]],"_int",sep='')),
                  c("pdi_p_med",paste(EN_MLVS[[4]],"_med",sep=''),paste(EN_MLVS[[4]],"_int",sep='')),
                  c("hpdi_p_med",paste(EN_MLVS[[5]],"_med",sep=''),paste(EN_MLVS[[5]],"_int",sep='')),
                  c("updi_p_med",paste(EN_MLVS[[6]],"_med",sep=''),paste(EN_MLVS[[6]],"_int",sep='')),
                  c("edip_p_med",paste(EN_MLVS[[7]],"_med",sep=''),paste(EN_MLVS[[7]],"_int",sep='')),
                  c("edih_p_med",paste(EN_MLVS[[8]],"_med",sep=''),paste(EN_MLVS[[8]],"_int",sep='')))

#for mbs, we get scale, median, carrier variable
for (i in c(EN_MBS[[1]],EN_MBS[[2]],EN_MBS[[3]],EN_MBS[[4]],EN_MBS[[5]],EN_MBS[[6]],EN_MBS[[7]],EN_MBS[[8]])){
  Taxon_MBS_Tss$new = asin(sqrt(Taxon_MBS_Tss[,i]))
  names(Taxon_MBS_Tss)[which(names(Taxon_MBS_Tss)=="new")] = paste(i,"_std",sep='')
  Taxon_MBS_Tss$new2 = ifelse(Taxon_MBS_Tss[,i]>median(Taxon_MBS_Tss[,i]),"1","0")
  names(Taxon_MBS_Tss)[which(names(Taxon_MBS_Tss)=="new2")] = paste(i,"_med",sep='')
  Taxon_MBS_Tss$new3 = ifelse(Taxon_MBS_Tss[,i]>0,"1","0")
  names(Taxon_MBS_Tss)[which(names(Taxon_MBS_Tss)=="new3")] = paste(i,"_int",sep='')
}

for (i in pre){
  Taxon_MBS_Tss$new = scale(Taxon_MBS_Tss[,i])
  names(Taxon_MBS_Tss)[which(names(Taxon_MBS_Tss)=="new")] = paste(i,"_std",sep='')
  Taxon_MBS_Tss$new2 = ifelse(Taxon_MBS_Tss[,i]>median(Taxon_MBS_Tss[,i],na.rm=T),"1","0")
  names(Taxon_MBS_Tss)[which(names(Taxon_MBS_Tss)=="new2")] = paste(i,"_med",sep='')
}

Int_MBS <- list(c("amed_p_std","amed_p_med"),c("ahei_p_std","ahei_p_med"), c("dash_p_std","dash_p_med"),c("pdi_p_std","pdi_p_med"),
                c("hpdi_p_std","hpdi_p_med"),c("updi_p_std","updi_p_med"),c("edip_p_std","edip_p_med"),c("edih_p_std","edih_p_med"))

Int_MBS2 <- list(c("amed_p_med"),c("ahei_p_med"), c("dash_p_med"),c("pdi_p_med"),c("hpdi_p_med"),c("updi_p_med"),c("edip_p_med"),c("edih_p_med"))

#save all data together: mlvs+mbs+sol (taxa+annotation)
save(Taxon_MLVS_Tss,Anno_MLVS,Bug_MLVS,EN_MLVS,Int_MLVS, #taxa,annotation file,taxa list in MLVS
     Taxon_MBS_Tss,Anno_MBS,Bug_MBS,EN_MBS,Int_MBS,      #taxa,annotation file,taxa list in MBS
     file="Data_For_Use_Microbiome.RData")

#--------------------------------------------------------------------------------------------
#
#           step6: calculate p for interaction
#
#--------------------------------------------------------------------------------------------
diet <- c("amed_score","ahei_score","dash_score","pdi_score","hpdi_score","updi_score","edip_score","edih_score")
ms <- c("amed","ahei","dash","pdi","hpdi","updi","edip","edih")

###mlvs
data = Taxon_MLVS_Tss
data[,c(diet,ms)] <- apply(data[,c(diet,ms)],2,scale)

#amed
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MLVS[[1]]){
  data$score=data[,i]
  fit1 <- glm(amed~score+amed_score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  fit2 <- glm(amed~amed_score*score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
amed_int <- res

#ahei in mlvs
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MLVS[[2]]){
  data$score=data[,i]
  fit1 <- glm(ahei~score+ahei_score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  fit2 <- glm(ahei~ahei_score*score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
ahei_int <- res

#dash in mlvs
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MLVS[[3]]){
  data$score=data[,i]
  fit1 <- glm(dash~score+dash_score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  fit2 <- glm(dash~dash_score*score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
dash_int <- res

#pdi in mlvs
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MLVS[[4]]){
  data$score=data[,i]
  fit1 <- glm(pdi~score+pdi_score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  fit2 <- glm(pdi~pdi_score*score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
pdi_int <- res

#hpdi in mlvs
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MLVS[[5]]){
  data$score=data[,i]
  fit1 <- glm(hpdi~score+hpdi_score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  fit2 <- glm(hpdi~hpdi_score*score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
hpdi_int <- res

#updi in mlvs
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MLVS[[6]]){
  data$score=data[,i]
  fit1 <- glm(updi~score+updi_score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  fit2 <- glm(updi~updi_score*score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
updi_int <- res

#edip in mlvs
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MLVS[[7]]){
  data$score=data[,i]
  fit1 <- glm(edip~score+edip_score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  fit2 <- glm(edip~edip_score*score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
edip_int <- res

#edih in mlvs
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MLVS[[8]]){
  data$score=data[,i]
  fit1 <- glm(edih~score+edih_score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  fit2 <- glm(edih~edih_score*score+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
edih_int <- res

#combine mlvs results
Res_MLVS <- rbind(amed_int,ahei_int,dash_int,pdi_int,hpdi_int,updi_int,edip_int,edih_int)
Res_MLVS$Study <- "MLVS" 

###mbs
data = Taxon_MBS_Tss
data[,c(diet,ms)] <- apply(data[,c(diet,ms)],2,scale)

#amed
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MBS[[1]]){
  data$score=data[,i]
  fit1 <- glm(amed~score+amed_score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  fit2 <- glm(amed~amed_score*score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
amed_int <- res

#ahei in MBS
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MBS[[2]]){
  data$score=data[,i]
  fit1 <- glm(ahei~score+ahei_score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  fit2 <- glm(ahei~ahei_score*score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
ahei_int <- res

#dash in MBS
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MBS[[3]]){
  data$score=data[,i]
  fit1 <- glm(dash~score+dash_score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  fit2 <- glm(dash~dash_score*score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
dash_int <- res

#pdi in MBS
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MBS[[4]]){
  data$score=data[,i]
  fit1 <- glm(pdi~score+pdi_score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  fit2 <- glm(pdi~pdi_score*score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
pdi_int <- res

#hpdi in MBS
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MBS[[5]]){
  data$score=data[,i]
  fit1 <- glm(hpdi~score+hpdi_score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  fit2 <- glm(hpdi~hpdi_score*score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
hpdi_int <- res

#updi in MBS
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MBS[[6]]){
  data$score=data[,i]
  fit1 <- glm(updi~score+updi_score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  fit2 <- glm(updi~updi_score*score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
updi_int <- res

#edip in MBS
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MBS[[7]]){
  data$score=data[,i]
  fit1 <- glm(edip~score+edip_score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  fit2 <- glm(edip~edip_score*score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
edip_int <- res

#edih in MBS
res = data.frame(Est=NA,SE=NA,Pinter=NA)
for(i in Int_MBS[[8]]){
  data$score=data[,i]
  fit1 <- glm(edih~score+edih_score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  fit2 <- glm(edih~edih_score*score+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data=data)
  res[i,"Est"] = coef(summary(fit2))[nrow(coef(summary(fit2))),1]
  res[i,"SE"] = coef(summary(fit2))[nrow(coef(summary(fit2))),2]
  res[i,"Pinter"] = lrtest(fit1,fit2)[2,5]
}

res <- res[-1,]
edih_int <- res

#combine MBS results
Res_MBS <- rbind(amed_int,ahei_int,dash_int,pdi_int,hpdi_int,updi_int,edip_int,edih_int)
Res_MBS$Study <- "MBS" 

#combine and save all results
Res_Int_All <- rbind(Res_MLVS,Res_MBS)
Res_Int_All2 <- subset(Res_Int_All,Pinter<=0.05)

#--------------------------------------------------------------------------------------------
#
#           step8: Interaction stratified by Median
#
#--------------------------------------------------------------------------------------------
#define function for stratified interaction analysis
intfunc1 = function(exposure1,exposure2,exposure3,outcome,dataframe){
  
  d = dataframe[,c(exposure1,exposure2,exposure3,outcome,"age_fec","bmi_bld","totMETs_paq","smoke_bld","probio_2m_fec","antibio_12m_fec","colsc_2m_fec","acid_2m_fec","stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6")]
  names(d)[1:4] = c("Score1","Score2","Score3","MS")
  
  # association
  res = rbind(
    
    coef(summary(glm(MS ~Score1*Score2+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data = d[which(d$Score3==0),])))["Score1:Score2",c(1,2,4)],
    coef(summary(glm(MS ~Score1*Score2+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data = d[which(d$Score3==1),])))["Score1:Score2",c(1,2,4)]
  )
  
  rownames(res) = c("Below","Above")
  colnames(res) = c("Est","sem","P")
  round(res,digit=12)
  
}

intfunc2 = function(exposure1,exposure2,exposure3,outcome,dataframe){
  
  d = dataframe[,c(exposure1,exposure2,exposure3,outcome,"ageyr","race","bmi","act","alco","smoke","calor","stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6","stooltype_fec.7","stooltype_fec.8","antibio_12m_fec","colsc_2m_fec","probio_2m_fec","acid_2m_fec")]
  names(d)[1:4] = c("Score1","Score2","Score3","MS")
  
  # association
  res = rbind(
    
    coef(summary(glm(MS ~Score1*Score2+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data = d[which(d$Score3==0),])))["Score1:Score2",c(1,2,4)],
    coef(summary(glm(MS ~Score1*Score2+ageyr+race+bmi+act+alco+smoke+calor+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6+stooltype_fec.7+stooltype_fec.8+antibio_12m_fec+colsc_2m_fec+probio_2m_fec+acid_2m_fec, data = d[which(d$Score3==1),])))["Score1:Score2",c(1,2,4)]
  )
  
  rownames(res) = c("Below","Above")
  colnames(res) = c("Est","sem","P")
  round(res,digit=12)
  
}

intfunc3 = function(exposure1,exposure2,exposure3,outcome,dataframe){
  
  d = dataframe[,c(exposure1,exposure2,exposure3,outcome,"AGE_V2","GENDER_V2","CIGARETTE_USE_V2","ALCOHOL_USE_V2","BMI_V2")]
  names(d)[1:4] = c("Score1","Score2","Score3","MS")
  
  # association
  res = rbind(
    
    coef(summary(glm(MS ~Score1*Score2+AGE_V2+GENDER_V2+CIGARETTE_USE_V2+ALCOHOL_USE_V2+BMI_V2, data = d[which(d$Score3==0),])))["Score1:Score2",c(1,2,4)],
    coef(summary(glm(MS ~Score1*Score2+AGE_V2+GENDER_V2+CIGARETTE_USE_V2+ALCOHOL_USE_V2+BMI_V2, data = d[which(d$Score3==1),])))["Score1:Score2",c(1,2,4)]
  )
  
  rownames(res) = c("Below","Above")
  colnames(res) = c("Est","sem","P")
  round(res,digit=12)
  
}

###mlvs
data = Taxon_MLVS_Tss
data[,c(diet,ms)] <- apply(data[,c(diet,ms)],2,scale)

t1 <- intfunc1("amed_score","amed_p_std","amed_p_med","amed",data) %>% as.data.frame(); t1$Ifmed = rownames(t1); t1$Trait="AMED"
t2 <- intfunc1("ahei_score","ahei_p_std","ahei_p_med","ahei",data) %>% as.data.frame(); t2$Ifmed = rownames(t2); t2$Trait="AHEI-2010"
t3 <- intfunc1("dash_score","dash_p_std","dash_p_med","dash",data) %>% as.data.frame(); t3$Ifmed = rownames(t3); t3$Trait="DASH"
t4 <- intfunc1("pdi_score","pdi_p_std","pdi_p_med","pdi",data) %>% as.data.frame(); t4$Ifmed = rownames(t4); t4$Trait="PDI"
t5 <- intfunc1("hpdi_score","hpdi_p_std","hpdi_p_med","hpdi",data) %>% as.data.frame(); t5$Ifmed = rownames(t5); t5$Trait="hPDI"
t6 <- intfunc1("updi_score","updi_p_std","updi_p_med","updi",data) %>% as.data.frame(); t6$Ifmed = rownames(t6); t6$Trait="uPDI"
t7 <- intfunc1("edip_score","edip_p_std","edip_p_med","edip",data) %>% as.data.frame(); t7$Ifmed = rownames(t7); t7$Trait="EDIP"
t8 <- intfunc1("edih_score","edih_p_std","edih_p_med","edih",data) %>% as.data.frame(); t8$Ifmed = rownames(t8); t8$Trait="EDIH"

Int_Med_MLVS <- rbind(t1,t2,t3,t4,t5,t6,t7,t8)
Int_Med_MLVS$Study = "MLVS"

###mbs
data = Taxon_MBS_Tss
data[,c(diet,ms)] <- apply(data[,c(diet,ms)],2,scale)

t1 <- intfunc2("amed_score","amed_p_std","amed_p_med","amed",data) %>% as.data.frame(); t1$Ifmed = rownames(t1); t1$Trait="AMED"
t2 <- intfunc2("ahei_score","ahei_p_std","ahei_p_med","ahei",data) %>% as.data.frame(); t2$Ifmed = rownames(t2); t2$Trait="AHEI-2010"
t3 <- intfunc2("dash_score","dash_p_std","dash_p_med","dash",data) %>% as.data.frame(); t3$Ifmed = rownames(t3); t3$Trait="DASH"
t4 <- intfunc2("pdi_score","pdi_p_std","pdi_p_med","pdi",data) %>% as.data.frame(); t4$Ifmed = rownames(t4); t4$Trait="PDI"
t5 <- intfunc2("hpdi_score","hpdi_p_std","hpdi_p_med","hpdi",data) %>% as.data.frame(); t5$Ifmed = rownames(t5); t5$Trait="hPDI"
t6 <- intfunc2("updi_score","updi_p_std","updi_p_med","updi",data) %>% as.data.frame(); t6$Ifmed = rownames(t6); t6$Trait="uPDI"
t7 <- intfunc2("edip_score","edip_p_std","edip_p_med","edip",data) %>% as.data.frame(); t7$Ifmed = rownames(t7); t7$Trait="EDIP"
t8 <- intfunc2("edih_score","edih_p_std","edih_p_med","edih",data) %>% as.data.frame(); t8$Ifmed = rownames(t8); t8$Trait="EDIH"

Int_Med_MBS <- rbind(t1,t2,t3,t4,t5,t6,t7,t8)
Int_Med_MBS$Study = "MBS"

#combine all results
Res_Int_Med = rbind(Int_Med_MLVS,Int_Med_MBS)

#--------------------------------------------------------------------------------------------
#
#           step9: Annotation 
#
#--------------------------------------------------------------------------------------------
#read results
Res_Int_All <- fread("Res_Int_All.csv") %>% as.data.frame()
Res_Int_All <- Res_Int_All[,-1]

names(Res_Int_All)[1] <- "Variable"
Res_Int_All$Raw_name <- Res_Int_All$Variable

#get the name of dietary pattern
Res_Int_All[c(1:23,336:337,128:159),"Trait"] <- "AMED"
Res_Int_All[c(24:46,338:339,160:199),"Trait"] <- "AHEI-2010"
Res_Int_All[c(47:60,340:341,200:233),"Trait"] <- "DASH"
Res_Int_All[c(61:74,342:343,234:243),"Trait"] <- "PDI"
Res_Int_All[c(75:85,344:345,244:271),"Trait"] <- "hPDI"
Res_Int_All[c(86:108,346:347,272:293),"Trait"] <- "uPDI"
Res_Int_All[c(109:116,348:349,294:319),"Trait"] <- "EDIP"
Res_Int_All[c(117:127,350:351,320:335),"Trait"] <- "EDIH"

#merge with annotation file in each study
var <- c("Variable","kingdom","phylum","class","order","family","genus","Species")

Anno_MLVS <- as.data.frame(Anno_MLVS)
Anno_MBS <- as.data.frame(Anno_MBS)

names(Anno_MBS)[2] <- "kingdom"

Res_Int_All <- rbind(left_join(Res_Int_All[which(Res_Int_All$Study=="MLVS"&Res_Int_All$Trait=="AMED"),],Anno_MLVS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MLVS"&Res_Int_All$Trait=="AHEI-2010"),],Anno_MLVS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MLVS"&Res_Int_All$Trait=="DASH"),],Anno_MLVS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MLVS"&Res_Int_All$Trait=="PDI"),],Anno_MLVS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MLVS"&Res_Int_All$Trait=="hPDI"),],Anno_MLVS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MLVS"&Res_Int_All$Trait=="uPDI"),],Anno_MLVS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MLVS"&Res_Int_All$Trait=="EDIP"),],Anno_MLVS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MLVS"&Res_Int_All$Trait=="EDIH"),],Anno_MLVS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MBS"&Res_Int_All$Trait=="AMED"),],Anno_MBS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MBS"&Res_Int_All$Trait=="AHEI-2010"),],Anno_MBS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MBS"&Res_Int_All$Trait=="DASH"),],Anno_MBS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MBS"&Res_Int_All$Trait=="PDI"),],Anno_MBS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MBS"&Res_Int_All$Trait=="hPDI"),],Anno_MBS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MBS"&Res_Int_All$Trait=="uPDI"),],Anno_MBS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MBS"&Res_Int_All$Trait=="EDIP"),],Anno_MBS[var],by="Variable"),
                     left_join(Res_Int_All[which(Res_Int_All$Study=="MBS"&Res_Int_All$Trait=="EDIH"),],Anno_MBS[var],by="Variable"))

#get significant species
Res_Int_Sig <- subset(Res_Int_All,Pinter<=0.05)

#save results
write.csv(Res_Int_Sig,file="Figure7.3.csv")

#remove duplicate columns
Taxon_MLVS_Tss <- Taxon_MLVS_Tss[,!duplicated(names(Taxon_MLVS_Tss))]
Taxon_MBS_Tss <- Taxon_MBS_Tss[,!duplicated(names(Taxon_MBS_Tss))]

#mlvs
p1 <- ggplot(Taxon_MLVS_Tss,aes(x=ahei_score, y=ahei,color=ahei_p_med))+geom_point(size=3,alpha=0.7)+geom_smooth(method="lm", fill=NA)+scale_color_manual(labels=c("Below median","Above median"),values=c('orange','steelblue'))+scale_y_continuous("")+scale_x_continuous("")+theme_classic()+theme(panel.grid=element_blank(),legend.title = element_blank(),legend.position = "bottom")
p2 <- ggplot(Taxon_MLVS_Tss,aes(x=dash_score, y=dash,color=dash_p_med))+geom_point(size=3,alpha=0.7)+geom_smooth(method="lm", fill=NA)+scale_color_manual(labels=c("Below median","Above median"),values=c('orange','steelblue'))+scale_y_continuous("")+scale_x_continuous("")+theme_classic()+theme(panel.grid=element_blank(),legend.title = element_blank(),legend.position = "bottom")
p3 <- ggplot(Taxon_MLVS_Tss,aes(x=hpdi_score, y=hpdi,color=s_45851_med))+geom_point(size=3,alpha=0.7)+geom_smooth(method="lm", fill=NA)+scale_color_manual(labels=c("Below median","Above median"),values=c('orange','steelblue'))+scale_y_continuous("")+scale_x_continuous("")+theme_classic()+theme(panel.grid=element_blank(),legend.title = element_blank(),legend.position = "bottom")
p4 <- ggplot(Taxon_MLVS_Tss,aes(x=updi_score, y=updi,color=s_457402_med))+geom_point(size=3,alpha=0.7)+geom_smooth(method="lm", fill=NA)+scale_color_manual(labels=c("Below median","Above median"),values=c('orange','steelblue'))+scale_y_continuous("")+scale_x_continuous("")+theme_classic()+theme(panel.grid=element_blank(),legend.title = element_blank(),legend.position = "bottom")

#mbs
p5 <- ggplot(Taxon_MBS_Tss,aes(x=pdi_score, y=pdi,color=pdi_p_med))+geom_point(size=3,alpha=0.7)+geom_smooth(method="lm", fill=NA)+scale_color_manual(labels=c("Below median","Above median"),values=c('orange','steelblue'))+scale_y_continuous("")+scale_x_continuous("")+theme_classic()+theme(panel.grid=element_blank(),legend.title = element_blank(),legend.position = "bottom")
p6 <- ggplot(Taxon_MBS_Tss,aes(x=edip_score, y=edip,color=edip_p_med))+geom_point(size=3,alpha=0.7)+geom_smooth(method="lm", fill=NA)+scale_color_manual(labels=c("Below median","Above median"),values=c('orange','steelblue'))+scale_y_continuous("")+scale_x_continuous("")+theme_classic()+theme(panel.grid=element_blank(),legend.title = element_blank(),legend.position = "bottom")

pdf("FigureS17.pdf",height=16,width=12)
plot_grid(p1,p2,p3,p4,p5,p6,ncol = 4)
dev.off()