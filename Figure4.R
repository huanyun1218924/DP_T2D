# --------------------------------------------------------------------------------
#Title: Dietary quality scores, metabolic signature, and cardiometabolic diseases
#Purpose: Association between diet/metabolic signatures and incident T2D
#Study: NHS/HPFS
#Path: /udd/nhhyu/DP_T2D/ProgramReview
#On: nantucket
#Programmer: Huan Yun (nhhyu)
#Date: 20230718
# --------------------------------------------------------------------------------

#load packages
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(coxme)

#set working dir
setwd("/udd/nhhyu/DP_T2D/FinalOutput")

#define functions for data analysis
int = function(var,dataframe) {
  y=dataframe[,var]
  inty=NA
  inty[!is.na(y)] = qnorm(rank(as.numeric(y[!is.na(y)]))/(length(as.numeric(y[!is.na(y)]))+1),mean=0,sd=1)
  inty
}

singlefun = function(exposure,outcome,dataframe){
  
  d = dataframe[,c(exposure,outcome,"ptime_diabetes","ageyr","studycaco","fast","aspirinuse","antihluse","smoking","actcat","fhxpremi","fhxdb","phxchol","phxhbp","BMIcata","caloravn","id")]
  names(d)[1:2] = c("Score","diabetes")
  
  # association
  res = rbind(
    
    coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Score + ageyr + strata(studycaco), data = d)))["Score",c(1,3,5)],
    coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Score + ageyr + strata(studycaco) + fast + aspirinuse + antihluse + as.factor(smoking) + actcat + fhxdb + phxchol + phxhbp + caloravn, data = d)))["Score",c(1,3,5)],
    coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Score + ageyr + strata(studycaco) + fast + aspirinuse + antihluse + as.factor(smoking) + actcat + fhxdb + phxchol + phxhbp + BMIcata + caloravn, data = d)))["Score",c(1,3,5)]
    
  )
  
  rownames(res) = c("Age-adj","MV","MV+BMI")
  colnames(res) = c("Est","sem","P")
  round(res,digit=100)
  
}

twofun = function(exposure1,exposure2,outcome,dataframe){
  
  d = dataframe[,c(exposure1,exposure2,outcome,"ptime_diabetes","ageyr","studycaco","fast","aspirinuse","antihluse","smoking","actcat","fhxpremi","fhxdb","phxchol","phxhbp","BMIcata","caloravn","id")]
  names(d)[1:3] = c("Score1","Score2","diabetes")
  
  # association
  fit1 = coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Score1 + Score2 + ageyr + strata(studycaco), data = d)))
  fit2 = coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Score1 + Score2 + ageyr + strata(studycaco) + fast + aspirinuse + antihluse + as.factor(smoking) + actcat + fhxdb + phxchol + phxhbp + caloravn, data = d)))
  fit3 = coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Score1 + Score2 + ageyr + strata(studycaco) + fast + aspirinuse + antihluse + as.factor(smoking) + actcat + fhxdb + phxchol + phxhbp + caloravn + BMIcata, data = d)))
  
  res = rbind( c(fit1["Score1",c(1,3,5)],fit1["Score2",c(1,3,5)]),
               c(fit2["Score1",c(1,3,5)],fit2["Score2",c(1,3,5)]),
               c(fit3["Score1",c(1,3,5)],fit3["Score2",c(1,3,5)]) )
  
  rownames(res) = c("Age-adj","MV","MV+BMI")
  colnames(res) = c("Score1_Est","Score1_sem","Score1_P","Score2_Est","Score2_sem","Score2_P")
  round(res,digit=100)
  
}

#load data
load("NHSHPFS_Final_T2D.RData")
t2d <- t2d4 %>% as.data.frame()
t2d$caloravn <- t2d$caloravn.x

covar <- c("ageyr","study","studycaco","fast","aspirinuse","mnpmh","smoking","actcat","fhxpremi","fhxdb",
           "phxchol", "phxhbp", "dbbase", "BMIcata","BMIcont","actcont",
           "amed","ahei","dash","opdi","hpdi","updi","edip","edih","diabetes","CVD","CVDdeath","CHD","Stroke",
           "amed_av","ahei2010_av","dashav","pdi_av","hpdi_av","updi_av","edipav","edihav","caloravn",
           "id","ptime_CVD","ptime_diabetes","ptime_CVDdeath","irt90","dtdxcvd","cabgbase","hrtbase",
           "strbase","canbase","dtdxdb2","dtdeath","lastq","cohort")

t2d$amed = int("amed2",t2d); t2d$amed_av = int("amed_av",t2d)
t2d$ahei = int("ahei2",t2d); t2d$ahei2010_av = int("ahei2010_av",t2d)
t2d$dash = int("dash2",t2d); t2d$dashav = int("dashav",t2d)
t2d$opdi = int("pdi2",t2d); t2d$pdi_av = int("pdi_av",t2d)
t2d$hpdi = int("hpdi2",t2d); t2d$hpdi_av = int("hpdi_av",t2d)
t2d$updi = int("updi2",t2d); t2d$updi_av = int("updi_av",t2d)
t2d$edip = int("edip2",t2d); t2d$edipav = int("edipav",t2d)
t2d$edih = int("edih2",t2d); t2d$edihav = int("edihav",t2d)

#diabetes
d11 <- singlefun("amed_av","diabetes",t2d)        
d21 <- singlefun("ahei2010_av","diabetes",t2d)     
d31 <- singlefun("dashav","diabetes",t2d)         
d41 <- singlefun("pdi_av","diabetes",t2d)          
d51 <- singlefun("hpdi_av","diabetes",t2d)         
d61 <- singlefun("updi_av","diabetes",t2d)         
d71 <- singlefun("edipav","diabetes",t2d)          
d81 <- singlefun("edihav","diabetes",t2d)        

d12 <- singlefun("amed2","diabetes",t2d)           
d22 <- singlefun("ahei2","diabetes",t2d)          
d32 <- singlefun("dash2","diabetes",t2d)          
d42 <- singlefun("pdi2","diabetes",t2d)           
d52 <- singlefun("hpdi2","diabetes",t2d)            
d62 <- singlefun("updi2","diabetes",t2d)           
d72 <- singlefun("edip2","diabetes",t2d)            
d82 <- singlefun("edih2","diabetes",t2d)           

d13 <- twofun("amed_av","amed2","diabetes",t2d)    
d23 <- twofun("ahei2010_av","ahei2","diabetes",t2d) 
d33 <- twofun("dashav","dash2","diabetes",t2d)      
d43 <- twofun("pdi_av","pdi2","diabetes",t2d)     
d53 <- twofun("hpdi_av","hpdi2","diabetes",t2d)     
d63 <- twofun("updi_av","updi2","diabetes",t2d)     
d73 <- twofun("edipav","edip2","diabetes",t2d)      
d83 <- twofun("edihav","edih2","diabetes",t2d)     

m <- rbind(cbind(d11,d12),cbind(d21,d22),cbind(d31,d32),cbind(d41,d42),
           cbind(d51,d52),cbind(d61,d62),cbind(d71,d72),cbind(d81,d82))

mm <- rbind(d13[3,],d23[3,],d33[3,],d43[3,],d53[3,],d63[3,],d73[3,],d83[3,])

#save results
list_of_sheets <- list("sheet_a" = m,"sheet_b" = mm) 
openxlsx::write.xlsx(list_of_sheets, "Figure4.xlsx",rowNames=T)