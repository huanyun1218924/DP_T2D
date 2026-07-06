library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(coxme)
library(metafor)
library(circlize)
library(RColorBrewer)
library(tidyr)
library(ComplexHeatmap)
library(gridExtra)

#define functions to be used
inormal <- function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

cat = function(prs,cutofftype,data) { # tertile, decile, top10, top20, topbut10, topbut20
  
  dat = data[c(prs)]
  names(dat) = c("prs")
  dat$cats = 0
  
  y.all = dat$prs
  
  
  if (cutofftype=="tertile") {
    dat$cats[ which( dat$prs>quantile(y.all, 0.33, na.rm=TRUE) ) ] = 1
    dat$cats[ which( dat$prs>quantile(y.all, 0.67, na.rm=TRUE) ) ] = 2
    
  } else if (cutofftype=="quartile") {
    dat$cats[ which( dat$prs>quantile(y.all, 0.25, na.rm=TRUE) ) ] = 1
    dat$cats[ which( dat$prs>quantile(y.all, 0.50, na.rm=TRUE) ) ] = 2
    dat$cats[ which( dat$prs>quantile(y.all, 0.75, na.rm=TRUE) ) ] = 3
    
  }else if (cutofftype=="quintile") {
    dat$cats[ which( dat$prs>quantile(y.all, 0.2, na.rm=TRUE) ) ] = 1
    dat$cats[ which( dat$prs>quantile(y.all, 0.4, na.rm=TRUE) ) ] = 2
    dat$cats[ which( dat$prs>quantile(y.all, 0.6, na.rm=TRUE) ) ] = 3
    dat$cats[ which( dat$prs>quantile(y.all, 0.8, na.rm=TRUE) ) ] = 4
    
  } else if (cutofftype=="decile") {
    dat$cats[ which( dat$prs>quantile(y.all, 0.1, na.rm=TRUE) ) ] = 1
    dat$cats[ which( dat$prs>quantile(y.all, 0.2, na.rm=TRUE) ) ] = 2
    dat$cats[ which( dat$prs>quantile(y.all, 0.3, na.rm=TRUE) ) ] = 3
    dat$cats[ which( dat$prs>quantile(y.all, 0.4, na.rm=TRUE) ) ] = 4
    dat$cats[ which( dat$prs>quantile(y.all, 0.5, na.rm=TRUE) ) ] = 5
    dat$cats[ which( dat$prs>quantile(y.all, 0.6, na.rm=TRUE) ) ] = 6
    dat$cats[ which( dat$prs>quantile(y.all, 0.7, na.rm=TRUE) ) ] = 7
    dat$cats[ which( dat$prs>quantile(y.all, 0.8, na.rm=TRUE) ) ] = 8
    dat$cats[ which( dat$prs>quantile(y.all, 0.9, na.rm=TRUE) ) ] = 9
    
  } else if (cutofftype=="top10" & methods=="all") {
    dat$cats[ which( dat$prs>quantile(y.all, 0.9, na.rm=TRUE) ) ] = 1
    
  } else if (cutofftype=="top20" & methods=="all") {
    dat$cats[ which( dat$prs>quantile(y.all, 0.8, na.rm=TRUE) ) ] = 1
    
    
  } else if (cutofftype=="topbut10" & methods=="all") {
    dat$cats[ which( dat$prs<quantile(y.all, 0.1, na.rm=TRUE) ) ] = -1
    dat$cats[ which( dat$prs>quantile(y.all, 0.9, na.rm=TRUE) ) ] = 1
    
  }  else if (cutofftype=="topbut20" & methods=="all") {
    dat$cats[ which( dat$prs<quantile(y.all, 0.2, na.rm=TRUE) ) ] = -1
    dat$cats[ which( dat$prs>quantile(y.all, 0.8, na.rm=TRUE) ) ] = 1
    
  } 
  
  dat$cats[is.na(dat$prs)]=NA
  return(dat$cats)
  
}

singlefun1 = function(exposure,outcome,dataframe){
  
  d = dataframe[,c(exposure,outcome,"ptime_diabetes","ageyr","fast","aspirinuse","smoking","act","fhxpremi","fhxdb","phxchol","phxhbp","BMIcont","id")]
  names(d)[1:2] = c("Score","diabetes")
  
  # association
  res = rbind(
    
    coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Score + ageyr, data = d)))["Score",c(1,3,5)],
    coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Score + ageyr+ fast + aspirinuse + as.factor(smoking) + act + fhxdb + phxchol + phxhbp, data = d)))["Score",c(1,3,5)],
    coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ Score + ageyr+ fast + aspirinuse + as.factor(smoking) + act + fhxdb + phxchol + phxhbp + BMIcont, data = d)))["Score",c(1,3,5)]
    
  )
  
  rownames(res) = c("Age-adj","MV","MV+BMI")
  colnames(res) = c("Est","sem","P")
  round(res,digit=100)
  
}

singlefun2 = function(exposure,outcome,dataframe){
  
  d = dataframe[,c(exposure,outcome,"ptime_diabetes","ageyr","fast","aspirinuse","smoking","act","fhxpremi","fhxdb","phxchol","phxhbp","BMIcont","id")]
  names(d)[1:2] = c("Score","diabetes")
  
  # association
  res = rbind(
    
    coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ relevel(as.factor(Score),ref='0') + ageyr, data = d)))[1:3,c(1,3,5)],
    coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ relevel(as.factor(Score),ref='0') + ageyr+ fast + aspirinuse + as.factor(smoking) + act + fhxdb + phxchol + phxhbp, data = d)))[1:3,c(1,3,5)],
    coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ relevel(as.factor(Score),ref='0') + ageyr+ fast + aspirinuse + as.factor(smoking) + act + fhxdb + phxchol + phxhbp + BMIcont, data = d)))[1:3,c(1,3,5)]
    
  )
  
  rownames(res) = rep(c("Age-adj","MV","MV+BMI"),each=3)
  colnames(res) = c("Est","sem","P")
  round(res,digit=100)
  
}

singlefun3 = function(exposure,outcome,dataframe){
  
  d = dataframe[,c(exposure,outcome,"ptime_diabetes","ageyr","fast","aspirinuse","smoking","act","fhxpremi","fhxdb","phxchol","phxhbp","BMIcont","id")]
  names(d)[1:2] = c("Score","diabetes")
  
  # association
  res = rbind(
    
    coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ as.numeric(Score) + ageyr, data = d)))[1,c(1,3,5)],
    coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ as.numeric(Score) + ageyr+ fast + aspirinuse + as.factor(smoking) + act + fhxdb + phxchol + phxhbp, data = d)))[1,c(1,3,5)],
    coef(summary(coxph(Surv(ptime_diabetes, diabetes) ~ as.numeric(Score) + ageyr+ fast + aspirinuse + as.factor(smoking) + act + fhxdb + phxchol + phxhbp + BMIcont, data = d)))[1,c(1,3,5)]
    
  )
  
  rownames(res) = rep(c("Age-adj","MV","MV+BMI"),each=1)
  colnames(res) = c("Est","sem","P")
  round(res,digit=100)
  
}

#--------------------------------------------------------------------------------------------
#
#           chunk1: PH assumption test 
#
#--------------------------------------------------------------------------------------------
#read sample data
load("t2d_sample.RData")

data_use <- t2d_sample

#get list of dietary patterns and metabolomic signatures
exposures <- c("amed_av","ahei_av","dash_av","pdi_av","hpdi_av","updi_av","edip_av","edih_av",
               "amed","ahei","dash","pdi","hpdi","updi","edip","edih")

var <- c("AMED","AHEI","DASH","PDI","hPDI","uPDI","EDIP","EDIH","AMED_Signature","AHEI_Signature","DASH_Signature","PDI_Signature","hPDI_Signature","uPDI_Signature","EDIP_Signature","EDIH_Signature")

covariates <-  c("ageyr","fast","aspirinuse","smoking","act","fhxpremi","fhxdb","phxchol","phxhbp","BMIcont")

#inverse-normal transformation for dietary patterns and metabolomic signatures
data_use[exposures] <- apply(data_use[exposures],2,inormal)

#plot PH assumption
plots<- list()

for (i in 1:16) {

  formula_text <- paste("Surv(ptime_diabetes, diabetes) ~",exposures[i],"+ ageyr + fast + aspirinuse + as.factor(smoking) + act + fhxdb + phxchol + phxhbp + BMIcont", sep='')
  fit <- coxph(as.formula(formula_text), data=data_use)
  ph <- cox.zph(fit,transform = "identity")
  ph_p <- ph$table[exposures[i],"p"]
  plot_data <- data.frame(time = ph$x,residual = ph$y[,exposures[i]])
 
  plots[i] <- ggplot(plot_data, aes(x=time, y=residual)) +
    geom_smooth(method = "loess", se = TRUE, color="#008B8B", fill="#AFEEEE",linetype = "solid",alpha=0.5) + 
    geom_hline(yintercept=0, linetype="dashed") +
    geom_vline(xintercept=10, linetype="dashed") +
    geom_vline(xintercept=15, linetype="dashed") +
    labs(
      title = paste(var[i]),
      subtitle = paste("PH test p =", signif(ph_p,3)),
      x = "Follow-up time",
      y = "Scaled Schoenfeld residuals") +
    theme_classic()+
    theme(axis.title = element_text(size=9,color="black",family="calibri"),
          plot.title = element_text(hjust=0.5,size=9,color="black",family="calibri"),
          plot.subtitle = element_text(hjust=0.5,size=9,color="black",family="calibri"),
          panel.grid=element_blank(),
          legend.position = "bottom",
          axis.text = element_text(size=9,color="black",family="calibri"))
}

p <- grid.arrange(grobs = plots, nrow = 4)

ggsave("PH_test.png", p, width = 12, height = 8, dpi = 600)

#--------------------------------------------------------------------------------------------
#
#           chunk2: association with T2D
#
#--------------------------------------------------------------------------------------------
data_use <- t2d_sample

var <- c("amed_av","ahei_av","dash_av","pdi_av","hpdi_av","updi_av","edip_av","edih_av",
         "amed","ahei","dash","pdi","hpdi","updi","edip","edih")

data_use[var] <- apply(data_use[var],2,inormal)

d11 <- singlefun1("amed_av","diabetes",data_use)        
d21 <- singlefun1("ahei_av","diabetes",data_use)     
d31 <- singlefun1("dash_av","diabetes",data_use)         
d41 <- singlefun1("pdi_av","diabetes",data_use)          
d51 <- singlefun1("hpdi_av","diabetes",data_use)         
d61 <- singlefun1("updi_av","diabetes",data_use)         
d71 <- singlefun1("edip_av","diabetes",data_use)          
d81 <- singlefun1("edih_av","diabetes",data_use)        

d12 <- singlefun1("amed","diabetes",data_use)           
d22 <- singlefun1("ahei","diabetes",data_use)          
d32 <- singlefun1("dash","diabetes",data_use)          
d42 <- singlefun1("pdi","diabetes",data_use)           
d52 <- singlefun1("hpdi","diabetes",data_use)            
d62 <- singlefun1("updi","diabetes",data_use)           
d72 <- singlefun1("edip","diabetes",data_use)            
d82 <- singlefun1("edih","diabetes",data_use) 

res1 <- rbind(cbind(d11,d12),cbind(d21,d22),cbind(d31,d32),cbind(d41,d42),
              cbind(d51,d52),cbind(d61,d62),cbind(d71,d72),cbind(d81,d82)) %>% as.data.frame()

res1$Trait <- rep(c("AMED","AHEI","DASH","PDI","hPDI","uPDI","EDIP","EDIH"),each=3)
res1$Type <- "Per SD"
res1$Model <- rep(c("Age","MV","MV+BMI"),times=8)

data_use$amed_av_Q4 = cat("amed_av","quartile",data_use)
data_use$ahei_av_Q4 = cat("ahei_av","quartile",data_use)
data_use$dash_av_Q4 = cat("dash_av","quartile",data_use)
data_use$pdi_av_Q4 = cat("pdi_av","quartile",data_use)
data_use$hpdi_av_Q4 = cat("hpdi_av","quartile",data_use)
data_use$updi_av_Q4 = cat("updi_av","quartile",data_use)
data_use$edip_av_Q4 = cat("edip_av","quartile",data_use)
data_use$edih_av_Q4 = cat("edih_av","quartile",data_use)

data_use$amed2_Q4 = cat("amed","quartile",data_use)
data_use$ahei2_Q4 = cat("ahei","quartile",data_use)
data_use$dash2_Q4 = cat("dash","quartile",data_use)
data_use$pdi2_Q4 = cat("pdi","quartile",data_use)
data_use$hpdi2_Q4 = cat("hpdi","quartile",data_use)
data_use$updi2_Q4 = cat("updi","quartile",data_use)
data_use$edip2_Q4 = cat("edip","quartile",data_use)
data_use$edih2_Q4 = cat("edih","quartile",data_use)

d11 <- singlefun2("amed_av_Q4","diabetes",data_use)        
d21 <- singlefun2("ahei_av_Q4","diabetes",data_use)     
d31 <- singlefun2("dash_av_Q4","diabetes",data_use)         
d41 <- singlefun2("pdi_av_Q4","diabetes",data_use)          
d51 <- singlefun2("hpdi_av_Q4","diabetes",data_use)         
d61 <- singlefun2("updi_av_Q4","diabetes",data_use)         
d71 <- singlefun2("edip_av_Q4","diabetes",data_use)          
d81 <- singlefun2("edih_av_Q4","diabetes",data_use)        

d12 <- singlefun2("amed2_Q4","diabetes",data_use)           
d22 <- singlefun2("ahei2_Q4","diabetes",data_use)          
d32 <- singlefun2("dash2_Q4","diabetes",data_use)          
d42 <- singlefun2("pdi2_Q4","diabetes",data_use)           
d52 <- singlefun2("hpdi2_Q4","diabetes",data_use)            
d62 <- singlefun2("updi2_Q4","diabetes",data_use)           
d72 <- singlefun2("edip2_Q4","diabetes",data_use)            
d82 <- singlefun2("edih2_Q4","diabetes",data_use) 

d13 <- singlefun3("amed_av_Q4","diabetes",data_use)           
d23 <- singlefun3("ahei_av_Q4","diabetes",data_use)          
d33 <- singlefun3("dash_av_Q4","diabetes",data_use)          
d43 <- singlefun3("pdi_av_Q4","diabetes",data_use)           
d53 <- singlefun3("hpdi_av_Q4","diabetes",data_use)            
d63 <- singlefun3("updi_av_Q4","diabetes",data_use)           
d73 <- singlefun3("edip_av_Q4","diabetes",data_use)            
d83 <- singlefun3("edih_av_Q4","diabetes",data_use) 

d14 <- singlefun3("amed2_Q4","diabetes",data_use)           
d24 <- singlefun3("ahei2_Q4","diabetes",data_use)          
d34 <- singlefun3("dash2_Q4","diabetes",data_use)          
d44 <- singlefun3("pdi2_Q4","diabetes",data_use)           
d54 <- singlefun3("hpdi2_Q4","diabetes",data_use)            
d64 <- singlefun3("updi2_Q4","diabetes",data_use)           
d74 <- singlefun3("edip2_Q4","diabetes",data_use)            
d84 <- singlefun3("edih2_Q4","diabetes",data_use) 

res2 <- rbind(cbind(d11,d12),cbind(d21,d22),cbind(d31,d32),cbind(d41,d42),
              cbind(d51,d52),cbind(d61,d62),cbind(d71,d72),cbind(d81,d82)) %>% as.data.frame()

res3 <- rbind(cbind(d13,d14),cbind(d23,d24),cbind(d33,d34),cbind(d43,d44),
              cbind(d53,d54),cbind(d63,d64),cbind(d73,d74),cbind(d83,d84)) %>% as.data.frame()

res2$Trait <- rep(c("AMED","AHEI","DASH","PDI","hPDI","uPDI","EDIP","EDIH"),each=9)
res2$Type <- rep(c("Q2","Q3","Q4","Q2","Q3","Q4","Q2","Q3","Q4"),times=8)
res2$Model <- rep(c("Age","Age","Age","MV","MV","MV","MV+BMI","MV+BMI","MV+BMI"),times=8)

res3$Trait <- rep(c("AMED","AHEI","DASH","PDI","hPDI","uPDI","EDIP","EDIH"),each=3)
res3$Type <- "Per SD"
res3$Model <- rep(c("Age","MV","MV+BMI"),times=8)

#--------------------------------------------------------------------------------------------
#
#           chunk3: mediation analysis
#
#--------------------------------------------------------------------------------------------
#define dataset for use
data <- t2d_sample

#define exposure for use. here we use AMED as an example
exp = "amed_av"

#define covariate
adj_covar <- c("ageyr","aspirinuse","fhxdb","antihluse","smoking","act","phxchol","phxhbp","energy","BMIcont")

#define mediator
amed0 <- na.omit(sig_list[,1:2]); amed <- data.frame(amed0[,-1]); rownames(amed) <- amed0[,1]
ahei0 <- na.omit(sig_list[,4:5]); ahei <- data.frame(ahei0[,-1]); rownames(ahei) <- ahei0[,1]
dash0 <- na.omit(sig_list[,7:8]); dash <- data.frame(dash0[,-1]); rownames(dash) <- dash0[,1]
opdi0 <- na.omit(sig_list[,10:11]); opdi <- data.frame(opdi0[,-1]); rownames(opdi) <- opdi0[,1]
hpdi0 <- na.omit(sig_list[,13:14]); hpdi <- data.frame(hpdi0[,-1]); rownames(hpdi) <- hpdi0[,1]
updi0 <- na.omit(sig_list[,16:17]); updi <- data.frame(updi0[,-1]); rownames(updi) <- updi0[,1]
edip0 <- na.omit(sig_list[,19:20]); edip <- data.frame(edip0[,-1]); rownames(edip) <- edip0[,1]
edih0 <- na.omit(sig_list[,22:23]); edih <- data.frame(edih0[,-1]); rownames(edih) <- edih0[,1]

names(amed0)[1] <- names(ahei0)[1] <- names(dash0)[1] <- names(opdi0)[1] <- names(hpdi0)[1] <- names(updi0)[1] <- names(edip0)[1] <- names(edih0)[1] <- "HMDB"

metab_names <- amed0[-1,]$HMDB #metabolites in the AMED signature

#inverse-normal transformation
data[,c(exp,metab_names)] <- apply(data[,c(exp,metab_names)],2,inormal)

#define the controlled values of mediators
metab_dataset <- data %>% dplyr::select(metab_names)

controlled_med<- metab_dataset%>%
  dplyr::select(metab_names)%>%reshape2::melt()%>%
  group_by(variable)%>%
  summarize(median=median(value,na.rm = T))

controlled_med<-as.vector(controlled_med$median)

#run model
nb <- 200

singlemed_test_result<-list()  

for(i in 1:length(metab_names)){
  
  dati =data[,c(exp_x,metab_names[i],"ptime_diabetes","diabetes",adj_covar)] #note: mediator,exposure,outcome,event can not be NA
  dati = dati[complete.cases(dati),]
  
  adj_covar_use = c("ageyr","act","energy","BMIcont")
  
  #add fast
  if(length(unique(dati$fast))>1) {
    adj_covar_use = c(adj_covar_use,"fast")
  } 
  
  #add aspirinuse
  if(length(unique(dati$aspirinuse))>1) {
    adj_covar_use = c(adj_covar_use,"aspirinuse")
  } 
  
  
  #add antihluse
  if(length(unique(dati$antihluse))>1) {
    adj_covar_use = c(adj_covar_use,"antihluse")
  }
  
  #add smoking
  if(length(unique(dati$smoking))>1) {
    adj_covar_use = c(adj_covar_use,"smoking")
  }
  
  
  #add fhxdb
  if(length(unique(dati$fhxdb))>1) {
    adj_covar_use = c(adj_covar_use,"fhxdb")
  }
  
  #add phxchol
  if(length(unique(dati$phxchol))>1) {
    adj_covar_use = c(adj_covar_use,"phxchol")
  }
  
  #add phxhbp
  if(length(unique(dati$phxhbp))>1) {
    adj_covar_use = c(adj_covar_use,"phxhbp")
  }
  
  med_list <- as.vector(metab_names[i]) 
  med_controll_value_list <- as.list(controlled_med[1])      #this is a list with single value
  mreg_list <- as.list(rep("linear",length(metab_names[i]))) #this for models in mediator analysis 
  
  model_cmest<- cmest(data = dati, 
                      model = "rb", 
                      outcome = "ptime_diabetes",          
                      event = "diabetes",                   
                      exposure = exp,                   
                      mediator = med_list,          
                      basec = c(adj_covar_use), 
                      EMint = FALSE,
                      mreg =mreg_list, 
                      yreg = "coxph",
                      # define the exposure change values. Results could be different if astar=1, a=2, or astar=-1,a=0
                      astar = 0, a = 1, 
                      mval = med_controll_value_list , # controlled values for mediators when calculating the controlled effects
                      estimation = "imputation",
                      inference = "bootstrap",
                      nboot=nb) 
  
  #print(paste0(a,"-",i)) #to monitoring the process of cmest model
  
  save_singlemed_result<-as.data.frame(as.matrix(summary(model_cmest)$summarydf))
  save_singlemed_result<-save_singlemed_result%>%rownames_to_column()
  
  save_singlemed_result$mediators<-paste0(metab_names[i],collapse = "+")
  singlemed_test_result[[i]]<-save_singlemed_result
}

singlemed_test_result_save <- do.call(rbind.data.frame,singlemed_test_result)

singlemed_test_result_save <- singlemed_test_result_save %>% group_by(rowname) %>% mutate(fdr=p.adjust(`P.val`,"BH"))

singlemed_test_result_save$exposure<-exp
singlemed_test_result_save$type<-"Single_mediator"

rm(list = ls())

#--------------------------------------------------------------------------------------------
#
#           chunk4: plotting
#
#--------------------------------------------------------------------------------------------
#plot quartile results
t2d_demo$Diet <- factor(t2d_demo$Diet, levels = c("AMED","AHEI","DASH","PDI","hPDI","uPDI","EDIP","EDIH"))

p1 <- ggplot(t2d_demo, aes(x =exp(Est_Meta), y = Quartile, color=Type, fill=Type)) +
  geom_point(size = 4, shape=16,position = position_dodge(width = 8), width=1.2) +
  geom_errorbarh(aes(xmin = exp(Est_Meta-1.96*SE_Meta), xmax = exp(Est_Meta+1.96*SE_Meta)), position = position_dodge(width = 8),width=1.2,height = 0) +
  scale_colour_manual(values = c("#006EBE","#FA5555")) +
  #scale_x_continuous(limits=c(0.5,1.02), breaks=seq(0.5,1.02,0.1))+
  scale_y_discrete(expand = expansion(mult = c(0, 0.01)))+
  geom_vline(aes(xintercept = 1.0),colour="black",size=0.4) +
  labs(title = "",x = "HR (95% CI)",y = "") +
  facet_wrap(~Diet,scales="free",nrow=1) +
  theme_classic()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text=element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size=11,hjust=0.5,color="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_text(size=11,color="black"),
        axis.title = element_text(size=11,color="black"),
        panel.spacing.x = unit(0.5, "cm"),
        panel.spacing.y = unit(3, "cm"),
        legend.title = element_blank()) +
  coord_flip()

ggsave("Figure4_1.png", plot = p1, width = 9, height = 2, dpi = 600)

#plot mediation results
coef <- med_demo[,c(1:1,3:4)] %>%
  pivot_wider(
    names_from = name,   # column to spread into new columns
    values_from = PM     # column containing values
  ) %>% as.data.frame()

rownames(coef) <- coef$Diet
coef <- coef[,-1]

col_fun = colorRamp2(
  breaks = seq(0, 20, length.out = 8),  # Adjust range based on your data
  colors = brewer.pal(8, "OrRd")
)

coef <- coef[c(3:3,6:6,2:2,4:4,1:1,5:5,8:8,7:7),]

p2 <- Heatmap(coef, 
              col=col_fun,
              rect_gp = gpar(col = "grey70", lwd = 1),
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8,fontfamily = "Calibri"),  
              column_names_gp = gpar(fontsize = 7,fontfamily = "Calibri"),
              row_gap = unit(2, "mm"),
              show_row_dend = FALSE,
              row_order = order(as.numeric(gsub("row", "", rownames(coef)))),
              row_title = NULL,
              column_names_side = "bottom", 
              column_names_rot = 60,
              width  = unit(18, "cm"),
              height = unit(5, "cm"),
              na_col = "white",
              show_column_dend = FALSE,
              #bottom_annotation = col_ha,
              column_order = order(as.numeric(gsub("column", "", colnames(coef)))),
              column_title = NULL,
              heatmap_legend_param = list(title = "Mediated proportion",labels_gp = gpar(fontfamily = "calibri", fontsize = 7),title_gp  = gpar(fontfamily = "calibri", fontface = "bold", fontsize = 7)),
              cell_fun = function(j, i, x, y, width, height, fill) {
                v <- coef[i, j]
                if (!is.na(v)) {
                  grid.text(print(v), x, y,gp = gpar(fontsize = 7, family="calibri",color="grey30"))
                }
              }
)

png("Figure4_2.png",width = 2600, height = 1200, res = 300)
p2
dev.off()

rm(list = ls())
