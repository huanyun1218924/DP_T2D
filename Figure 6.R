library(data.table)
library(dplyr)
library(readxl)
library(VennDiagram)
library(venn)
library(Maaslin2)
library(lme4)
library(lmerTest)
library(CMAverse)
library(RColorBrewer)
library(ComplexHeatmap)


#--------------------------------------------------------------------------------------------
#
#           chunk1: identify signature-related microbial species
#
#--------------------------------------------------------------------------------------------
#read sample data
load("micro_sample.RData")

data_use <- micro_sample

#define species for use
Microbiome<-anno_sample$microNamePhenotypes = setdiff(names(data_use),Microbiome)

#define exposures for use
var_use <- c("amed","ahei","dash","pdi","hpdi","updi","edip","edih")

met_name <- c("AMED","AHEI","DASH","PDI","hPDI","uPDI","EDIP","EDIH")

#define covariates for uese
AdjVars = c("age_fec","bmi_bld","totMETs_paq","smoke_bld","probio_2m_fec","antibio_12m_fec","colsc_2m_fec","acid_2m_fec","stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6")

#perform inverse-normal transformation for signatures
inormal <- function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

data_use[var_use] <- apply(data_use[var_use],2,inormal) 

#run micro-was for each signature using maaslin2
for(i in 1:length(var_use)) {
  
  vari = var_use[i]
  namei= met_name[i]
  
  Maaslin2(input_data=data_use[,Microbiome],
           input_metadata=data_use[,Phenotypes],
           output=paste("~dir/",namei,sep=''), 
           min_abundance = 0.0001,
           min_prevalence = 0.1, 
           normalization = "None",  #Should set to NONE because species data were already normalized by TSS
           transform = "AST", 
           analysis_method = "LM",
           fixed_effects = c(vari,AdjVars), 
           correction = "BH",
           standardize = TRUE, 
           cores = 8,
           plot_heatmap=FALSE,
           plot_scatter=FALSE)
  
  print(i)
  
}

#--------------------------------------------------------------------------------------------
#
#           chunk2: association between individual metabolites and microbial species
#
#--------------------------------------------------------------------------------------------
#data transformation
data_use[sp] <- apply(data_use[sp],2,function(x)asin(sqrt(x)))  #species associated with signatures
data_use[var] <- apply(data_use[var],2,inormal)                 #metabolites included in the signatures
data_use[exp] <- apply(data_use[exp],2,inormal)                 #dietary pattern scores

#association between diet (exposure) and GMB (mediator)
res = data.frame(Diet=NA,microName=NA,Est=NA,se=NA,P=NA)

x = 0

for (i in exp){
  for (j in sp){
    
    data_use[,"Score1"] <- data_use[,i] #diet
    data_use[,"Score2"] <- data_use[,j] #GMB
    
    tryCatch({
      fit = glm(Score2~Score1+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data_use)
      
      x = x+1
      res[x,"Diet"] = i
      res[x,"microName"] = j
      res[x,"Est"] = coef(summary(fit))["Score1",1:1]
      res[x,"se"] = coef(summary(fit))["Score1",2:2]
      res[x,"P"] = coef(summary(fit))["Score1",4:4]
    }
    , error = function(e) {
      
      # Handle error: skip this variable and move on
      print(paste("Error running model for", i, ":", e$message))
      
    })
  }
}

diet_gmb <- left_join(res,anno_sample[,c("microName","species")],by="microName")

#association between metabolites (outcome) and GMB (mediator)
res = data.frame(microName=NA,HMDB=NA,Est=NA,se=NA,P=NA)

x = 0

for (i in var){
  for (j in sp){
    
    data_use[,"Met"] <- data_use[,i]
    data_use[,"GMB"] <- data_use[,j]
    
    tryCatch({
      fit = glm(GMB~Met+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data_use)
      
      x = x+1
      res[x,"microName"] = j
      res[x,"HMDB"] = i
      res[x,"Est"] = coef(summary(fit))["Met",1:1]
      res[x,"se"] = coef(summary(fit))["Met",2:2]
      res[x,"P"] = coef(summary(fit))["Met",4:4]
    }
    , error = function(e) {
      
      # Handle error: skip this variable and move on
      print(paste("Error running model for", i, ":", e$message))
      
    })
  }
}

gmb_met <- left_join(res,anno_sample[,c("microName","species")],by="microName") %>% left_join(anno[,c("HMDB","Name")],by="HMDB")

#association between diet (exposure) and metabolite (outcome)
res = data.frame(Diet=NA,HMDB=NA,Est=NA,se=NA,P=NA)

x = 0

for (i in exp){
  for (j in var){
    
    data_use[,"Score1"] <- data_use[,i] #diet
    data_use[,"Score2"] <- data_use[,j] #metabolite
    
    tryCatch({
      fit = glm(Score2~Score1+age_fec+bmi_bld+totMETs_paq+smoke_bld+probio_2m_fec+antibio_12m_fec+colsc_2m_fec+acid_2m_fec+stooltype_fec.1+stooltype_fec.2+stooltype_fec.3+stooltype_fec.4+stooltype_fec.5+stooltype_fec.6, data=data_use)
      
      x = x+1
      res[x,"Diet"] = i
      res[x,"HMDB"] = j
      res[x,"Est"] = coef(summary(fit))["Score1",1:1]
      res[x,"se"] = coef(summary(fit))["Score1",2:2]
      res[x,"P"] = coef(summary(fit))["Score1",4:4]
    }
    , error = function(e) {
      
      # Handle error: skip this variable and move on
      print(paste("Error running model for", i, ":", e$message))
      
    })
  }
}

diet_met <- left_join(res,anno[,c("HMDB","Name")],by="HMDB")

#--------------------------------------------------------------------------------------------
#
#           chunk3: mediation analysis
#
#--------------------------------------------------------------------------------------------
#define dataset
data_use <- micro_sample

#here we take AMED as an example
covar <- c("ageyr","race","bmi","act","alco","smoke","calor","stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6","stooltype_fec.7","stooltype_fec.8","antibio_12m_fec","colsc_2m_fec","probio_2m_fec","acid_2m_fec")

#define exposure
exp1 <- "amed_av"

#define mediator
med1 <- sp 

#define outcome
out1 <- "amed"

#data transformation
data_use[exp1] <- apply(data_use[exp1],2,inormal)
data_use[med1] <- apply(data_use[med1],2,function(x)asin(sqrt(x)))
data_use[out1] <- apply(data_use[out1],2,inormal)

metab_names <- med1

metab_dataset <- data_use %>% dplyr::select(metab_names)

controlled_med<- metab_dataset%>%
  dplyr::select(metab_names)%>%reshape2::melt()%>%
  group_by(variable)%>%
  summarize(median=median(value,na.rm = T))

controlled_med<-as.vector(controlled_med$median)

#run mediation analysis
singlemed_test_result<-list()  

exp_x <- exp1

for(i in 1:length(metab_names)){
  
  dati =data_use[,c(exp_x,metab_names[i],out1,covar)] #note: mediator,exposure,outcome,event can not be NA
  dati = dati[complete.cases(dati),]
  dati[,exp_x] <- as.numeric(dati[,exp_x])
  dati[,metab_names[i]] <- as.numeric(dati[,metab_names[i]])
  
  med_list <- as.vector(metab_names[i]) 
  med_controll_value_list <- as.list(controlled_med[1])      #this is a list with single value
  mreg_list <- as.list(rep("linear",length(metab_names[i]))) #this for models in mediator analysis 
  
  model_cmest<- cmest(data = dati,model = "rb",outcome = out1[1],exposure = exp_x,mediator = med_list,basec = c(covar),EMint = FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=100) 
  
  #print(paste0(a,"-",i)) #to monitoring the process of cmest model
  
  save_singlemed_result<-as.data.frame(as.matrix(summary(model_cmest)$summarydf))
  save_singlemed_result<-save_singlemed_result%>%rownames_to_column()
  
  save_singlemed_result$mediators<-paste0(metab_names[i],collapse = "+")
  singlemed_test_result[[i]]<-save_singlemed_result
}

singlemed_test_result_save <- do.call(rbind.data.frame,singlemed_test_result)

singlemed_test_result_save <- singlemed_test_result_save %>% group_by(rowname) %>% mutate(fdr=p.adjust(`P.val`,"BH"))

singlemed_test_result_save$exposure<-exp_x
singlemed_test_result_save$outcome<-out1[1]

#--------------------------------------------------------------------------------------------
#
#           chunk4: plotting
#
#--------------------------------------------------------------------------------------------
col_fun = colorRamp2(
  breaks = seq(0.01, -0.01, length.out = 9),  # Adjust range based on your data
  colors = brewer.pal(9, "RdBu")
)

ht <- Heatmap(coef, 
        col=col_fun,
        rect_gp = gpar(col = "grey50", lwd = 1),
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8.5,fontfamily = "Calibri"),  
        column_names_gp = gpar(fontsize = 8.5,fontfamily = "Calibri"),
        row_gap = unit(2, "mm"),
        show_row_dend = FALSE,
        row_order = order(as.numeric(gsub("row", "", rownames(coef)))),
        row_title = NULL,
        column_names_side = "bottom", 
        column_names_rot = 50,
        na_col = "white",
        width  = unit(11.5, "cm"),
        height = unit(15.5, "cm"),
        show_column_dend = FALSE,
        column_order = order(as.numeric(gsub("column", "", colnames(coef)))),
        column_title = NULL,
        heatmap_legend_param = list(title = "Coefficient",labels_gp = gpar(fontfamily = "calibri", fontsize = 9),title_gp  = gpar(fontfamily = "calibri", fontface = "bold", fontsize = 10))
)

png("Diet_GMB_Met.png",width = 2500, height = 3000, res = 300)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()
