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
Taxon_use_new_inv <- Taxon_use_new

Phenotypes = setdiff(names(Taxon_use_new_inv),Microbiome)

var_use <- c("amed","ahei","dash","pdi","hpdi","updi","edip","edih")

met_name <- c("AMED","AHEI","DASH","PDI","hPDI","uPDI","EDIP","EDIH")

AdjVars = c("AGE_V2","GENDER_V2","CIGARETTE_USE_V2","ALCOHOL_USE_V2","GPAQ_TOTAL_MET","DIAB_FAMHIST","HIGH_TOTAL_CHOL2_V2","HYPERTENSION2_V2","BMI_V2")

Taxon_use_new_inv[var_use] <- apply(Taxon_use_new_inv[var_use],2,inormal) 

#of note, the value is between 0 and 100 and should be transformed into between 0 and 1
Taxon_use_new_inv[,anno_sol$microName] <- Taxon_use_new_inv[,anno_sol$microName]/100

for(i in 1:length(var_use)) {
  
  vari = var_use[i]
  namei= met_name[i]
  
  Maaslin2(input_data=Taxon_use_new_inv[,Microbiome],
           input_metadata=Taxon_use_new_inv[,Phenotypes],
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
load("Data_For_Med.RData")

data_use[sp] <- apply(data_use[sp],2,function(x)asin(sqrt(x)))
data_use[var2] <- apply(data_use[var2],2,inormal)
data_use[exp] <- apply(data_use[exp],2,inormal)

#association between diet (exposure) and GMB (mediator)
res = data.frame(Diet=NA,microName=NA,Est=NA,se=NA,P=NA)

x = 0

for (i in exp){
  for (j in sp){
    
    data_use[,"Score1"] <- data_use[,i] #diet
    data_use[,"Score2"] <- data_use[,j] #GMB
    
    tryCatch({
      fit = glm(Score2~Score1+AGE_V2+GENDER_V2+ALCOHOL_USE_V2+CIGARETTE_USE_V2+HIGH_TOTAL_CHOL2_V2+HYPERTENSION2_V2+DYSLIPIDEMIA_V2+BMI_V2, data=data_use)
      
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

diet_gmb_sol <- left_join(res,anno_sol[,c("microName","species")],by="microName")

#association between metabolites (outcome) and GMB (mediator)
res = data.frame(microName=NA,HMDB=NA,Est=NA,se=NA,P=NA)

x = 0

for (i in var2){
  for (j in sp){
    
    data_use[,"Met"] <- data_use[,i]
    data_use[,"GMB"] <- data_use[,j]
    
    tryCatch({
      fit = glm(GMB~Met+AGE_V2+GENDER_V2+ALCOHOL_USE_V2+CIGARETTE_USE_V2+HIGH_TOTAL_CHOL2_V2+HYPERTENSION2_V2+DYSLIPIDEMIA_V2+BMI_V2, data=data_use)
      
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

gmb_met_sol <- left_join(res,anno_sol[,c("microName","species")],by="microName") %>% left_join(met_t2d[,c("HMDB","name")],by="HMDB")

#association between diet (exposure) and metabolite (outcome)
res = data.frame(Diet=NA,HMDB=NA,Est=NA,se=NA,P=NA)

x = 0

for (i in exp){
  for (j in var2){
    
    data_use[,"Score1"] <- data_use[,i] #diet
    data_use[,"Score2"] <- data_use[,j] #metabolite
    
    tryCatch({
      fit = glm(Score2~Score1+AGE_V2+GENDER_V2+ALCOHOL_USE_V2+CIGARETTE_USE_V2+HIGH_TOTAL_CHOL2_V2+HYPERTENSION2_V2+DYSLIPIDEMIA_V2+BMI_V2, data=data_use)
      
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

diet_met_sol <- left_join(res,met_t2d[,c("HMDB","name")],by="HMDB")

#--------------------------------------------------------------------------------------------
#
#           chunk3: mediation analysis
#
#--------------------------------------------------------------------------------------------
covar <- c("AGE_V2","GENDER_V2","ALCOHOL_USE_V2","CIGARETTE_USE_V2","HIGH_TOTAL_CHOL2_V2","HYPERTENSION2_V2","DYSLIPIDEMIA_V2","BMI_V2")

exp1 <- "amed_score"
med1 <- unique(res_sol[which(res_sol$species %in% unique(res_mwas_sig[which(res_mwas_sig$Diet=="AMED"),]$species)),]$microName) 
out1 <- "amed"

all.int.use <- sol_use

all.int.use[exp1] <- apply(all.int.use[exp1],2,inormal)
all.int.use[med1] <- apply(all.int.use[med1],2,function(x)asin(sqrt(x)))
all.int.use[out1] <- apply(all.int.use[out1],2,inormal)

metab_names <- med1

metab_dataset <- all.int.use %>% dplyr::select(metab_names)

controlled_med<- metab_dataset%>%
  dplyr::select(metab_names)%>%data.table::melt()%>%
  group_by(variable)%>%
  summarize(median=median(value,na.rm = T))

controlled_med<-as.vector(controlled_med$median)

singlemed_test_result<-list()  

exp_x <- exp1

for(i in 1:length(metab_names)){
  
  dati =all.int.use[,c(exp_x,metab_names[i],out1,covar)] #note: mediator,exposure,outcome,event can not be NA
  dati = dati[complete.cases(dati),]
  dati[,exp_x] <- as.numeric(dati[,exp_x])
  dati[,metab_names[i]] <- as.numeric(dati[,metab_names[i]])
  
  med_list <- as.vector(metab_names[i]) 
  med_controll_value_list <- as.list(controlled_med[1])      #this is a list with single value
  mreg_list <- as.list(rep("linear",length(metab_names[i]))) #this for models in mediator analysis 
  
  model_cmest<- cmest(data = dati,model = "rb",outcome = out1[1],exposure = exp_x,mediator = med_list,basec = c(covar),EMint = FALSE,mreg=mreg_list <- as.list(rep("linear",1)),yreg="linear",astar=0,a=1,mval=as.list(controlled_med[1]),estimation="imputation",inference="bootstrap",nboot=100,yvar=1) 
  
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

dat1 <- singlemed_test_result_save

#--------------------------------------------------------------------------------------------
#
#           chunk4: plotting
#
#--------------------------------------------------------------------------------------------
col_fun = colorRamp2(
  breaks = seq(0.01, -0.01, length.out = 9),  # Adjust range based on your data
  colors = brewer.pal(9, "RdBu")
)

row_ha<-rowAnnotation(foo = substr(t2d_lab$sig, 1 , 1),
                      col = list(foo = c("a" = "#EE6A50", "b" = "#4F81BD", "c" = "grey"),gp = gpar(col = "black")),
                      annotation_name_side = "bottom",
                      annotation_name_rot=90,
                      annotation_legend_param = list(foo = list(title = "With T2D",labels = c("Positive","Negative","Non-significant"),labels_gp = gpar(fontfamily = "calibri", fontsize = 9),title_gp  = gpar(fontfamily = "calibri", fontface = "bold", fontsize = 10))),
                      annotation_label = "T2D",
                      annotation_name_gp = gpar(fontmaily="calibri",fontsize=9,fontface="bold"),
                      simple_anno_size = unit(0.2, "cm")) #trait classification

col_ha<-columnAnnotation(foo = substr(anno_lab$phylum, 1 , 1),
                      col = list(foo = c("a" = "#00b0f0", "b" = "#95a2ff", "f" = "#00b050", "p" = "#ffc000"),gp = gpar(col = "black")),
                      annotation_name_side = "right",annotation_name_rot=0,
                      annotation_legend_param = list(foo = list(title = "Phylum",labels = c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria"),labels_gp = gpar(fontfamily = "calibri", fontsize = 9),title_gp  = gpar(fontfamily = "calibri", fontface = "bold", fontsize = 10))),
                      annotation_label = "Phylum",
                      annotation_name_gp = gpar(fontmaily="calibri",fontsize=9,fontface="bold"),
                      simple_anno_size = unit(0.2, "cm")) #trait classification

col_ha2 <- columnAnnotation("Prevalence" = anno_barplot(pre_lab$Prevalance,border=TRUE,gp = gpar(fill = "#95a2ff",col="#95a2ff"),bar_width = 0.8),
                            annotation_name_side = "left",
                            annotation_name_rot=0,
                            annotation_label = NULL,
                            annotation_name_gp = gpar(fontmaily="calibri",fontsize=10,fontface="bold"),
                            col=list(c = ("grey90"),gp = gpar(col = "white"))) #overall h2

Heatmap(coef_f, 
        col=col_fun,
        rect_gp = gpar(col = "grey50", lwd = 1),
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8.5,fontfamily = "Calibri"),  
        column_names_gp = gpar(fontsize = 8.5,fontfamily = "Calibri"),
        row_split = c(rep('A',each=8),rep('B',each=8),rep('C',each=ncol(coef)-16)),
        row_gap = unit(2, "mm"),
        show_row_dend = FALSE,
        row_order = order(as.numeric(gsub("row", "", rownames(coef_f)))),
        row_title = NULL,
        column_names_side = "bottom", 
        column_names_rot = 50,
        na_col = "white",
        width  = unit(11.5, "cm"),
        height = unit(15.5, "cm"),
        show_column_dend = FALSE,
        bottom_annotation = col_ha,
        top_annotation = col_ha2,
        left_annotation = row_ha,
        column_order = order(as.numeric(gsub("column", "", colnames(coef_f)))),
        #column_split = LETTERS[1:dim(htmap_coef)[2]],
        column_title = NULL,
        heatmap_legend_param = list(title = "Coefficient",labels_gp = gpar(fontfamily = "calibri", fontsize = 9),title_gp  = gpar(fontfamily = "calibri", fontface = "bold", fontsize = 10)),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(pval[i, j] <= 0.05) {
            grid.text("**", x, y,gp = gpar(fontsize = 8, fontfamily="calibri",color="grey30"))
          } else if(pval2[i, j] <= 0.05) {
            grid.text("*", x, y,gp = gpar(fontsize = 8, fontfamily="calibri",color="grey30"))
          }
        }
)

png("Diet_GMB_Met.png",width = 2500, height = 3000, res = 300)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

ggplot(df, aes(x = group, y = item, size = PM, fill = group)) +
  geom_point(shape = 21, color = "black", alpha = 0.7) +
  geom_text(aes(label = paste(value,"%",sep='')), family="Calibri",size = 2, color = "black") +
  scale_fill_manual(values = c("AMED" = "#6A5ACD", "AHEI" = "#9370DB", "DASH" = "#BA55D3", "hPDI" = "#DDA0DD", "uPDI" = "#D8BFD8"))+
  scale_size_continuous(range = c(3, 10)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid=element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size=9,color="black",family="calibri"))