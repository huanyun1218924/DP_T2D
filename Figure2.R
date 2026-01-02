# --------------------------------------------------------------------------------
#Title: Dietary quality scores, metabolic signature, and cardiometabolic diseases
#Purpose: Correlation between dietary pattern and corresponding metabolic signatures (sensitivity analysis)
#Study: LVS, NHS/HPFS
#Path: /udd/nhhyu/DP_T2D/ProgramReview
#On: nantucket
#Programmer: Huan Yun (nhhyu)
#Date: 20230718
# --------------------------------------------------------------------------------

#load packages
library(data.table)
library(readxl)
library(dplyr)
library(plyr)
library(gmodels)
library(ggplot2)
library(cowplot)
require(tidyr)
require(vegan)
require(ape)

#set working dir
setwd("/udd/nhhyu/DP_T2D/FinalOutput")

#--------------------------------------------------------------------------------------------
#
#           Step1: estimate pearson correlation between diet and corresponding metabolic signature
#
#--------------------------------------------------------------------------------------------
#load data
load("Final_use_ms.RData")

#calculate correlation between dietary scores by ffq and metabolic signature
dqs <- c("id","study","diabetes","amed_av","ahei_av","dash_av","pdi_av","hpdi_av","updi_av","edip_av","edih_av")
sig <- c("amed2","ahei2","dash2","pdi2","hpdi2","updi2","edip2","edih2")

all <- rbind(nhs1_ms_use[,c(dqs,sig)],nhs2_ms_use[,c(dqs,sig)],hpfs_ms_use[,c(dqs,sig)])
all <- all[!duplicated(all$id),]

all$diabetes <- as.factor(all$diabetes)
all$study <- as.factor(all$study)

dp <- c("amed_av","ahei_av","dash_av","pdi_av","hpdi_av","updi_av","edip_av","edih_av")
rs <- cor(all[dp], all[sig], method = "pearson", use = "pairwise") %>% as.data.frame()

#save results
write.csv(rs, file="Figure2.csv")

#--------------------------------------------------------------------------------------------
#
#            step2 - plot the number of metabolites composing of each metabolic signature
#
#--------------------------------------------------------------------------------------------
#load signature profile
sig_list <- read_excel("Metabolic signature of DQSs_Feb 13.xlsx", sheet = "Cross-platform_SOL (n=122)") %>% as.data.frame()

#get signature for each dietary score
amed0 <- na.omit(sig_list[,1:2]); amed <- data.frame(amed0[,-1]); rownames(amed) <- amed0[,1]
ahei0 <- na.omit(sig_list[,4:5]); ahei <- data.frame(ahei0[,-1]); rownames(ahei) <- ahei0[,1]
dash0 <- na.omit(sig_list[,7:8]); dash <- data.frame(dash0[,-1]); rownames(dash) <- dash0[,1]
opdi0 <- na.omit(sig_list[,10:11]); opdi <- data.frame(opdi0[,-1]); rownames(opdi) <- opdi0[,1]
hpdi0 <- na.omit(sig_list[,13:14]); hpdi <- data.frame(hpdi0[,-1]); rownames(hpdi) <- hpdi0[,1]
updi0 <- na.omit(sig_list[,16:17]); updi <- data.frame(updi0[,-1]); rownames(updi) <- updi0[,1]
edip0 <- na.omit(sig_list[,19:20]); edip <- data.frame(edip0[,-1]); rownames(edip) <- edip0[,1]
edih0 <- na.omit(sig_list[,22:23]); edih <- data.frame(edih0[,-1]); rownames(edih) <- edih0[,1]

names(amed0)[1] <- names(ahei0)[1] <- names(dash0)[1] <- names(opdi0)[1] <- names(hpdi0)[1] <- names(updi0)[1] <- names(edip0)[1] <- names(edih0)[1] <- "HMDB"

#clean the annotation dataset
anno <- fread("Annotation_final.csv")
anno[218:227,3] <- "Other lipids"

#get the number for each subclass
amed1 <- subset(anno,HMDB %in% amed0$HMDB); dim(amed1) #77
ahei1 <- subset(anno,HMDB %in% ahei0$HMDB); dim(ahei1) #68
dash1 <- subset(anno,HMDB %in% dash0$HMDB); dim(dash1) #73
opdi1 <- subset(anno,HMDB %in% opdi0$HMDB); dim(opdi1) #46
hpdi1 <- subset(anno,HMDB %in% hpdi0$HMDB); dim(hpdi1) #58
updi1 <- subset(anno,HMDB %in% updi0$HMDB); dim(updi1) #51
edip1 <- subset(anno,HMDB %in% edip0$HMDB); dim(edip1) #97
edih1 <- subset(anno,HMDB %in% edih0$HMDB); dim(edih1) #59

amed2 <- subset(anno,HMDB %in% amed0$HMDB); dim(amed2) 
ahei2 <- subset(anno,HMDB %in% ahei0$HMDB); dim(ahei2) 
dash2 <- subset(anno,HMDB %in% dash0$HMDB); dim(dash2) 
opdi2 <- subset(anno,HMDB %in% opdi0$HMDB); dim(opdi2) 
hpdi2 <- subset(anno,HMDB %in% hpdi0$HMDB); dim(hpdi2) 
updi2 <- subset(anno,HMDB %in% updi0$HMDB); dim(updi2) 
edip2 <- subset(anno,HMDB %in% edip0$HMDB); dim(edip2)
edih2 <- subset(anno,HMDB %in% edih0$HMDB); dim(edih2) 

#deal with lipid issue
met_num <- dplyr::bind_rows(table(amed2$Subclass),table(ahei2$Subclass),
                            table(dash2$Subclass),table(opdi2$Subclass),
                            table(hpdi2$Subclass),table(updi2$Subclass),
                            table(edip2$Subclass),table(edih2$Subclass)) %>% t() %>% as.data.frame()

names(met_num) <- c("AMED","AHEI-2010","DASH","PDI","hPDI","uPDI","EDIP","EDIH")
met_num$name <- rownames(met_num)

melted_met_num <- melt(met_num, id.vars = "name")
melted_met_num <- drop_na(melted_met_num)
melted_met_num$value <- as.factor(melted_met_num$value)

melted_met_num$name <- factor(melted_met_num$name,levels = c("Carbohydrates","Amino acids","Acylcarnitines","Glycerophospholipids","Glycerolipids",
                                                             "Sphingolipids","Fatty acids","Other lipids","Cofactors and vitamins","Nucleotides","Xenobiotics"))
#plot stacked bar chart
pdf("Figure2.pdf", width = 7.5, height = 10)
ggplot(melted_met_num, aes(fill=name, y=value, x=variable)) + 
  geom_bar(position = "stack",stat='identity') +
  scale_fill_manual(values=c("#8cd2c8","#4F81BD","#beb9dc","#fa826e","#fab464","#b4dc64","#facde6","#95a2ff","#be82be","#B15928","#5F7530"),guide = guide_legend(reverse = TRUE))+
  geom_text(aes(label = value, y=value), position = position_stack(vjust = 0.5), size = 3) +
  #scale_y_continuous(breaks=seq(0,70,10),expand = c(0,0))+
  expand_limits(y=c(0, 70))+
  #scale_fill_discrete(limits = c("Amino acids and amines","Carbohydrates","Glycerolipids","Glycerophospholipids","Sphingolipids","Acylcarnitines","Fatty acids","Other lipids","Nucleotides","Cofactors and vitamins","Xenobiotics")) +
  theme_classic() +
  theme(legend.position = "right",
        axis.text.y = element_text(size=9,color="black"),
        axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size=9,color="black")) +
  labs(x= "Metabolic signature", y = "Number of metabolites", color = "Subclass") +
  guides(fill=guide_legend(title = "Subclass",ncol=1,byrow = FALSE)) 
dev.off()