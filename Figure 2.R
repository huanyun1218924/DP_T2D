#load packages
library(plyr)
library(glmnet)
library(ggplot2)
library(glmnet)
library(cvTools)

#--------------------------------------------------------------------------------------------
#
#           chunk1: model training
#
#--------------------------------------------------------------------------------------------
#read data
load("model_training_sample.RData")

train_sample[1:5,1:5] #training set
test_sample[1:5,1:5]  #testing set

dim(train_sample) 
dim(test_sample) 

trainingset = as.data.frame(train_sample)
testingset = as.data.frame(test_sample)
met = as.character(names(train_sample)[10:295])

#repeat elastic net for 500 times with all the predictions
record = data.frame(repNo=NA,NoMetabs=NA,sol_cor=NA)
pdf("TestingPlots_AMED.pdf",width=8, height=24)

for (inrep in 1:500) { 
  
  # traning set on the testing sets
  
  Training_CV = cv.glmnet(as.matrix(trainingset[,met]), trainingset[,"amed1"], nfolds=10, alpha=0.5, family="gaussian")
  lambda_1se_10F = Training_CV$lambda.1se
  lambda_min_10F = Training_CV$lambda.min
  
  Training_M = glmnet(as.matrix(trainingset[,met]), trainingset[,"amed1"], family="gaussian", alpha=0.5)
  
  cmin = coef(Training_M, s=lambda_min_10F)
  metabsmin = data.frame(cmin[which(cmin[,1]!=0),])
  names(metabsmin) = "Test.min"
  dim(metabsmin)
  
  write.table(metabsmin,paste("AMED_inrep",inrep,".txt",sep=''))
  
  #### apply
  testingset[,dim(testingset)[2]+1] = as.numeric(predict(Training_M, as.matrix(testingset[,met]), type="response", s=lambda_min_10F)[,1])
  names(testingset)[dim(testingset)[2]] = paste("AMED_inrep",inrep,sep='')
  
  #### test statistics
  record[inrep,"repNo"] = paste("inrep",inrep,sep='')
  record[inrep,"NoMetabs"] = dim(metabsmin)[1]-1
  record[inrep,"sol_cor"] = cor(testingset[,"amed1"][which(testingset[,dim(testingset)[2]]>0 & testingset[,dim(testingset)[2]]<9)],testingset[,dim(testingset)[2]][which(testingset[,dim(testingset)[2]]>0 & testingset[,dim(testingset)[2]]<9)])

  #### Plot performance
  par(mfrow=c(1,2))
  
  plot(testingset[,"amed1"][which(testingset[,dim(testingset)[2]]>0 & testingset[,dim(testingset)[2]]<9)],testingset[,dim(testingset)[2]][which(testingset[,dim(testingset)[2]]>0 & testingset[,dim(testingset)[2]]<9)])
  abline(lm(predicted ~ selfreported,data=data.frame(selfreported=testingset[,"amed1"][which(testingset[,dim(testingset)[2]]>0 & testingset[,dim(testingset)[2]]<9)],predicted=testingset[,dim(testingset)[2]][which(testingset[,dim(testingset)[2]]>0 & testingset[,dim(testingset)[2]]<9)])),col="darkred",lwd=2)
  plot(as.factor(testingset[,"amed1"][which(testingset[,dim(testingset)[2]]>0 & testingset[,dim(testingset)[2]]<9)]),testingset[,dim(testingset)[2]][which(testingset[,dim(testingset)[2]]>0 & testingset[,dim(testingset)[2]]<9)])
  
}

dev.off()

#--------------------------------------------------------------------------------------------
#
#           chunk2: estimate pearson correlation
#
#--------------------------------------------------------------------------------------------
#calculate correlation between dietary scores by ffq and metabolic signature
dqs <- c("id","study","diabetes","amed_av","ahei_av","dash_av","pdi_av","hpdi_av","updi_av","edip_av","edih_av")
sig <- c("amed2","ahei2","dash2","pdi2","hpdi2","updi2","edip2","edih2")

dp <- c("amed_av","ahei_av","dash_av","pdi_av","hpdi_av","updi_av","edip_av","edih_av")
rs <- cor(all[dp], all[sig], method = "pearson", use = "pairwise") %>% as.data.frame()

#--------------------------------------------------------------------------------------------
#
#            chunk3 - plotting
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
