#load packages
library(data.table)
library(dplyr)
library(plyr)
library(glmnet)
library(ggplot2)
library(glmnet)
library(cvTools)

#--------------------------------------------------------------------------------------------
#
#           chunk1: develop metabolomic signatures for dietary pattern scores using elastic net regression
#           note: here we use AMED as an example
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

rm(list = ls())

#--------------------------------------------------------------------------------------------
#
#           chunk2: metabolomic signature calculation in cohorts 
#           Note: here we calculate metabolomic signature in the test samples, with available metabolites used for calculation 
#
#--------------------------------------------------------------------------------------------
#read sample metabolome data
load("signature_sample.RData")  #it contains train sample data, test sample data, annotation file, and signature information

data_use <- test_sample

#read signature information: 
amed0 <- na.omit(sig_list[,1:2]); amed <- data.frame(amed0[,-1]); rownames(amed) <- amed0[,1]
ahei0 <- na.omit(sig_list[,3:4]); ahei <- data.frame(ahei0[,-1]); rownames(ahei) <- ahei0[,1]
dash0 <- na.omit(sig_list[,5:6]); dash <- data.frame(dash0[,-1]); rownames(dash) <- dash0[,1]
opdi0 <- na.omit(sig_list[,7:8]); opdi <- data.frame(opdi0[,-1]); rownames(opdi) <- opdi0[,1]
hpdi0 <- na.omit(sig_list[,9:10]); hpdi <- data.frame(hpdi0[,-1]); rownames(hpdi) <- hpdi0[,1]
updi0 <- na.omit(sig_list[,11:12]); updi <- data.frame(updi0[,-1]); rownames(updi) <- updi0[,1]
edip0 <- na.omit(sig_list[,13:14]); edip <- data.frame(edip0[,-1]); rownames(edip) <- edip0[,1]
edih0 <- na.omit(sig_list[,15:16]); edih <- data.frame(edih0[,-1]); rownames(edih) <- edih0[,1]

#make sure the coefficients from the elastic net model are numeric type
amed$Coefficient <- as.numeric(amed$Coefficient)
ahei$Coefficient <- as.numeric(ahei$Coefficient)
dash$Coefficient <- as.numeric(dash$Coefficient)
opdi$Coefficient <- as.numeric(opdi$Coefficient)
hpdi$Coefficient <- as.numeric(hpdi$Coefficient)
updi$Coefficient <- as.numeric(updi$Coefficient)
edip$Coefficient <- as.numeric(edip$Coefficient)
edih$Coefficient <- as.numeric(edih$Coefficient)

#extract metabolite name information
var <- names(data_use)[10:ncol(test_sample)]  

#standarize metabolites
data_use_std <- data_use
data_use_std[var] <- apply(data_use_std[var],2,scale)

#calculate metabolomic signature (formula:weighted sum of selected standarized metabolite plus intercept)
data_use_std$amed = apply(mapply(`*`, data_use_std[,which(colnames(data_use_std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(data_use_std)),]))), 1, sum)+as.numeric(amed[1,1])
data_use_std$ahei = apply(mapply(`*`, data_use_std[,which(colnames(data_use_std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(data_use_std)),]))), 1, sum)+as.numeric(ahei[1,1])
data_use_std$dash = apply(mapply(`*`, data_use_std[,which(colnames(data_use_std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(data_use_std)),]))), 1, sum)+as.numeric(dash[1,1])
data_use_std$pdi  = apply(mapply(`*`, data_use_std[,which(colnames(data_use_std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(data_use_std)),]))), 1, sum)+as.numeric(opdi[1,1])
data_use_std$hpdi = apply(mapply(`*`, data_use_std[,which(colnames(data_use_std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(data_use_std)),]))), 1, sum)+as.numeric(hpdi[1,1])
data_use_std$updi = apply(mapply(`*`, data_use_std[,which(colnames(data_use_std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(data_use_std)),]))), 1, sum)+as.numeric(updi[1,1])
data_use_std$edip = apply(mapply(`*`, data_use_std[,which(colnames(data_use_std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(data_use_std)),]))), 1, sum)+as.numeric(edip[1,1])
data_use_std$edih = apply(mapply(`*`, data_use_std[,which(colnames(data_use_std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(data_use_std)),]))), 1, sum)+as.numeric(edih[1,1])

#combine the metabolomic signature with raw dataset for subsequent analysis
data_use_sig <- merge(data_use,data_use_std[,c("id","amed",ahei","dash","pdi","hpdi","updi","edip","edih")],by="id")

#--------------------------------------------------------------------------------------------
#
#           chunk3: estimate pearson correlation
#           Note: here we used the derived metabolomic signature in the test sample as an example
#
#--------------------------------------------------------------------------------------------
#define dietary pattern score and metabolomic signature variables
dp <- c("amed_score","ahei_score","dash_score","opdi_score","hpdi_score","updi_score","edip_score","edih_score")       #dietary pattern score
ms <- c("amed","ahei","dash","opdi","hpdi","updi","edip","edih")                                                       #metabolomic signature

#inverse-normal transformation for dietary pattern scores and metabolomic signatures
inormal <- function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

data_use_sig[,c(dp,ms)] <- apply(data_use_sig[,c(dp,ms)],2,inormal)

#calculate correlation between dietary scores by ffq and metabolic signature
data <- data_use_sig
pairs <- data.frame(x = dp,y = ms)

res <- pairs %>%
  mutate(result = map2(x, y,~ cor.test(data[[.x]],data[[.y]],method = "pearson")),r = map_dbl(result, ~ unname(.x$estimate)),LCI = map_dbl(result, ~ .x$conf.int[1]),UCI = map_dbl(result, ~ .x$conf.int[2]),p = map_dbl(result, ~ .x$p.value)) %>%
  select(-result) %>%
  mutate(result = sprintf("r = %.2f (95%% CI: %.2f–%.2f), P = %.2e",r, LCI, UCI, p))

res5$Study <- "WHI"
res5$FDR <- p.adjust(res5$p,method="fdr")

rm(list = ls())

#--------------------------------------------------------------------------------------------
#
#            chunk4 - count the subclass information and plotting 
#
#--------------------------------------------------------------------------------------------
#read signature information
load("signature_sample.RData")

amed0 <- na.omit(sig_list[,1:2]); amed <- data.frame(amed0[,-1]); rownames(amed) <- amed0[,1]
ahei0 <- na.omit(sig_list[,4:5]); ahei <- data.frame(ahei0[,-1]); rownames(ahei) <- ahei0[,1]
dash0 <- na.omit(sig_list[,7:8]); dash <- data.frame(dash0[,-1]); rownames(dash) <- dash0[,1]
opdi0 <- na.omit(sig_list[,10:11]); opdi <- data.frame(opdi0[,-1]); rownames(opdi) <- opdi0[,1]
hpdi0 <- na.omit(sig_list[,13:14]); hpdi <- data.frame(hpdi0[,-1]); rownames(hpdi) <- hpdi0[,1]
updi0 <- na.omit(sig_list[,16:17]); updi <- data.frame(updi0[,-1]); rownames(updi) <- updi0[,1]
edip0 <- na.omit(sig_list[,19:20]); edip <- data.frame(edip0[,-1]); rownames(edip) <- edip0[,1]
edih0 <- na.omit(sig_list[,22:23]); edih <- data.frame(edih0[,-1]); rownames(edih) <- edih0[,1]

names(amed0)[1] <- names(ahei0)[1] <- names(dash0)[1] <- names(opdi0)[1] <- names(hpdi0)[1] <- names(updi0)[1] <- names(edip0)[1] <- names(edih0)[1] <- "HMDB"

#get the number for each subclass
amed1 <- subset(anno,HMDB %in% amed0$HMDB); dim(amed1) 
ahei1 <- subset(anno,HMDB %in% ahei0$HMDB); dim(ahei1) 
dash1 <- subset(anno,HMDB %in% dash0$HMDB); dim(dash1) 
opdi1 <- subset(anno,HMDB %in% opdi0$HMDB); dim(opdi1) 
hpdi1 <- subset(anno,HMDB %in% hpdi0$HMDB); dim(hpdi1) 
updi1 <- subset(anno,HMDB %in% updi0$HMDB); dim(updi1) 
edip1 <- subset(anno,HMDB %in% edip0$HMDB); dim(edip1) 
edih1 <- subset(anno,HMDB %in% edih0$HMDB); dim(edih1)

amed2 <- subset(anno,HMDB %in% amed0$HMDB); dim(amed2) 
ahei2 <- subset(anno,HMDB %in% ahei0$HMDB); dim(ahei2) 
dash2 <- subset(anno,HMDB %in% dash0$HMDB); dim(dash2) 
opdi2 <- subset(anno,HMDB %in% opdi0$HMDB); dim(opdi2) 
hpdi2 <- subset(anno,HMDB %in% hpdi0$HMDB); dim(hpdi2) 
updi2 <- subset(anno,HMDB %in% updi0$HMDB); dim(updi2) 
edip2 <- subset(anno,HMDB %in% edip0$HMDB); dim(edip2)
edih2 <- subset(anno,HMDB %in% edih0$HMDB); dim(edih2) 

met_num <- dplyr::bind_rows(table(amed2$Subclass),table(ahei2$Subclass),
                            table(dash2$Subclass),table(opdi2$Subclass),
                            table(hpdi2$Subclass),table(updi2$Subclass),
                            table(edip2$Subclass),table(edih2$Subclass)) %>% t() %>% as.data.frame()

names(met_num) <- c("AMED","AHEI","DASH","PDI","hPDI","uPDI","EDIP","EDIH")
met_num$name <- rownames(met_num)

melted_met_num <- reshape2::melt(met_num, id.vars = "name")
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

#sessionInfo()
#R version 4.3.3 (2024-02-29)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Rocky Linux 9.7 (Blue Onyx)

#Matrix products: default
#BLAS:   /app/R-4.3.3@i86-rhel9.0/lib64/R/lib/libRblas.so 
#LAPACK: FlexiBLAS OPENBLAS-OPENMP;  LAPACK version 3.9.0

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] cvTools_0.3.3       robustbase_0.99-7   lattice_0.22-7      ggplot2_4.0.2       glmnet_4.1-10      
#[6] Matrix_1.6-5        plyr_1.8.9          dplyr_1.2.0         data.table_1.18.2.1
