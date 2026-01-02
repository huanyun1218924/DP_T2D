# --------------------------------------------------------------------------------
#Title: Dietary quality scores, metabolic signature, and cardiometabolic diseases
#Purpose: Explore shared and distinct mediation role of metabolic signature
#Study: NHS/HPFS
#Path: /udd/nhhyu/DP_T2D/ProgramReview
#On: nantucket
#Programmer: Huan Yun
#Date: 20231011
#Note: for healthy dietary pattern and plant-based dietary patterns, we will have four signatures
# --------------------------------------------------------------------------------

#load packages
library(data.table)
library(readxl)
library(writexl)
library(readxl)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(coxme)
library(CMAverse)

#set working dir
setwd("/udd/nhhyu/DP_T2D/FinalOutput")

#--------------------------------------------------------------------------------------------
#
#           Step1: get list of shared metabolic signature and distinct metabolic signature
#
#--------------------------------------------------------------------------------------------
#load signature profile
sig_list <- read_excel("Metabolic signature of DQSs_Feb 13.xlsx", sheet = "Cross-platform_SOL (n=122)") %>% as.data.frame()

#split metabolic signature for each dietary score
amed1 <- na.omit(sig_list[,1:2]); amed <- data.frame(amed1[,-1]); rownames(amed) <- amed1[,1]
ahei1 <- na.omit(sig_list[,4:5]); ahei <- data.frame(ahei1[,-1]); rownames(ahei) <- ahei1[,1]
dash1 <- na.omit(sig_list[,7:8]); dash <- data.frame(dash1[,-1]); rownames(dash) <- dash1[,1]
opdi1 <- na.omit(sig_list[,10:11]); opdi <- data.frame(opdi1[,-1]); rownames(opdi) <- opdi1[,1]
hpdi1 <- na.omit(sig_list[,13:14]); hpdi <- data.frame(hpdi1[,-1]); rownames(hpdi) <- hpdi1[,1]
updi1 <- na.omit(sig_list[,16:17]); updi <- data.frame(updi1[,-1]); rownames(updi) <- updi1[,1]
edip1 <- na.omit(sig_list[,19:20]); edip <- data.frame(edip1[,-1]); rownames(edip) <- edip1[,1]
edih1 <- na.omit(sig_list[,22:23]); edih <- data.frame(edih1[,-1]); rownames(edih) <- edih1[,1]

#get shared and distinct metabolites among AMED, AHEI, and DASH
hdp_1 <- intersect(amed1[-1,]$`AMED (n=48)`,ahei1[-1,]$`AHEI-2010 (n=61)`) %>% intersect(dash1[-1,]$`DASH (n=55)`) #32:amed&ahei&dash
hdp_2 <- intersect(amed1[-1,]$`AMED (n=48)`,ahei1[-1,]$`AHEI-2010 (n=61)`) %>% setdiff(hdp_1) #4:amed&ahei
hdp_3 <- intersect(amed1[-1,]$`AMED (n=48)`,dash1[-1,]$`DASH (n=55)`) %>% setdiff(hdp_1)      #6:amed&dash
hdp_4 <- intersect(ahei1[-1,]$`AHEI-2010 (n=61)`,dash1[-1,]$`DASH (n=55)`) %>% setdiff(hdp_1) #12:ahei&dash
hdp_5 <- setdiff(amed1[-1,]$`AMED (n=48)`,c(hdp_1,hdp_2,hdp_3))      #6:amed
hdp_6 <- setdiff(ahei1[-1,]$`AHEI-2010 (n=61)`,c(hdp_1,hdp_2,hdp_4)) #13:ahei
hdp_7 <- setdiff(dash1[-1,]$`DASH (n=55)`,c(hdp_1,hdp_3,hdp_4))      #5:dash

#get shared and distinct metabolites among three PDIs
pdp_1 <- intersect(opdi1[-1,]$`oPDI (n=36)`,hpdi1[-1,]$`hPDI (n=39)`) %>% intersect(updi1[-1,]$`uPDI (n=43)`) #10:pdi&hpdi&updi
pdp_2 <- intersect(opdi1[-1,]$`oPDI (n=36)`,hpdi1[-1,]$`hPDI (n=39)`) %>% setdiff(pdp_1) #7:pdi&hpdi
pdp_3 <- intersect(opdi1[-1,]$`oPDI (n=36)`,updi1[-1,]$`uPDI (n=43)`) %>% setdiff(pdp_1) #4:pdi&updi
pdp_4 <- intersect(hpdi1[-1,]$`hPDI (n=39)`,updi1[-1,]$`uPDI (n=43)`) %>% setdiff(pdp_1) #11:hpdi&updi
pdp_5 <- setdiff(opdi1[-1,]$`oPDI (n=36)`,c(pdp_1,pdp_2,pdp_3)) #15:pdi
pdp_6 <- setdiff(hpdi1[-1,]$`hPDI (n=39)`,c(pdp_1,pdp_2,pdp_4)) #10:hpdi
pdp_7 <- setdiff(updi1[-1,]$`uPDI (n=43)`,c(pdp_1,pdp_3,pdp_4)) #18:updi

#get shared and distinct metabolite between EDIP and EDIH
udp_1 <- intersect(edip1[-1,]$`EDIP (n=66)`,edih1[-1,]$`EDIH (n=37)`) #27:edip&edih
udp_2 <- setdiff(edip1[-1,]$`EDIP (n=66)`,udp_1) #39:edip
udp_3 <- setdiff(edih1[-1,]$`EDIH (n=37)`,udp_1) #10:edih

#get coefficient for each signature
amed_0 <- amed1
amed_1 <- amed1[which(amed1$`AMED (n=48)` %in% c("(Intercept)",hdp_1)),]; rownames(amed_1) <- amed_1$`AMED (n=48)`
amed_2 <- amed1[which(amed1$`AMED (n=48)` %in% c("(Intercept)",hdp_2)),]; rownames(amed_2) <- amed_2$`AMED (n=48)`
amed_3 <- amed1[which(amed1$`AMED (n=48)` %in% c("(Intercept)",hdp_3)),]; rownames(amed_3) <- amed_3$`AMED (n=48)`
amed_4 <- amed1[which(amed1$`AMED (n=48)` %in% c("(Intercept)",hdp_5)),]; rownames(amed_4) <- amed_4$`AMED (n=48)`

ahei_0 <- ahei1
ahei_1 <- ahei1[which(ahei1$`AHEI-2010 (n=61)` %in% c("(Intercept)",hdp_1)),]; rownames(ahei_1) <- ahei_1$`AHEI-2010 (n=61)`
ahei_2 <- ahei1[which(ahei1$`AHEI-2010 (n=61)` %in% c("(Intercept)",hdp_2)),]; rownames(ahei_2) <- ahei_2$`AHEI-2010 (n=61)`
ahei_3 <- ahei1[which(ahei1$`AHEI-2010 (n=61)` %in% c("(Intercept)",hdp_4)),]; rownames(ahei_3) <- ahei_3$`AHEI-2010 (n=61)`
ahei_4 <- ahei1[which(ahei1$`AHEI-2010 (n=61)` %in% c("(Intercept)",hdp_6)),]; rownames(ahei_4) <- ahei_4$`AHEI-2010 (n=61)`

dash_0 <- dash1
dash_1 <- dash1[which(dash1$`DASH (n=55)` %in% c("(Intercept)",hdp_1)),]; rownames(dash_1) <- dash_1$`DASH (n=55)`
dash_2 <- dash1[which(dash1$`DASH (n=55)` %in% c("(Intercept)",hdp_3)),]; rownames(dash_2) <- dash_2$`DASH (n=55)`
dash_3 <- dash1[which(dash1$`DASH (n=55)` %in% c("(Intercept)",hdp_4)),]; rownames(dash_3) <- dash_3$`DASH (n=55)`
dash_4 <- dash1[which(dash1$`DASH (n=55)` %in% c("(Intercept)",hdp_7)),]; rownames(dash_4) <- dash_4$`DASH (n=55)`

pdi_0 <- opdi1
pdi_1 <- opdi1[which(opdi1$`oPDI (n=36)` %in% c("(Intercept)",pdp_1)),]; rownames(pdi_1) <- pdi_1$`oPDI (n=36)`
pdi_2 <- opdi1[which(opdi1$`oPDI (n=36)` %in% c("(Intercept)",pdp_2)),]; rownames(pdi_2) <- pdi_2$`oPDI (n=36)`
pdi_3 <- opdi1[which(opdi1$`oPDI (n=36)` %in% c("(Intercept)",pdp_3)),]; rownames(pdi_3) <- pdi_3$`oPDI (n=36)`
pdi_4 <- opdi1[which(opdi1$`oPDI (n=36)` %in% c("(Intercept)",pdp_5)),]; rownames(pdi_4) <- pdi_4$`oPDI (n=36)`

hpdi_0 <- hpdi1
hpdi_1 <- hpdi1[which(hpdi1$`hPDI (n=39)` %in% c("(Intercept)",pdp_1)),]; rownames(hpdi_1) <- hpdi_1$`hPDI (n=39)` 
hpdi_2 <- hpdi1[which(hpdi1$`hPDI (n=39)` %in% c("(Intercept)",pdp_2)),]; rownames(hpdi_2) <- hpdi_2$`hPDI (n=39)` 
hpdi_3 <- hpdi1[which(hpdi1$`hPDI (n=39)` %in% c("(Intercept)",pdp_4)),]; rownames(hpdi_3) <- hpdi_3$`hPDI (n=39)` 
hpdi_4 <- hpdi1[which(hpdi1$`hPDI (n=39)` %in% c("(Intercept)",pdp_6)),]; rownames(hpdi_4) <- hpdi_4$`hPDI (n=39)` 

updi_0 <- updi1
updi_1 <- updi1[which(updi1$`uPDI (n=43)` %in% c("(Intercept)",pdp_1)),]; rownames(updi_1) <-  updi_1$`uPDI (n=43)`
updi_2 <- updi1[which(updi1$`uPDI (n=43)` %in% c("(Intercept)",pdp_3)),]; rownames(updi_2) <-  updi_2$`uPDI (n=43)`
updi_3 <- updi1[which(updi1$`uPDI (n=43)` %in% c("(Intercept)",pdp_4)),]; rownames(updi_3) <-  updi_3$`uPDI (n=43)`
updi_4 <- updi1[which(updi1$`uPDI (n=43)` %in% c("(Intercept)",pdp_7)),]; rownames(updi_4) <-  updi_4$`uPDI (n=43)` 

edip_0 <- edip1
edip_1 <- edip1[which(edip1$`EDIP (n=66)` %in% c("(Intercept)",udp_1)),]; rownames(edip_1) <- edip_1$`EDIP (n=66)` 
edip_2 <- edip1[which(edip1$`EDIP (n=66)` %in% c("(Intercept)",udp_2)),]; rownames(edip_2) <- edip_2$`EDIP (n=66)` 

edih_0 <- edih1
edih_1 <- edih1[which(edih1$`EDIH (n=37)` %in% c("(Intercept)",udp_1)),]; rownames(edih_1) <- edih_1$`EDIH (n=37)` 
edih_2 <- edih1[which(edih1$`EDIH (n=37)` %in% c("(Intercept)",udp_3)),]; rownames(edih_1) <- edih_1$`EDIH (n=37)` 

#--------------------------------------------------------------------------------------------
#
#           Step2: calculate shared metabolic signature and distinct metabolic signature in NHS/HPFS
#
#--------------------------------------------------------------------------------------------
#load nhs/hpfs metabolome data
load("Extracted_NHS_HPFS_20221225.RData")

#split metabolic signature for each dietary score
amed1 <- na.omit(sig_list[,1:2]); amed <- data.frame(amed1[,-1]); rownames(amed) <- amed1[,1]
ahei1 <- na.omit(sig_list[,4:5]); ahei <- data.frame(ahei1[,-1]); rownames(ahei) <- ahei1[,1]
dash1 <- na.omit(sig_list[,7:8]); dash <- data.frame(dash1[,-1]); rownames(dash) <- dash1[,1]
opdi1 <- na.omit(sig_list[,10:11]); opdi <- data.frame(opdi1[,-1]); rownames(opdi) <- opdi1[,1]
hpdi1 <- na.omit(sig_list[,13:14]); hpdi <- data.frame(hpdi1[,-1]); rownames(hpdi) <- hpdi1[,1]
updi1 <- na.omit(sig_list[,16:17]); updi <- data.frame(updi1[,-1]); rownames(updi) <- updi1[,1]
edip1 <- na.omit(sig_list[,19:20]); edip <- data.frame(edip1[,-1]); rownames(edip) <- edip1[,1]
edih1 <- na.omit(sig_list[,22:23]); edih <- data.frame(edih1[,-1]); rownames(edih) <- edih1[,1]

#met3:cross-platfrom
amed <- data.frame(amed_0[,2]); rownames(amed) <- amed_0$`AMED (n=48)`
ahei <- data.frame(ahei_0[,2]); rownames(ahei) <- ahei_0$`AHEI-2010 (n=61)`
dash <- data.frame(dash_0[,2]); rownames(dash) <- dash_0$`DASH (n=55)`
opdi <- data.frame(pdi_0[,2]); rownames(opdi) <- pdi_0$`oPDI (n=36)`
hpdi <- data.frame(hpdi_0[,2]); rownames(hpdi) <- hpdi_0$`hPDI (n=39)`
updi <- data.frame(updi_0[,2]); rownames(updi) <- updi_0$`uPDI (n=43)`
edip <- data.frame(edip_0[,2]); rownames(edip) <- edip_0$`EDIP (n=66)`
edih <- data.frame(edih_0[,2]); rownames(edih) <- edih_0$`EDIH (n=37)`

#calculate metabolic signature in each sub-study
#nhs1.stroke.std
nhs1.stroke.std$amed = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.stroke.std$ahei = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.stroke.std$dash = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.stroke.std$opdi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.stroke.std$hpdi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.stroke.std$updi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.stroke.std$edip = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.stroke.std$edih = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.als.std
nhs1.als.std$amed = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.als.std$ahei = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.als.std$dash = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.als.std$opdi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.als.std$hpdi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.als.std$updi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.als.std$edip = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.als.std$edih = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.breast.std
nhs1.breast.std$amed = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.breast.std$ahei = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.breast.std$dash = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.breast.std$opdi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.breast.std$hpdi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.breast.std$updi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.breast.std$edip = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.breast.std$edih = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.colon.std
nhs1.colon.std$amed = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.colon.std$ahei = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.colon.std$dash = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.colon.std$opdi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.colon.std$hpdi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.colon.std$updi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.colon.std$edip = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.colon.std$edih = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.diabetes.std
nhs1.diabetes.std$amed = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.diabetes.std$ahei = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.diabetes.std$dash = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.diabetes.std$opdi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.diabetes.std$hpdi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.diabetes.std$updi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.diabetes.std$edip = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.diabetes.std$edih = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.exfoliation.glaucoma.std
nhs1.exfoliation.glaucoma.std$amed = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.exfoliation.glaucoma.std$ahei = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.exfoliation.glaucoma.std$dash = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.exfoliation.glaucoma.std$opdi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.exfoliation.glaucoma.std$hpdi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.exfoliation.glaucoma.std$updi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.exfoliation.glaucoma.std$edip = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.exfoliation.glaucoma.std$edih = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.ibd.std
nhs1.ibd.std$amed = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.ibd.std$ahei = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.ibd.std$dash = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.ibd.std$opdi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.ibd.std$hpdi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.ibd.std$updi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.ibd.std$edip = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.ibd.std$edih = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.ovarian.std
nhs1.ovarian.std$amed = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.ovarian.std$ahei = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.ovarian.std$dash = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.ovarian.std$opdi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.ovarian.std$hpdi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.ovarian.std$updi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.ovarian.std$edip = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.ovarian.std$edih = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.parkinsons.std
nhs1.parkinsons.std$amed = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.parkinsons.std$ahei = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.parkinsons.std$dash = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.parkinsons.std$opdi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.parkinsons.std$hpdi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.parkinsons.std$updi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.parkinsons.std$edip = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.parkinsons.std$edih = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.poag.std
nhs1.poag.std$amed = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.poag.std$ahei = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.poag.std$dash = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.poag.std$opdi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.poag.std$hpdi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.poag.std$updi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.poag.std$edip = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.poag.std$edih = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.racial.diff.std
nhs1.racial.diff.std$amed = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.racial.diff.std$ahei = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.racial.diff.std$dash = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.racial.diff.std$opdi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.racial.diff.std$hpdi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.racial.diff.std$updi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.racial.diff.std$edip = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.racial.diff.std$edih = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.rheumatoid.std
nhs1.rheumatoid.std$amed = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.rheumatoid.std$ahei = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.rheumatoid.std$dash = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.rheumatoid.std$opdi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.rheumatoid.std$hpdi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.rheumatoid.std$updi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.rheumatoid.std$edip = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.rheumatoid.std$edih = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.breast.std
nhs2.breast.std$amed = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.breast.std$ahei = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.breast.std$dash = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.breast.std$opdi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.breast.std$hpdi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.breast.std$updi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.breast.std$edip = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.breast.std$edih = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.ibd.std
nhs2.ibd.std$amed = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.ibd.std$ahei = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.ibd.std$dash = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.ibd.std$opdi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.ibd.std$hpdi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.ibd.std$updi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.ibd.std$edip = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.ibd.std$edih = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.diabetes.std
nhs2.diabetes.std$amed = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.diabetes.std$ahei = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.diabetes.std$dash = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.diabetes.std$opdi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.diabetes.std$hpdi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.diabetes.std$updi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.diabetes.std$edip = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.diabetes.std$edih = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.ovarian.std
nhs2.ovarian.std$amed = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.ovarian.std$ahei = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.ovarian.std$dash = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.ovarian.std$opdi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.ovarian.std$hpdi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.ovarian.std$updi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.ovarian.std$edip = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.ovarian.std$edih = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.poag.std
nhs2.poag.std$amed = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.poag.std$ahei = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.poag.std$dash = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.poag.std$opdi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.poag.std$hpdi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.poag.std$updi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.poag.std$edip = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.poag.std$edih = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.rheumatoid.std
nhs2.rheumatoid.std$amed = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.rheumatoid.std$ahei = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.rheumatoid.std$dash = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.rheumatoid.std$opdi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.rheumatoid.std$hpdi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.rheumatoid.std$updi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.rheumatoid.std$edip = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.rheumatoid.std$edih = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.stress.std
nhs2.stress.std$amed = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.stress.std$ahei = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.stress.std$dash = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.stress.std$opdi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.stress.std$hpdi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.stress.std$updi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.stress.std$edip = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.stress.std$edih = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.stroke.std
nhs2.stroke.std$amed = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.stroke.std$ahei = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.stroke.std$dash = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.stroke.std$opdi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.stroke.std$hpdi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.stroke.std$updi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.stroke.std$edip = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.stroke.std$edih = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.als.std
hpfs.als.std$amed = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.als.std$ahei = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.als.std$dash = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.als.std$opdi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.als.std$hpdi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.als.std$updi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.als.std$edip = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.als.std$edih = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.colon.std
hpfs.colon.std$amed = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.colon.std$ahei = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.colon.std$dash = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.colon.std$opdi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.colon.std$hpdi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.colon.std$updi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.colon.std$edip = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.colon.std$edih = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.exfoliation.glaucoma.std
hpfs.exfoliation.glaucoma.std$amed = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.exfoliation.glaucoma.std$ahei = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.exfoliation.glaucoma.std$dash = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.exfoliation.glaucoma.std$opdi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.exfoliation.glaucoma.std$hpdi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.exfoliation.glaucoma.std$updi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.exfoliation.glaucoma.std$edip = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.exfoliation.glaucoma.std$edih = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.parkinsons.std
hpfs.parkinsons.std$amed = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.parkinsons.std$ahei = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.parkinsons.std$dash = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.parkinsons.std$opdi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.parkinsons.std$hpdi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.parkinsons.std$updi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.parkinsons.std$edip = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.parkinsons.std$edih = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.poag.std
hpfs.poag.std$amed = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.poag.std$ahei = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.poag.std$dash = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.poag.std$opdi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.poag.std$hpdi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.poag.std$updi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.poag.std$edip = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.poag.std$edih = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#merge nhs1, nhs2, and hpfs sub-studies (n=13,258)
sig <- c("amed","ahei","dash","opdi","hpdi","updi","edip","edih")

nhs1_met3 <- rbind(nhs1.stroke.std[sig], nhs1.als.std[sig], nhs1.breast.std[sig], nhs1.colon.std[sig], nhs1.diabetes.std[sig],
                   nhs1.exfoliation.glaucoma.std[sig], nhs1.ibd.std[sig], nhs1.ovarian.std[sig], nhs1.parkinsons.std[sig],
                   nhs1.poag.std[sig], nhs1.racial.diff.std[sig],nhs1.rheumatoid.std[sig]) #12 sub-studies, 8134

nhs2_met3 <- rbind(nhs2.breast.std[sig], nhs2.ibd.std[sig], nhs2.diabetes.std[sig], nhs2.ovarian.std[sig], nhs2.poag.std[sig],
                   nhs2.rheumatoid.std[sig], nhs2.stress.std[sig], nhs2.stroke.std[sig]) #8 sub-studies, 3473

hpfs_met3 <- rbind(hpfs.als.std[sig], hpfs.colon.std[sig], hpfs.exfoliation.glaucoma.std[sig], hpfs.parkinsons.std[sig],
                   hpfs.poag.std[sig]) #5 sub-studies, 1651

nhs1_met3$id <- rownames(nhs1_met3); nhs2_met3$id <- rownames(nhs2_met3); hpfs_met3$id <- rownames(hpfs_met3)

#met4:overlap among three
amed <- data.frame(amed_1[,2]); rownames(amed) <- amed_1$`AMED (n=48)`
ahei <- data.frame(ahei_1[,2]); rownames(ahei) <- ahei_1$`AHEI-2010 (n=61)`
dash <- data.frame(dash_1[,2]); rownames(dash) <- dash_1$`DASH (n=55)`
opdi <- data.frame(pdi_1[,2]); rownames(opdi) <- pdi_1$`oPDI (n=36)`
hpdi <- data.frame(hpdi_1[,2]); rownames(hpdi) <- hpdi_1$`hPDI (n=39)`
updi <- data.frame(updi_1[,2]); rownames(updi) <- updi_1$`uPDI (n=43)`
edip <- data.frame(edip_1[,2]); rownames(edip) <- edip_1$`EDIP (n=66)`
edih <- data.frame(edih_1[,2]); rownames(edih) <- edih_1$`EDIH (n=37)`

#calculate metabolic signature in each sub-study
#nhs1.stroke.std
nhs1.stroke.std$amed = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.stroke.std$ahei = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.stroke.std$dash = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.stroke.std$opdi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.stroke.std$hpdi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.stroke.std$updi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.stroke.std$edip = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.stroke.std$edih = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.als.std
nhs1.als.std$amed = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.als.std$ahei = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.als.std$dash = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.als.std$opdi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.als.std$hpdi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.als.std$updi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.als.std$edip = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.als.std$edih = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.breast.std
nhs1.breast.std$amed = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.breast.std$ahei = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.breast.std$dash = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.breast.std$opdi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.breast.std$hpdi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.breast.std$updi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.breast.std$edip = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.breast.std$edih = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.colon.std
nhs1.colon.std$amed = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.colon.std$ahei = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.colon.std$dash = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.colon.std$opdi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.colon.std$hpdi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.colon.std$updi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.colon.std$edip = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.colon.std$edih = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.diabetes.std
nhs1.diabetes.std$amed = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.diabetes.std$ahei = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.diabetes.std$dash = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.diabetes.std$opdi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.diabetes.std$hpdi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.diabetes.std$updi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.diabetes.std$edip = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.diabetes.std$edih = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.exfoliation.glaucoma.std
nhs1.exfoliation.glaucoma.std$amed = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.exfoliation.glaucoma.std$ahei = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.exfoliation.glaucoma.std$dash = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.exfoliation.glaucoma.std$opdi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.exfoliation.glaucoma.std$hpdi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.exfoliation.glaucoma.std$updi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.exfoliation.glaucoma.std$edip = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.exfoliation.glaucoma.std$edih = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.ibd.std
nhs1.ibd.std$amed = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.ibd.std$ahei = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.ibd.std$dash = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.ibd.std$opdi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.ibd.std$hpdi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.ibd.std$updi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.ibd.std$edip = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.ibd.std$edih = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.ovarian.std
nhs1.ovarian.std$amed = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.ovarian.std$ahei = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.ovarian.std$dash = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.ovarian.std$opdi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.ovarian.std$hpdi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.ovarian.std$updi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.ovarian.std$edip = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.ovarian.std$edih = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.parkinsons.std
nhs1.parkinsons.std$amed = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.parkinsons.std$ahei = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.parkinsons.std$dash = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.parkinsons.std$opdi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.parkinsons.std$hpdi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.parkinsons.std$updi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.parkinsons.std$edip = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.parkinsons.std$edih = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.poag.std
nhs1.poag.std$amed = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.poag.std$ahei = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.poag.std$dash = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.poag.std$opdi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.poag.std$hpdi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.poag.std$updi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.poag.std$edip = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.poag.std$edih = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.racial.diff.std
nhs1.racial.diff.std$amed = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.racial.diff.std$ahei = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.racial.diff.std$dash = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.racial.diff.std$opdi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.racial.diff.std$hpdi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.racial.diff.std$updi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.racial.diff.std$edip = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.racial.diff.std$edih = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.rheumatoid.std
nhs1.rheumatoid.std$amed = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.rheumatoid.std$ahei = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.rheumatoid.std$dash = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.rheumatoid.std$opdi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.rheumatoid.std$hpdi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.rheumatoid.std$updi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.rheumatoid.std$edip = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.rheumatoid.std$edih = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.breast.std
nhs2.breast.std$amed = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.breast.std$ahei = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.breast.std$dash = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.breast.std$opdi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.breast.std$hpdi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.breast.std$updi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.breast.std$edip = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.breast.std$edih = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.ibd.std
nhs2.ibd.std$amed = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.ibd.std$ahei = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.ibd.std$dash = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.ibd.std$opdi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.ibd.std$hpdi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.ibd.std$updi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.ibd.std$edip = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.ibd.std$edih = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.diabetes.std
nhs2.diabetes.std$amed = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.diabetes.std$ahei = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.diabetes.std$dash = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.diabetes.std$opdi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.diabetes.std$hpdi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.diabetes.std$updi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.diabetes.std$edip = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.diabetes.std$edih = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.ovarian.std
nhs2.ovarian.std$amed = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.ovarian.std$ahei = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.ovarian.std$dash = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.ovarian.std$opdi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.ovarian.std$hpdi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.ovarian.std$updi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.ovarian.std$edip = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.ovarian.std$edih = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.poag.std
nhs2.poag.std$amed = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.poag.std$ahei = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.poag.std$dash = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.poag.std$opdi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.poag.std$hpdi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.poag.std$updi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.poag.std$edip = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.poag.std$edih = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.rheumatoid.std
nhs2.rheumatoid.std$amed = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.rheumatoid.std$ahei = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.rheumatoid.std$dash = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.rheumatoid.std$opdi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.rheumatoid.std$hpdi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.rheumatoid.std$updi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.rheumatoid.std$edip = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.rheumatoid.std$edih = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.stress.std
nhs2.stress.std$amed = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.stress.std$ahei = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.stress.std$dash = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.stress.std$opdi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.stress.std$hpdi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.stress.std$updi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.stress.std$edip = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.stress.std$edih = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.stroke.std
nhs2.stroke.std$amed = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.stroke.std$ahei = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.stroke.std$dash = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.stroke.std$opdi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.stroke.std$hpdi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.stroke.std$updi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.stroke.std$edip = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.stroke.std$edih = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.als.std
hpfs.als.std$amed = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.als.std$ahei = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.als.std$dash = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.als.std$opdi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.als.std$hpdi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.als.std$updi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.als.std$edip = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.als.std$edih = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.colon.std
hpfs.colon.std$amed = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.colon.std$ahei = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.colon.std$dash = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.colon.std$opdi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.colon.std$hpdi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.colon.std$updi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.colon.std$edip = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.colon.std$edih = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.exfoliation.glaucoma.std
hpfs.exfoliation.glaucoma.std$amed = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.exfoliation.glaucoma.std$ahei = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.exfoliation.glaucoma.std$dash = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.exfoliation.glaucoma.std$opdi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.exfoliation.glaucoma.std$hpdi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.exfoliation.glaucoma.std$updi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.exfoliation.glaucoma.std$edip = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.exfoliation.glaucoma.std$edih = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.parkinsons.std
hpfs.parkinsons.std$amed = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.parkinsons.std$ahei = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.parkinsons.std$dash = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.parkinsons.std$opdi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.parkinsons.std$hpdi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.parkinsons.std$updi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.parkinsons.std$edip = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.parkinsons.std$edih = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.poag.std
hpfs.poag.std$amed = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.poag.std$ahei = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.poag.std$dash = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.poag.std$opdi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.poag.std$hpdi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.poag.std$updi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.poag.std$edip = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.poag.std$edih = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#merge nhs1, nhs2, and hpfs sub-studies (n=13,258)
sig <- c("amed","ahei","dash","opdi","hpdi","updi","edip","edih")

nhs1_met4 <- rbind(nhs1.stroke.std[sig], nhs1.als.std[sig], nhs1.breast.std[sig], nhs1.colon.std[sig], nhs1.diabetes.std[sig],
                   nhs1.exfoliation.glaucoma.std[sig], nhs1.ibd.std[sig], nhs1.ovarian.std[sig], nhs1.parkinsons.std[sig],
                   nhs1.poag.std[sig], nhs1.racial.diff.std[sig],nhs1.rheumatoid.std[sig]) #12 sub-studies, 8134

nhs2_met4 <- rbind(nhs2.breast.std[sig], nhs2.ibd.std[sig], nhs2.diabetes.std[sig], nhs2.ovarian.std[sig], nhs2.poag.std[sig],
                   nhs2.rheumatoid.std[sig], nhs2.stress.std[sig], nhs2.stroke.std[sig]) #8 sub-studies, 3473

hpfs_met4 <- rbind(hpfs.als.std[sig], hpfs.colon.std[sig], hpfs.exfoliation.glaucoma.std[sig], hpfs.parkinsons.std[sig],
                   hpfs.poag.std[sig]) #5 sub-studies, 1651

nhs1_met4$id <- rownames(nhs1_met4); nhs2_met4$id <- rownames(nhs2_met4); hpfs_met4$id <- rownames(hpfs_met4)

#met5:overlap between two
amed <- data.frame(amed_2[,2]); rownames(amed) <- amed_2$`AMED (n=48)`
ahei <- data.frame(ahei_2[,2]); rownames(ahei) <- ahei_2$`AHEI-2010 (n=61)`
dash <- data.frame(dash_2[,2]); rownames(dash) <- dash_2$`DASH (n=55)`
opdi <- data.frame(pdi_2[,2]); rownames(opdi) <- pdi_2$`oPDI (n=36)`
hpdi <- data.frame(hpdi_2[,2]); rownames(hpdi) <- hpdi_2$`hPDI (n=39)`
updi <- data.frame(updi_2[,2]); rownames(updi) <- updi_2$`uPDI (n=43)`
edip <- data.frame(edip_2[,2]); rownames(edip) <- edip_2$`EDIP (n=66)`
edih <- data.frame(edih_2[,2]); rownames(edih) <- edih_2$`EDIH (n=37)`

#calculate metabolic signature in each sub-study
#nhs1.stroke.std
nhs1.stroke.std$amed = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.stroke.std$ahei = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.stroke.std$dash = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.stroke.std$opdi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.stroke.std$hpdi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.stroke.std$updi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.stroke.std$edip = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.stroke.std$edih = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.als.std
nhs1.als.std$amed = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.als.std$ahei = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.als.std$dash = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.als.std$opdi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.als.std$hpdi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.als.std$updi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.als.std$edip = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.als.std$edih = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.breast.std
nhs1.breast.std$amed = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.breast.std$ahei = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.breast.std$dash = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.breast.std$opdi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.breast.std$hpdi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.breast.std$updi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.breast.std$edip = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.breast.std$edih = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.colon.std
nhs1.colon.std$amed = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.colon.std$ahei = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.colon.std$dash = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.colon.std$opdi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.colon.std$hpdi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.colon.std$updi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.colon.std$edip = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.colon.std$edih = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.diabetes.std
nhs1.diabetes.std$amed = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.diabetes.std$ahei = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.diabetes.std$dash = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.diabetes.std$opdi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.diabetes.std$hpdi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.diabetes.std$updi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.diabetes.std$edip = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.diabetes.std$edih = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.exfoliation.glaucoma.std
nhs1.exfoliation.glaucoma.std$amed = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.exfoliation.glaucoma.std$ahei = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.exfoliation.glaucoma.std$dash = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.exfoliation.glaucoma.std$opdi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.exfoliation.glaucoma.std$hpdi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.exfoliation.glaucoma.std$updi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.exfoliation.glaucoma.std$edip = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.exfoliation.glaucoma.std$edih = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.ibd.std
nhs1.ibd.std$amed = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.ibd.std$ahei = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.ibd.std$dash = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.ibd.std$opdi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.ibd.std$hpdi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.ibd.std$updi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.ibd.std$edip = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.ibd.std$edih = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.ovarian.std
nhs1.ovarian.std$amed = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.ovarian.std$ahei = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.ovarian.std$dash = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.ovarian.std$opdi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.ovarian.std$hpdi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.ovarian.std$updi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.ovarian.std$edip = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.ovarian.std$edih = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.parkinsons.std
nhs1.parkinsons.std$amed = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.parkinsons.std$ahei = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.parkinsons.std$dash = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.parkinsons.std$opdi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.parkinsons.std$hpdi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.parkinsons.std$updi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.parkinsons.std$edip = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.parkinsons.std$edih = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.poag.std
nhs1.poag.std$amed = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.poag.std$ahei = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.poag.std$dash = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.poag.std$opdi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.poag.std$hpdi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.poag.std$updi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.poag.std$edip = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.poag.std$edih = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.racial.diff.std
nhs1.racial.diff.std$amed = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.racial.diff.std$ahei = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.racial.diff.std$dash = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.racial.diff.std$opdi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.racial.diff.std$hpdi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.racial.diff.std$updi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.racial.diff.std$edip = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.racial.diff.std$edih = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.rheumatoid.std
nhs1.rheumatoid.std$amed = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.rheumatoid.std$ahei = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.rheumatoid.std$dash = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.rheumatoid.std$opdi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.rheumatoid.std$hpdi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.rheumatoid.std$updi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.rheumatoid.std$edip = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.rheumatoid.std$edih = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.breast.std
nhs2.breast.std$amed = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.breast.std$ahei = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.breast.std$dash = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.breast.std$opdi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.breast.std$hpdi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.breast.std$updi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.breast.std$edip = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.breast.std$edih = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.ibd.std
nhs2.ibd.std$amed = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.ibd.std$ahei = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.ibd.std$dash = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.ibd.std$opdi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.ibd.std$hpdi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.ibd.std$updi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.ibd.std$edip = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.ibd.std$edih = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.diabetes.std
nhs2.diabetes.std$amed = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.diabetes.std$ahei = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.diabetes.std$dash = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.diabetes.std$opdi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.diabetes.std$hpdi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.diabetes.std$updi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.diabetes.std$edip = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.diabetes.std$edih = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.ovarian.std
nhs2.ovarian.std$amed = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.ovarian.std$ahei = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.ovarian.std$dash = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.ovarian.std$opdi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.ovarian.std$hpdi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.ovarian.std$updi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.ovarian.std$edip = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.ovarian.std$edih = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.poag.std
nhs2.poag.std$amed = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.poag.std$ahei = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.poag.std$dash = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.poag.std$opdi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.poag.std$hpdi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.poag.std$updi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.poag.std$edip = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.poag.std$edih = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.rheumatoid.std
nhs2.rheumatoid.std$amed = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.rheumatoid.std$ahei = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.rheumatoid.std$dash = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.rheumatoid.std$opdi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.rheumatoid.std$hpdi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.rheumatoid.std$updi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.rheumatoid.std$edip = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.rheumatoid.std$edih = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.stress.std
nhs2.stress.std$amed = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.stress.std$ahei = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.stress.std$dash = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.stress.std$opdi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.stress.std$hpdi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.stress.std$updi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.stress.std$edip = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.stress.std$edih = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.stroke.std
nhs2.stroke.std$amed = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.stroke.std$ahei = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.stroke.std$dash = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.stroke.std$opdi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.stroke.std$hpdi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.stroke.std$updi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.stroke.std$edip = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.stroke.std$edih = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.als.std
hpfs.als.std$amed = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.als.std$ahei = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.als.std$dash = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.als.std$opdi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.als.std$hpdi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.als.std$updi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.als.std$edip = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.als.std$edih = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.colon.std
hpfs.colon.std$amed = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.colon.std$ahei = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.colon.std$dash = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.colon.std$opdi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.colon.std$hpdi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.colon.std$updi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.colon.std$edip = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.colon.std$edih = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.exfoliation.glaucoma.std
hpfs.exfoliation.glaucoma.std$amed = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.exfoliation.glaucoma.std$ahei = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.exfoliation.glaucoma.std$dash = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.exfoliation.glaucoma.std$opdi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.exfoliation.glaucoma.std$hpdi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.exfoliation.glaucoma.std$updi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.exfoliation.glaucoma.std$edip = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.exfoliation.glaucoma.std$edih = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.parkinsons.std
hpfs.parkinsons.std$amed = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.parkinsons.std$ahei = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.parkinsons.std$dash = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.parkinsons.std$opdi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.parkinsons.std$hpdi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.parkinsons.std$updi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.parkinsons.std$edip = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.parkinsons.std$edih = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.poag.std
hpfs.poag.std$amed = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.poag.std$ahei = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.poag.std$dash = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.poag.std$opdi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.poag.std$hpdi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.poag.std$updi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.poag.std$edip = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.poag.std$edih = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#merge nhs1, nhs2, and hpfs sub-studies (n=13,258)
sig <- c("amed","ahei","dash","opdi","hpdi","updi","edip","edih")

nhs1_met5 <- rbind(nhs1.stroke.std[sig], nhs1.als.std[sig], nhs1.breast.std[sig], nhs1.colon.std[sig], nhs1.diabetes.std[sig],
                   nhs1.exfoliation.glaucoma.std[sig], nhs1.ibd.std[sig], nhs1.ovarian.std[sig], nhs1.parkinsons.std[sig],
                   nhs1.poag.std[sig], nhs1.racial.diff.std[sig],nhs1.rheumatoid.std[sig]) #12 sub-studies, 8134

nhs2_met5 <- rbind(nhs2.breast.std[sig], nhs2.ibd.std[sig], nhs2.diabetes.std[sig], nhs2.ovarian.std[sig], nhs2.poag.std[sig],
                   nhs2.rheumatoid.std[sig], nhs2.stress.std[sig], nhs2.stroke.std[sig]) #8 sub-studies, 3473

hpfs_met5 <- rbind(hpfs.als.std[sig], hpfs.colon.std[sig], hpfs.exfoliation.glaucoma.std[sig], hpfs.parkinsons.std[sig],
                   hpfs.poag.std[sig]) #5 sub-studies, 1651

nhs1_met5$id <- rownames(nhs1_met5); nhs2_met5$id <- rownames(nhs2_met5); hpfs_met5$id <- rownames(hpfs_met5)

#met6:overlap between two
amed <- data.frame(amed_3[,2]); rownames(amed) <- amed_3$`AMED (n=48)`
ahei <- data.frame(ahei_3[,2]); rownames(ahei) <- ahei_3$`AHEI-2010 (n=61)`
dash <- data.frame(dash_3[,2]); rownames(dash) <- dash_3$`DASH (n=55)`
opdi <- data.frame(pdi_3[,2]); rownames(opdi) <- pdi_3$`oPDI (n=36)`
hpdi <- data.frame(hpdi_3[,2]); rownames(hpdi) <- hpdi_3$`hPDI (n=39)`
updi <- data.frame(updi_3[,2]); rownames(updi) <- updi_3$`uPDI (n=43)`
edip <- data.frame(edip_1[,2]); rownames(edip) <- edip_1$`EDIP (n=66)`
edih <- data.frame(edih_1[,2]); rownames(edih) <- edih_1$`EDIH (n=37)`

#calculate metabolic signature in each sub-study
#nhs1.stroke.std
nhs1.stroke.std$amed = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.stroke.std$ahei = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.stroke.std$dash = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.stroke.std$opdi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.stroke.std$hpdi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.stroke.std$updi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.stroke.std$edip = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.stroke.std$edih = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.als.std
nhs1.als.std$amed = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.als.std$ahei = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.als.std$dash = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.als.std$opdi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.als.std$hpdi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.als.std$updi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.als.std$edip = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.als.std$edih = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.breast.std
nhs1.breast.std$amed = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.breast.std$ahei = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.breast.std$dash = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.breast.std$opdi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.breast.std$hpdi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.breast.std$updi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.breast.std$edip = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.breast.std$edih = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.colon.std
nhs1.colon.std$amed = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.colon.std$ahei = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.colon.std$dash = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.colon.std$opdi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.colon.std$hpdi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.colon.std$updi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.colon.std$edip = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.colon.std$edih = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.diabetes.std
nhs1.diabetes.std$amed = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.diabetes.std$ahei = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.diabetes.std$dash = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.diabetes.std$opdi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.diabetes.std$hpdi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.diabetes.std$updi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.diabetes.std$edip = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.diabetes.std$edih = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.exfoliation.glaucoma.std
nhs1.exfoliation.glaucoma.std$amed = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.exfoliation.glaucoma.std$ahei = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.exfoliation.glaucoma.std$dash = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.exfoliation.glaucoma.std$opdi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.exfoliation.glaucoma.std$hpdi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.exfoliation.glaucoma.std$updi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.exfoliation.glaucoma.std$edip = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.exfoliation.glaucoma.std$edih = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.ibd.std
nhs1.ibd.std$amed = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.ibd.std$ahei = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.ibd.std$dash = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.ibd.std$opdi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.ibd.std$hpdi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.ibd.std$updi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.ibd.std$edip = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.ibd.std$edih = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.ovarian.std
nhs1.ovarian.std$amed = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.ovarian.std$ahei = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.ovarian.std$dash = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.ovarian.std$opdi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.ovarian.std$hpdi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.ovarian.std$updi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.ovarian.std$edip = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.ovarian.std$edih = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.parkinsons.std
nhs1.parkinsons.std$amed = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.parkinsons.std$ahei = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.parkinsons.std$dash = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.parkinsons.std$opdi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.parkinsons.std$hpdi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.parkinsons.std$updi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.parkinsons.std$edip = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.parkinsons.std$edih = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.poag.std
nhs1.poag.std$amed = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.poag.std$ahei = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.poag.std$dash = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.poag.std$opdi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.poag.std$hpdi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.poag.std$updi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.poag.std$edip = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.poag.std$edih = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.racial.diff.std
nhs1.racial.diff.std$amed = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.racial.diff.std$ahei = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.racial.diff.std$dash = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.racial.diff.std$opdi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.racial.diff.std$hpdi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.racial.diff.std$updi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.racial.diff.std$edip = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.racial.diff.std$edih = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.rheumatoid.std
nhs1.rheumatoid.std$amed = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.rheumatoid.std$ahei = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.rheumatoid.std$dash = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.rheumatoid.std$opdi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.rheumatoid.std$hpdi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.rheumatoid.std$updi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.rheumatoid.std$edip = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.rheumatoid.std$edih = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.breast.std
nhs2.breast.std$amed = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.breast.std$ahei = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.breast.std$dash = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.breast.std$opdi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.breast.std$hpdi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.breast.std$updi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.breast.std$edip = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.breast.std$edih = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.ibd.std
nhs2.ibd.std$amed = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.ibd.std$ahei = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.ibd.std$dash = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.ibd.std$opdi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.ibd.std$hpdi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.ibd.std$updi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.ibd.std$edip = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.ibd.std$edih = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.diabetes.std
nhs2.diabetes.std$amed = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.diabetes.std$ahei = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.diabetes.std$dash = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.diabetes.std$opdi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.diabetes.std$hpdi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.diabetes.std$updi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.diabetes.std$edip = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.diabetes.std$edih = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.ovarian.std
nhs2.ovarian.std$amed = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.ovarian.std$ahei = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.ovarian.std$dash = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.ovarian.std$opdi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.ovarian.std$hpdi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.ovarian.std$updi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.ovarian.std$edip = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.ovarian.std$edih = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.poag.std
nhs2.poag.std$amed = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.poag.std$ahei = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.poag.std$dash = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.poag.std$opdi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.poag.std$hpdi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.poag.std$updi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.poag.std$edip = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.poag.std$edih = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.rheumatoid.std
nhs2.rheumatoid.std$amed = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.rheumatoid.std$ahei = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.rheumatoid.std$dash = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.rheumatoid.std$opdi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.rheumatoid.std$hpdi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.rheumatoid.std$updi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.rheumatoid.std$edip = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.rheumatoid.std$edih = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.stress.std
nhs2.stress.std$amed = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.stress.std$ahei = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.stress.std$dash = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.stress.std$opdi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.stress.std$hpdi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.stress.std$updi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.stress.std$edip = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.stress.std$edih = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.stroke.std
nhs2.stroke.std$amed = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.stroke.std$ahei = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.stroke.std$dash = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.stroke.std$opdi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.stroke.std$hpdi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.stroke.std$updi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.stroke.std$edip = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.stroke.std$edih = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.als.std
hpfs.als.std$amed = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.als.std$ahei = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.als.std$dash = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.als.std$opdi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.als.std$hpdi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.als.std$updi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.als.std$edip = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.als.std$edih = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.colon.std
hpfs.colon.std$amed = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.colon.std$ahei = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.colon.std$dash = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.colon.std$opdi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.colon.std$hpdi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.colon.std$updi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.colon.std$edip = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.colon.std$edih = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.exfoliation.glaucoma.std
hpfs.exfoliation.glaucoma.std$amed = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.exfoliation.glaucoma.std$ahei = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.exfoliation.glaucoma.std$dash = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.exfoliation.glaucoma.std$opdi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.exfoliation.glaucoma.std$hpdi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.exfoliation.glaucoma.std$updi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.exfoliation.glaucoma.std$edip = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.exfoliation.glaucoma.std$edih = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.parkinsons.std
hpfs.parkinsons.std$amed = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.parkinsons.std$ahei = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.parkinsons.std$dash = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.parkinsons.std$opdi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.parkinsons.std$hpdi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.parkinsons.std$updi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.parkinsons.std$edip = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.parkinsons.std$edih = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.poag.std
hpfs.poag.std$amed = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.poag.std$ahei = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.poag.std$dash = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.poag.std$opdi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.poag.std$hpdi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.poag.std$updi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.poag.std$edip = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.poag.std$edih = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#merge nhs1, nhs2, and hpfs sub-studies (n=13,258)
sig <- c("amed","ahei","dash","opdi","hpdi","updi","edip","edih")

nhs1_met6 <- rbind(nhs1.stroke.std[sig], nhs1.als.std[sig], nhs1.breast.std[sig], nhs1.colon.std[sig], nhs1.diabetes.std[sig],
                   nhs1.exfoliation.glaucoma.std[sig], nhs1.ibd.std[sig], nhs1.ovarian.std[sig], nhs1.parkinsons.std[sig],
                   nhs1.poag.std[sig], nhs1.racial.diff.std[sig],nhs1.rheumatoid.std[sig]) #12 sub-studies, 8134

nhs2_met6 <- rbind(nhs2.breast.std[sig], nhs2.ibd.std[sig], nhs2.diabetes.std[sig], nhs2.ovarian.std[sig], nhs2.poag.std[sig],
                   nhs2.rheumatoid.std[sig], nhs2.stress.std[sig], nhs2.stroke.std[sig]) #8 sub-studies, 3473

hpfs_met6 <- rbind(hpfs.als.std[sig], hpfs.colon.std[sig], hpfs.exfoliation.glaucoma.std[sig], hpfs.parkinsons.std[sig],
                   hpfs.poag.std[sig]) #5 sub-studies, 1651

nhs1_met6$id <- rownames(nhs1_met6); nhs2_met6$id <- rownames(nhs2_met6); hpfs_met6$id <- rownames(hpfs_met6)

#overlap between two
amed <- data.frame(amed_4[,2]); rownames(amed) <- amed_4$`AMED (n=48)`
ahei <- data.frame(ahei_4[,2]); rownames(ahei) <- ahei_4$`AHEI-2010 (n=61)`
dash <- data.frame(dash_4[,2]); rownames(dash) <- dash_4$`DASH (n=55)`
opdi <- data.frame(pdi_4[,2]); rownames(opdi) <- pdi_4$`oPDI (n=36)`
hpdi <- data.frame(hpdi_4[,2]); rownames(hpdi) <- hpdi_4$`hPDI (n=39)`
updi <- data.frame(updi_4[,2]); rownames(updi) <- updi_4$`uPDI (n=43)`
edip <- data.frame(edip_2[,2]); rownames(edip) <- edip_2$`EDIP (n=66)`
edih <- data.frame(edih_2[,2]); rownames(edih) <- edih_2$`EDIH (n=37)`

#calculate metabolic signature in each sub-study
#nhs1.stroke.std
nhs1.stroke.std$amed = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.stroke.std$ahei = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.stroke.std$dash = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.stroke.std$opdi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.stroke.std$hpdi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.stroke.std$updi = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.stroke.std$edip = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.stroke.std$edih = apply(mapply(`*`, nhs1.stroke.std[,which(colnames(nhs1.stroke.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.stroke.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.als.std
nhs1.als.std$amed = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.als.std$ahei = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.als.std$dash = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.als.std$opdi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.als.std$hpdi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.als.std$updi = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.als.std$edip = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.als.std$edih = apply(mapply(`*`, nhs1.als.std[,which(colnames(nhs1.als.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.als.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.breast.std
nhs1.breast.std$amed = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.breast.std$ahei = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.breast.std$dash = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.breast.std$opdi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.breast.std$hpdi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.breast.std$updi = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.breast.std$edip = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.breast.std$edih = apply(mapply(`*`, nhs1.breast.std[,which(colnames(nhs1.breast.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.breast.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.colon.std
nhs1.colon.std$amed = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.colon.std$ahei = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.colon.std$dash = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.colon.std$opdi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.colon.std$hpdi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.colon.std$updi = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.colon.std$edip = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.colon.std$edih = apply(mapply(`*`, nhs1.colon.std[,which(colnames(nhs1.colon.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.colon.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.diabetes.std
nhs1.diabetes.std$amed = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.diabetes.std$ahei = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.diabetes.std$dash = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.diabetes.std$opdi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.diabetes.std$hpdi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.diabetes.std$updi = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.diabetes.std$edip = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.diabetes.std$edih = apply(mapply(`*`, nhs1.diabetes.std[,which(colnames(nhs1.diabetes.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.diabetes.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.exfoliation.glaucoma.std
nhs1.exfoliation.glaucoma.std$amed = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.exfoliation.glaucoma.std$ahei = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.exfoliation.glaucoma.std$dash = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.exfoliation.glaucoma.std$opdi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.exfoliation.glaucoma.std$hpdi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.exfoliation.glaucoma.std$updi = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.exfoliation.glaucoma.std$edip = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.exfoliation.glaucoma.std$edih = apply(mapply(`*`, nhs1.exfoliation.glaucoma.std[,which(colnames(nhs1.exfoliation.glaucoma.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.ibd.std
nhs1.ibd.std$amed = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.ibd.std$ahei = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.ibd.std$dash = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.ibd.std$opdi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.ibd.std$hpdi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.ibd.std$updi = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.ibd.std$edip = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.ibd.std$edih = apply(mapply(`*`, nhs1.ibd.std[,which(colnames(nhs1.ibd.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.ibd.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.ovarian.std
nhs1.ovarian.std$amed = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.ovarian.std$ahei = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.ovarian.std$dash = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.ovarian.std$opdi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.ovarian.std$hpdi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.ovarian.std$updi = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.ovarian.std$edip = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.ovarian.std$edih = apply(mapply(`*`, nhs1.ovarian.std[,which(colnames(nhs1.ovarian.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.ovarian.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.parkinsons.std
nhs1.parkinsons.std$amed = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.parkinsons.std$ahei = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.parkinsons.std$dash = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.parkinsons.std$opdi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.parkinsons.std$hpdi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.parkinsons.std$updi = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.parkinsons.std$edip = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.parkinsons.std$edih = apply(mapply(`*`, nhs1.parkinsons.std[,which(colnames(nhs1.parkinsons.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.parkinsons.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.poag.std
nhs1.poag.std$amed = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.poag.std$ahei = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.poag.std$dash = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.poag.std$opdi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.poag.std$hpdi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.poag.std$updi = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.poag.std$edip = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.poag.std$edih = apply(mapply(`*`, nhs1.poag.std[,which(colnames(nhs1.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.racial.diff.std
nhs1.racial.diff.std$amed = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.racial.diff.std$ahei = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.racial.diff.std$dash = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.racial.diff.std$opdi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.racial.diff.std$hpdi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.racial.diff.std$updi = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.racial.diff.std$edip = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.racial.diff.std$edih = apply(mapply(`*`, nhs1.racial.diff.std[,which(colnames(nhs1.racial.diff.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.racial.diff.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs1.rheumatoid.std
nhs1.rheumatoid.std$amed = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs1.rheumatoid.std$ahei = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs1.rheumatoid.std$dash = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs1.rheumatoid.std$opdi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs1.rheumatoid.std$hpdi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs1.rheumatoid.std$updi = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs1.rheumatoid.std$edip = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs1.rheumatoid.std$edih = apply(mapply(`*`, nhs1.rheumatoid.std[,which(colnames(nhs1.rheumatoid.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs1.rheumatoid.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.breast.std
nhs2.breast.std$amed = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.breast.std$ahei = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.breast.std$dash = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.breast.std$opdi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.breast.std$hpdi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.breast.std$updi = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.breast.std$edip = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.breast.std$edih = apply(mapply(`*`, nhs2.breast.std[,which(colnames(nhs2.breast.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.breast.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.ibd.std
nhs2.ibd.std$amed = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.ibd.std$ahei = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.ibd.std$dash = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.ibd.std$opdi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.ibd.std$hpdi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.ibd.std$updi = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.ibd.std$edip = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.ibd.std$edih = apply(mapply(`*`, nhs2.ibd.std[,which(colnames(nhs2.ibd.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.ibd.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.diabetes.std
nhs2.diabetes.std$amed = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.diabetes.std$ahei = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.diabetes.std$dash = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.diabetes.std$opdi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.diabetes.std$hpdi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.diabetes.std$updi = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.diabetes.std$edip = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.diabetes.std$edih = apply(mapply(`*`, nhs2.diabetes.std[,which(colnames(nhs2.diabetes.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.diabetes.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.ovarian.std
nhs2.ovarian.std$amed = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.ovarian.std$ahei = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.ovarian.std$dash = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.ovarian.std$opdi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.ovarian.std$hpdi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.ovarian.std$updi = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.ovarian.std$edip = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.ovarian.std$edih = apply(mapply(`*`, nhs2.ovarian.std[,which(colnames(nhs2.ovarian.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.ovarian.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.poag.std
nhs2.poag.std$amed = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.poag.std$ahei = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.poag.std$dash = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.poag.std$opdi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.poag.std$hpdi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.poag.std$updi = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.poag.std$edip = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.poag.std$edih = apply(mapply(`*`, nhs2.poag.std[,which(colnames(nhs2.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.rheumatoid.std
nhs2.rheumatoid.std$amed = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.rheumatoid.std$ahei = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.rheumatoid.std$dash = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.rheumatoid.std$opdi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.rheumatoid.std$hpdi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.rheumatoid.std$updi = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.rheumatoid.std$edip = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.rheumatoid.std$edih = apply(mapply(`*`, nhs2.rheumatoid.std[,which(colnames(nhs2.rheumatoid.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.rheumatoid.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.stress.std
nhs2.stress.std$amed = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.stress.std$ahei = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.stress.std$dash = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.stress.std$opdi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.stress.std$hpdi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.stress.std$updi = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.stress.std$edip = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.stress.std$edih = apply(mapply(`*`, nhs2.stress.std[,which(colnames(nhs2.stress.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.stress.std)),]))), 1, sum)+as.numeric(edih[1,1])

#nhs2.stroke.std
nhs2.stroke.std$amed = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(amed[1,1])
nhs2.stroke.std$ahei = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(ahei[1,1])
nhs2.stroke.std$dash = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(dash[1,1])
nhs2.stroke.std$opdi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(opdi[1,1])
nhs2.stroke.std$hpdi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
nhs2.stroke.std$updi = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(updi[1,1])
nhs2.stroke.std$edip = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(edip[1,1])
nhs2.stroke.std$edih = apply(mapply(`*`, nhs2.stroke.std[,which(colnames(nhs2.stroke.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(nhs2.stroke.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.als.std
hpfs.als.std$amed = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.als.std$ahei = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.als.std$dash = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.als.std$opdi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.als.std$hpdi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.als.std$updi = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.als.std$edip = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.als.std$edih = apply(mapply(`*`, hpfs.als.std[,which(colnames(hpfs.als.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.als.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.colon.std
hpfs.colon.std$amed = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.colon.std$ahei = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.colon.std$dash = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.colon.std$opdi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.colon.std$hpdi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.colon.std$updi = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.colon.std$edip = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.colon.std$edih = apply(mapply(`*`, hpfs.colon.std[,which(colnames(hpfs.colon.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.colon.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.exfoliation.glaucoma.std
hpfs.exfoliation.glaucoma.std$amed = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.exfoliation.glaucoma.std$ahei = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.exfoliation.glaucoma.std$dash = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.exfoliation.glaucoma.std$opdi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.exfoliation.glaucoma.std$hpdi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.exfoliation.glaucoma.std$updi = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.exfoliation.glaucoma.std$edip = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.exfoliation.glaucoma.std$edih = apply(mapply(`*`, hpfs.exfoliation.glaucoma.std[,which(colnames(hpfs.exfoliation.glaucoma.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.exfoliation.glaucoma.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.parkinsons.std
hpfs.parkinsons.std$amed = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.parkinsons.std$ahei = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.parkinsons.std$dash = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.parkinsons.std$opdi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.parkinsons.std$hpdi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.parkinsons.std$updi = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.parkinsons.std$edip = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.parkinsons.std$edih = apply(mapply(`*`, hpfs.parkinsons.std[,which(colnames(hpfs.parkinsons.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.parkinsons.std)),]))), 1, sum)+as.numeric(edih[1,1])

#hpfs.poag.std
hpfs.poag.std$amed = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(amed))], as.numeric(t(amed[which(rownames(amed) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(amed[1,1])
hpfs.poag.std$ahei = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(ahei))], as.numeric(t(ahei[which(rownames(ahei) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(ahei[1,1])
hpfs.poag.std$dash = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(dash))], as.numeric(t(dash[which(rownames(dash) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(dash[1,1])
hpfs.poag.std$opdi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(opdi))], as.numeric(t(opdi[which(rownames(opdi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(opdi[1,1])
hpfs.poag.std$hpdi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(hpdi))], as.numeric(t(hpdi[which(rownames(hpdi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(hpdi[1,1])
hpfs.poag.std$updi = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(updi))], as.numeric(t(updi[which(rownames(updi) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(updi[1,1])
hpfs.poag.std$edip = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(edip))], as.numeric(t(edip[which(rownames(edip) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(edip[1,1])
hpfs.poag.std$edih = apply(mapply(`*`, hpfs.poag.std[,which(colnames(hpfs.poag.std) %in% rownames(edih))], as.numeric(t(edih[which(rownames(edih) %in% colnames(hpfs.poag.std)),]))), 1, sum)+as.numeric(edih[1,1])

#merge nhs1, nhs2, and hpfs sub-studies (n=13,258)
sig <- c("amed","ahei","dash","opdi","hpdi","updi","edip","edih")

nhs1_met7 <- rbind(nhs1.stroke.std[sig], nhs1.als.std[sig], nhs1.breast.std[sig], nhs1.colon.std[sig], nhs1.diabetes.std[sig],
                   nhs1.exfoliation.glaucoma.std[sig], nhs1.ibd.std[sig], nhs1.ovarian.std[sig], nhs1.parkinsons.std[sig],
                   nhs1.poag.std[sig], nhs1.racial.diff.std[sig],nhs1.rheumatoid.std[sig]) #12 sub-studies, 8134

nhs2_met7 <- rbind(nhs2.breast.std[sig], nhs2.ibd.std[sig], nhs2.diabetes.std[sig], nhs2.ovarian.std[sig], nhs2.poag.std[sig],
                   nhs2.rheumatoid.std[sig], nhs2.stress.std[sig], nhs2.stroke.std[sig]) #8 sub-studies, 3473

hpfs_met7 <- rbind(hpfs.als.std[sig], hpfs.colon.std[sig], hpfs.exfoliation.glaucoma.std[sig], hpfs.parkinsons.std[sig],
                   hpfs.poag.std[sig]) #5 sub-studies, 1651

nhs1_met7$id <- rownames(nhs1_met7); nhs2_met7$id <- rownames(nhs2_met7); hpfs_met7$id <- rownames(hpfs_met7)

#rename all signatures
names(nhs1_met3) <- names(nhs2_met3) <- names(hpfs_met3) <- c("amed3","ahei3","dash3","pdi3","hpdi3","updi3","edip3","edih3","id")
names(nhs1_met4) <- names(nhs2_met4) <- names(hpfs_met4) <- c("amed4","ahei4","dash4","pdi4","hpdi4","updi4","edip4","edih4","id")
names(nhs1_met5) <- names(nhs2_met5) <- names(hpfs_met5) <- c("amed5","ahei5","dash5","pdi5","hpdi5","updi5","edip5","edih5","id")
names(nhs1_met6) <- names(nhs2_met6) <- names(hpfs_met6) <- c("amed6","ahei6","dash6","pdi6","hpdi6","updi6","edip6","edih6","id")
names(nhs1_met7) <- names(nhs2_met7) <- names(hpfs_met7) <- c("amed7","ahei7","dash7","pdi7","hpdi7","updi7","edip7","edih7","id")

#combine all signatures
nhs1_met <- merge(nhs1_met3,nhs1_met4, by="id") %>% merge(nhs1_met5, by="id") %>% merge(nhs1_met6[,-c(7:8)], by="id") %>% merge(nhs1_met7[,-c(7:8)], by="id")
nhs2_met <- merge(nhs2_met3,nhs2_met4, by="id") %>% merge(nhs2_met5, by="id") %>% merge(nhs2_met6[,-c(7:8)], by="id") %>% merge(nhs2_met7[,-c(7:8)], by="id")
hpfs_met <- merge(hpfs_met3,hpfs_met4, by="id") %>% merge(hpfs_met5, by="id") %>% merge(hpfs_met6[,-c(7:8)], by="id") %>% merge(hpfs_met7[,-c(7:8)], by="id")

#combine all data
all <- rbind(nhs1_met,nhs2_met,hpfs_met) #13258
all$id <- substr(all$id,1,6)

#remove duplicates
all2 <- all[!duplicated(all$id),] #12453

#--------------------------------------------------------------------------------------------
#
#           Step3: mediation analysis to explore mediated proportion of shared and distinct metabolic signatures in diet-T2D associations
#
#--------------------------------------------------------------------------------------------
#load phenotype data
load("NHSHPFS_Final_T2D.RData")

all.int.use = t2d4[!duplicated(t2d4$id),]

all.int.use$caloravn = all.int.use$caloravn.x

#merge with newly calculated signature-mediator
all.int.use <- merge(all.int.use,all2,by="id")

#check covariates
covar = c("ageyr","studycaco","fast","aspirinuse","smoking","actcat","fhxdb","phxchol","phxhbp","caloravn","BMIcata")

#impute missing continuous covariate using median value
medfunc <- function(x){
  missplace <- which(is.na(x))
  x[missplace] <- median(x, na.rm = T)/2
  x
}

all.int.use[,c("ageyr","caloravn")] <- apply(all.int.use[,c("ageyr","caloravn")],2,medfunc)

#define mediators
amed_l <- c("amed2","amed3","amed4","amed5","amed6","amed7")
ahei_l <- c("ahei2","ahei3","ahei4","ahei5","ahei6","ahei7")
dash_l <- c("dash2","dash3","dash4","dash5","dash6","dash7")
pdi_l <- c("pdi2","pdi3","pdi4","pdi5","pdi6","pdi7")
hpdi_l <- c("hpdi2","hpdi3","hpdi4","hpdi5","hpdi6","hpdi7")
updi_l <- c("updi2","updi3","updi4","updi5","updi6","updi7")
edip_l <- c("edip2","edip3","edip4","edip5")
edih_l <- c("edih2","edih3","edih4","edih5")

#define exposure
share_exposure_list = c("amed_av","ahei_av","dash_av","pdi_av","hpdi_av","updi_av","edip_av","edih_av")

#define covariates
adj_covar = covar

###amed
#define mediator-cross-platform signature
metab_names <- amed_l

#define the controlled values of mediators
metab_dataset <- all.int.use %>% dplyr::select(c(metab_names,"id"))

controlled_med<- metab_dataset%>%
  dplyr::select(metab_names)%>%data.table::melt()%>%
  group_by(variable)%>%
  summarize(median=median(value,na.rm = T))

controlled_med<-as.vector(controlled_med$median)

#run model
nb <- 200

singlemed_test_result<-list()  

exp_x = "amed_av"

for(i in 1:length(metab_names)){
  
  dati =all.int.use[,c(exp_x,metab_names[i],"ptime_diabetes","diabetes",adj_covar)] #note: mediator,exposure,outcome,event can not be NA
  dati = dati[complete.cases(dati),]
  
  adj_covar_use = c("ageyr","caloravn")
  
  #add smoking
  if(length(unique(dati$current.smoking))>1) {
    adj_covar_use = c(adj_covar_use,"smoking")
  }  
  
  #add fast
  if(length(unique(dati$fast))>1) {
    adj_covar_use = c(adj_covar_use,"fast")
  }  
  
  #add studycaco
  if(length(unique(dati$caco))>1) {
    adj_covar_use = c(adj_covar_use,"studycaco")
  }
  
  #add aspirinuse
  if(length(unique(dati$antihluse))>1) {
    adj_covar_use = c(adj_covar_use,"aspirinuse")
  }
  
  #add fhxdb
  if(length(unique(dati$fhxdb))>1) {
    adj_covar_use = c(adj_covar_use,"fhxdb")
  }
  
  #add phxhbp
  if(length(unique(dati$phxhbp))>1) {
    adj_covar_use = c(adj_covar_use,"phxhbp")
  }
  
  #add phxchol
  if(length(unique(dati$phxchol))>1) {
    adj_covar_use = c(adj_covar_use,"phxchol")
  }
  
  #add actcat
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"actcat")
  }
  
  #add BMIcata
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"BMIcata")
  }  
  
  med_list <- as.vector(metab_names[i]) 
  med_controll_value_list <- as.list(controlled_med[1])      #this is a list with single value
  mreg_list <- as.list(rep("linear",length(metab_names[i]))) #this for models in mediator analysis 
  
  model_cmest<- cmest(data = dati, 
                      model = "rb", 
                      outcome = "ptime_diabetes",          
                      event = "diabetes",                   
                      exposure = exp_x,                   
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
                      nboot=nb,
                      yvar=1) 
  
  #print(paste0(a,"-",i)) #to monitoring the process of cmest model
  
  save_singlemed_result<-as.data.frame(as.matrix(summary(model_cmest)$summarydf))
  save_singlemed_result<-save_singlemed_result%>%rownames_to_column()
  
  save_singlemed_result$mediators<-paste0(metab_names[i],collapse = "+")
  singlemed_test_result[[i]]<-save_singlemed_result
}

singlemed_test_result_save <- do.call(rbind.data.frame,singlemed_test_result)

singlemed_test_result_save <- singlemed_test_result_save %>% group_by(rowname) %>% mutate(fdr=p.adjust(`P.val`,"BH"))

singlemed_test_result_save$exposure<-exp_x
singlemed_test_result_save$type<-"Single_mediator"

#save results
amed <- singlemed_test_result_save

###ahei
#define mediator-cross-platform signature
metab_names <- ahei_l

#define the controlled values of mediators
metab_dataset <- all.int.use %>% dplyr::select(metab_names)

controlled_med<- metab_dataset%>%
  dplyr::select(metab_names)%>%data.table::melt()%>%
  group_by(variable)%>%
  summarize(median=median(value,na.rm = T))

controlled_med<-as.vector(controlled_med$median)

#run model
nb <- 200

singlemed_test_result<-list()  

exp_x = "ahei2010_av"

for(i in 1:length(metab_names)){
  
  dati =all.int.use[,c(exp_x,metab_names[i],"ptime_diabetes","diabetes",adj_covar)] #note: mediator,exposure,outcome,event can not be NA
  dati = dati[complete.cases(dati),]
  
  adj_covar_use = c("ageyr","caloravn")
  
  #add smoking
  if(length(unique(dati$current.smoking))>1) {
    adj_covar_use = c(adj_covar_use,"smoking")
  }  
  
  #add fast
  if(length(unique(dati$fast))>1) {
    adj_covar_use = c(adj_covar_use,"fast")
  }  
  
  #add studycaco
  if(length(unique(dati$caco))>1) {
    adj_covar_use = c(adj_covar_use,"studycaco")
  }
  
  #add aspirinuse
  if(length(unique(dati$antihluse))>1) {
    adj_covar_use = c(adj_covar_use,"aspirinuse")
  }
  
  #add fhxdb
  if(length(unique(dati$fhxdb))>1) {
    adj_covar_use = c(adj_covar_use,"fhxdb")
  }
  
  #add phxhbp
  if(length(unique(dati$phxhbp))>1) {
    adj_covar_use = c(adj_covar_use,"phxhbp")
  }
  
  #add phxchol
  if(length(unique(dati$phxchol))>1) {
    adj_covar_use = c(adj_covar_use,"phxchol")
  }
  
  #add actcat
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"actcat")
  }
  
  #add BMIcata
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"BMIcata")
  }  
  
  med_list <- as.vector(metab_names[i]) 
  med_controll_value_list <- as.list(controlled_med[1])      #this is a list with single value
  mreg_list <- as.list(rep("linear",length(metab_names[i]))) #this for models in mediator analysis 
  
  model_cmest<- cmest(data = dati, 
                      model = "rb", 
                      outcome = "ptime_diabetes",          
                      event = "diabetes",                   
                      exposure = exp_x,                   
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
                      nboot=nb,
                      yvar=1) 
  
  #print(paste0(a,"-",i)) #to monitoring the process of cmest model
  
  save_singlemed_result<-as.data.frame(as.matrix(summary(model_cmest)$summarydf))
  save_singlemed_result<-save_singlemed_result%>%rownames_to_column()
  
  save_singlemed_result$mediators<-paste0(metab_names[i],collapse = "+")
  singlemed_test_result[[i]]<-save_singlemed_result
}

singlemed_test_result_save <- do.call(rbind.data.frame,singlemed_test_result)

singlemed_test_result_save <- singlemed_test_result_save %>% group_by(rowname) %>% mutate(fdr=p.adjust(`P.val`,"BH"))

singlemed_test_result_save$exposure<-exp_x
singlemed_test_result_save$type<-"Single_mediator"

#save results
ahei <- singlemed_test_result_save

###dash
#define mediator-cross-platform signature
metab_names <- dash_l

#define the controlled values of mediators
metab_dataset <- all.int.use %>% dplyr::select(metab_names)

controlled_med<- metab_dataset%>%
  dplyr::select(metab_names)%>%data.table::melt()%>%
  group_by(variable)%>%
  summarize(median=median(value,na.rm = T))

controlled_med<-as.vector(controlled_med$median)

#run model
nb <- 200

singlemed_test_result<-list()  

exp_x = "dashav"

for(i in 1:length(metab_names)){
  
  dati =all.int.use[,c(exp_x,metab_names[i],"ptime_diabetes","diabetes",adj_covar)] #note: mediator,exposure,outcome,event can not be NA
  dati = dati[complete.cases(dati),]
  
  adj_covar_use = c("ageyr","caloravn")
  
  #add smoking
  if(length(unique(dati$current.smoking))>1) {
    adj_covar_use = c(adj_covar_use,"smoking")
  }  
  
  #add fast
  if(length(unique(dati$fast))>1) {
    adj_covar_use = c(adj_covar_use,"fast")
  }  
  
  #add studycaco
  if(length(unique(dati$caco))>1) {
    adj_covar_use = c(adj_covar_use,"studycaco")
  }
  
  #add aspirinuse
  if(length(unique(dati$antihluse))>1) {
    adj_covar_use = c(adj_covar_use,"aspirinuse")
  }
  
  #add fhxdb
  if(length(unique(dati$fhxdb))>1) {
    adj_covar_use = c(adj_covar_use,"fhxdb")
  }
  
  #add phxhbp
  if(length(unique(dati$phxhbp))>1) {
    adj_covar_use = c(adj_covar_use,"phxhbp")
  }
  
  #add phxchol
  if(length(unique(dati$phxchol))>1) {
    adj_covar_use = c(adj_covar_use,"phxchol")
  }
  
  #add actcat
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"actcat")
  }
  
  #add BMIcata
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"BMIcata")
  }  
  
  med_list <- as.vector(metab_names[i]) 
  med_controll_value_list <- as.list(controlled_med[1])      #this is a list with single value
  mreg_list <- as.list(rep("linear",length(metab_names[i]))) #this for models in mediator analysis 
  
  model_cmest<- cmest(data = dati, 
                      model = "rb", 
                      outcome = "ptime_diabetes",          
                      event = "diabetes",                   
                      exposure = exp_x,                   
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
                      nboot=nb,
                      yvar=1) 
  
  #print(paste0(a,"-",i)) #to monitoring the process of cmest model
  
  save_singlemed_result<-as.data.frame(as.matrix(summary(model_cmest)$summarydf))
  save_singlemed_result<-save_singlemed_result%>%rownames_to_column()
  
  save_singlemed_result$mediators<-paste0(metab_names[i],collapse = "+")
  singlemed_test_result[[i]]<-save_singlemed_result
}

singlemed_test_result_save <- do.call(rbind.data.frame,singlemed_test_result)

singlemed_test_result_save <- singlemed_test_result_save %>% group_by(rowname) %>% mutate(fdr=p.adjust(`P.val`,"BH"))

singlemed_test_result_save$exposure<-exp_x
singlemed_test_result_save$type<-"Single_mediator"

#save results
dash <- singlemed_test_result_save

###pdi
#define mediator-cross-platform signature
metab_names <- pdi_l

#define the controlled values of mediators
metab_dataset <- all.int.use %>% dplyr::select(metab_names)

controlled_med<- metab_dataset%>%
  dplyr::select(metab_names)%>%data.table::melt()%>%
  group_by(variable)%>%
  summarize(median=median(value,na.rm = T))

controlled_med<-as.vector(controlled_med$median)

#run model
nb <- 200

singlemed_test_result<-list()  

exp_x = "pdi_av"

for(i in 1:length(metab_names)){
  
  dati =all.int.use[,c(exp_x,metab_names[i],"ptime_diabetes","diabetes",adj_covar)] #note: mediator,exposure,outcome,event can not be NA
  dati = dati[complete.cases(dati),]
  
  adj_covar_use = c("ageyr","caloravn")
  
  #add smoking
  if(length(unique(dati$current.smoking))>1) {
    adj_covar_use = c(adj_covar_use,"smoking")
  }  
  
  #add fast
  if(length(unique(dati$fast))>1) {
    adj_covar_use = c(adj_covar_use,"fast")
  }  
  
  #add studycaco
  if(length(unique(dati$caco))>1) {
    adj_covar_use = c(adj_covar_use,"studycaco")
  }
  
  #add aspirinuse
  if(length(unique(dati$antihluse))>1) {
    adj_covar_use = c(adj_covar_use,"aspirinuse")
  }
  
  #add fhxdb
  if(length(unique(dati$fhxdb))>1) {
    adj_covar_use = c(adj_covar_use,"fhxdb")
  }
  
  #add phxhbp
  if(length(unique(dati$phxhbp))>1) {
    adj_covar_use = c(adj_covar_use,"phxhbp")
  }
  
  #add phxchol
  if(length(unique(dati$phxchol))>1) {
    adj_covar_use = c(adj_covar_use,"phxchol")
  }
  
  #add actcat
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"actcat")
  }
  
  #add BMIcata
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"BMIcata")
  }  
  
  med_list <- as.vector(metab_names[i]) 
  med_controll_value_list <- as.list(controlled_med[1])      #this is a list with single value
  mreg_list <- as.list(rep("linear",length(metab_names[i]))) #this for models in mediator analysis 
  
  model_cmest<- cmest(data = dati, 
                      model = "rb", 
                      outcome = "ptime_diabetes",          
                      event = "diabetes",                   
                      exposure = exp_x,                   
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
                      nboot=nb,
                      yvar=1) 
  
  #print(paste0(a,"-",i)) #to monitoring the process of cmest model
  
  save_singlemed_result<-as.data.frame(as.matrix(summary(model_cmest)$summarydf))
  save_singlemed_result<-save_singlemed_result%>%rownames_to_column()
  
  save_singlemed_result$mediators<-paste0(metab_names[i],collapse = "+")
  singlemed_test_result[[i]]<-save_singlemed_result
}

singlemed_test_result_save <- do.call(rbind.data.frame,singlemed_test_result)

singlemed_test_result_save <- singlemed_test_result_save %>% group_by(rowname) %>% mutate(fdr=p.adjust(`P.val`,"BH"))

singlemed_test_result_save$exposure<-exp_x
singlemed_test_result_save$type<-"Single_mediator"

#save results
pdi <- singlemed_test_result_save

###hpdi
#define mediator-cross-platform signature
metab_names <- hpdi_l

#define the controlled values of mediators
metab_dataset <- all.int.use %>% dplyr::select(metab_names)

controlled_med<- metab_dataset%>%
  dplyr::select(metab_names)%>%data.table::melt()%>%
  group_by(variable)%>%
  summarize(median=median(value,na.rm = T))

controlled_med<-as.vector(controlled_med$median)

#run model
nb <- 200

singlemed_test_result<-list()  

exp_x = "hpdi_av"

for(i in 1:length(metab_names)){
  
  dati =all.int.use[,c(exp_x,metab_names[i],"ptime_diabetes","diabetes",adj_covar)] #note: mediator,exposure,outcome,event can not be NA
  dati = dati[complete.cases(dati),]
  
  adj_covar_use = c("ageyr","caloravn")
  
  #add smoking
  if(length(unique(dati$current.smoking))>1) {
    adj_covar_use = c(adj_covar_use,"smoking")
  }  
  
  #add fast
  if(length(unique(dati$fast))>1) {
    adj_covar_use = c(adj_covar_use,"fast")
  }  
  
  #add studycaco
  if(length(unique(dati$caco))>1) {
    adj_covar_use = c(adj_covar_use,"studycaco")
  }
  
  #add aspirinuse
  if(length(unique(dati$antihluse))>1) {
    adj_covar_use = c(adj_covar_use,"aspirinuse")
  }
  
  #add fhxdb
  if(length(unique(dati$fhxdb))>1) {
    adj_covar_use = c(adj_covar_use,"fhxdb")
  }
  
  #add phxhbp
  if(length(unique(dati$phxhbp))>1) {
    adj_covar_use = c(adj_covar_use,"phxhbp")
  }
  
  #add phxchol
  if(length(unique(dati$phxchol))>1) {
    adj_covar_use = c(adj_covar_use,"phxchol")
  }
  
  #add actcat
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"actcat")
  }
  
  #add BMIcata
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"BMIcata")
  }  
  
  med_list <- as.vector(metab_names[i]) 
  med_controll_value_list <- as.list(controlled_med[1])      #this is a list with single value
  mreg_list <- as.list(rep("linear",length(metab_names[i]))) #this for models in mediator analysis 
  
  model_cmest<- cmest(data = dati, 
                      model = "rb", 
                      outcome = "ptime_diabetes",          
                      event = "diabetes",                   
                      exposure = exp_x,                   
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
                      nboot=nb,
                      yvar=1) 
  
  #print(paste0(a,"-",i)) #to monitoring the process of cmest model
  
  save_singlemed_result<-as.data.frame(as.matrix(summary(model_cmest)$summarydf))
  save_singlemed_result<-save_singlemed_result%>%rownames_to_column()
  
  save_singlemed_result$mediators<-paste0(metab_names[i],collapse = "+")
  singlemed_test_result[[i]]<-save_singlemed_result
}

singlemed_test_result_save <- do.call(rbind.data.frame,singlemed_test_result)

singlemed_test_result_save <- singlemed_test_result_save %>% group_by(rowname) %>% mutate(fdr=p.adjust(`P.val`,"BH"))

singlemed_test_result_save$exposure<-exp_x
singlemed_test_result_save$type<-"Single_mediator"

#save results
hpdi <- singlemed_test_result_save

###updi
#define mediator-cross-platform signature
metab_names <- updi_l

#define the controlled values of mediators
metab_dataset <- all.int.use %>% dplyr::select(metab_names)

controlled_med<- metab_dataset%>%
  dplyr::select(metab_names)%>%data.table::melt()%>%
  group_by(variable)%>%
  summarize(median=median(value,na.rm = T))

controlled_med<-as.vector(controlled_med$median)

#run model
nb <- 200

singlemed_test_result<-list()  

exp_x = "updi_av"

for(i in 1:length(metab_names)){
  
  dati =all.int.use[,c(exp_x,metab_names[i],"ptime_diabetes","diabetes",adj_covar)] #note: mediator,exposure,outcome,event can not be NA
  dati = dati[complete.cases(dati),]
  
  adj_covar_use = c("ageyr","caloravn")
  
  #add smoking
  if(length(unique(dati$current.smoking))>1) {
    adj_covar_use = c(adj_covar_use,"smoking")
  }  
  
  #add fast
  if(length(unique(dati$fast))>1) {
    adj_covar_use = c(adj_covar_use,"fast")
  }  
  
  #add studycaco
  if(length(unique(dati$caco))>1) {
    adj_covar_use = c(adj_covar_use,"studycaco")
  }
  
  #add aspirinuse
  if(length(unique(dati$antihluse))>1) {
    adj_covar_use = c(adj_covar_use,"aspirinuse")
  }
  
  #add fhxdb
  if(length(unique(dati$fhxdb))>1) {
    adj_covar_use = c(adj_covar_use,"fhxdb")
  }
  
  #add phxhbp
  if(length(unique(dati$phxhbp))>1) {
    adj_covar_use = c(adj_covar_use,"phxhbp")
  }
  
  #add phxchol
  if(length(unique(dati$phxchol))>1) {
    adj_covar_use = c(adj_covar_use,"phxchol")
  }
  
  #add actcat
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"actcat")
  }
  
  #add BMIcata
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"BMIcata")
  }  
  
  med_list <- as.vector(metab_names[i]) 
  med_controll_value_list <- as.list(controlled_med[1])      #this is a list with single value
  mreg_list <- as.list(rep("linear",length(metab_names[i]))) #this for models in mediator analysis 
  
  model_cmest<- cmest(data = dati, 
                      model = "rb", 
                      outcome = "ptime_diabetes",          
                      event = "diabetes",                   
                      exposure = exp_x,                   
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
                      nboot=nb,
                      yvar=1) 
  
  #print(paste0(a,"-",i)) #to monitoring the process of cmest model
  
  save_singlemed_result<-as.data.frame(as.matrix(summary(model_cmest)$summarydf))
  save_singlemed_result<-save_singlemed_result%>%rownames_to_column()
  
  save_singlemed_result$mediators<-paste0(metab_names[i],collapse = "+")
  singlemed_test_result[[i]]<-save_singlemed_result
}

singlemed_test_result_save <- do.call(rbind.data.frame,singlemed_test_result)

singlemed_test_result_save <- singlemed_test_result_save %>% group_by(rowname) %>% mutate(fdr=p.adjust(`P.val`,"BH"))

singlemed_test_result_save$exposure<-exp_x
singlemed_test_result_save$type<-"Single_mediator"

#save results
updi <- singlemed_test_result_save

###edip
#define mediator-cross-platform signature
metab_names <- edip_l

#define the controlled values of mediators
metab_dataset <- all.int.use %>% dplyr::select(metab_names)

controlled_med<- metab_dataset%>%
  dplyr::select(metab_names)%>%data.table::melt()%>%
  group_by(variable)%>%
  summarize(median=median(value,na.rm = T))

controlled_med<-as.vector(controlled_med$median)

#run model
nb <- 200

singlemed_test_result<-list()  

exp_x = "edipav"

for(i in 1:length(metab_names)){
  
  dati =all.int.use[,c(exp_x,metab_names[i],"ptime_diabetes","diabetes",adj_covar)] #note: mediator,exposure,outcome,event can not be NA
  dati = dati[complete.cases(dati),]
  
  adj_covar_use = c("ageyr","caloravn")
  
  #add smoking
  if(length(unique(dati$current.smoking))>1) {
    adj_covar_use = c(adj_covar_use,"smoking")
  }  
  
  #add fast
  if(length(unique(dati$fast))>1) {
    adj_covar_use = c(adj_covar_use,"fast")
  }  
  
  #add studycaco
  if(length(unique(dati$caco))>1) {
    adj_covar_use = c(adj_covar_use,"studycaco")
  }
  
  #add aspirinuse
  if(length(unique(dati$antihluse))>1) {
    adj_covar_use = c(adj_covar_use,"aspirinuse")
  }
  
  #add fhxdb
  if(length(unique(dati$fhxdb))>1) {
    adj_covar_use = c(adj_covar_use,"fhxdb")
  }
  
  #add phxhbp
  if(length(unique(dati$phxhbp))>1) {
    adj_covar_use = c(adj_covar_use,"phxhbp")
  }
  
  #add phxchol
  if(length(unique(dati$phxchol))>1) {
    adj_covar_use = c(adj_covar_use,"phxchol")
  }
  
  #add actcat
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"actcat")
  }
  
  #add BMIcata
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"BMIcata")
  }  
  
  med_list <- as.vector(metab_names[i]) 
  med_controll_value_list <- as.list(controlled_med[1])      #this is a list with single value
  mreg_list <- as.list(rep("linear",length(metab_names[i]))) #this for models in mediator analysis 
  
  model_cmest<- cmest(data = dati, 
                      model = "rb", 
                      outcome = "ptime_diabetes",          
                      event = "diabetes",                   
                      exposure = exp_x,                   
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
                      nboot=nb,
                      yvar=1) 
  
  #print(paste0(a,"-",i)) #to monitoring the process of cmest model
  
  save_singlemed_result<-as.data.frame(as.matrix(summary(model_cmest)$summarydf))
  save_singlemed_result<-save_singlemed_result%>%rownames_to_column()
  
  save_singlemed_result$mediators<-paste0(metab_names[i],collapse = "+")
  singlemed_test_result[[i]]<-save_singlemed_result
}

singlemed_test_result_save <- do.call(rbind.data.frame,singlemed_test_result)

singlemed_test_result_save <- singlemed_test_result_save %>% group_by(rowname) %>% mutate(fdr=p.adjust(`P.val`,"BH"))

singlemed_test_result_save$exposure<-exp_x
singlemed_test_result_save$type<-"Single_mediator"

#save results
edip <- singlemed_test_result_save

###edih
#define mediator-cross-platform signature
metab_names <- edih_l

#define the controlled values of mediators
metab_dataset <- all.int.use %>% dplyr::select(metab_names)

controlled_med<- metab_dataset%>%
  dplyr::select(metab_names)%>%data.table::melt()%>%
  group_by(variable)%>%
  summarize(median=median(value,na.rm = T))

controlled_med<-as.vector(controlled_med$median)

#run model
nb <- 200

singlemed_test_result<-list()  

exp_x = "edihav"

for(i in 1:length(metab_names)){
  
  dati =all.int.use[,c(exp_x,metab_names[i],"ptime_diabetes","diabetes",adj_covar)] #note: mediator,exposure,outcome,event can not be NA
  dati = dati[complete.cases(dati),]
  
  adj_covar_use = c("ageyr","caloravn")
  
  #add smoking
  if(length(unique(dati$current.smoking))>1) {
    adj_covar_use = c(adj_covar_use,"smoking")
  }  
  
  #add fast
  if(length(unique(dati$fast))>1) {
    adj_covar_use = c(adj_covar_use,"fast")
  }  
  
  #add studycaco
  if(length(unique(dati$caco))>1) {
    adj_covar_use = c(adj_covar_use,"studycaco")
  }
  
  #add aspirinuse
  if(length(unique(dati$antihluse))>1) {
    adj_covar_use = c(adj_covar_use,"aspirinuse")
  }
  
  #add fhxdb
  if(length(unique(dati$fhxdb))>1) {
    adj_covar_use = c(adj_covar_use,"fhxdb")
  }
  
  #add phxhbp
  if(length(unique(dati$phxhbp))>1) {
    adj_covar_use = c(adj_covar_use,"phxhbp")
  }
  
  #add phxchol
  if(length(unique(dati$phxchol))>1) {
    adj_covar_use = c(adj_covar_use,"phxchol")
  }
  
  #add actcat
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"actcat")
  }
  
  #add BMIcata
  if(exp_x %in% adj_covar){
    adj_covar_use = adj_covar_use
  } else{
    adj_covar_use = c(adj_covar_use,"BMIcata")
  }  
  
  med_list <- as.vector(metab_names[i]) 
  med_controll_value_list <- as.list(controlled_med[1])      #this is a list with single value
  mreg_list <- as.list(rep("linear",length(metab_names[i]))) #this for models in mediator analysis 
  
  model_cmest<- cmest(data = dati, 
                      model = "rb", 
                      outcome = "ptime_diabetes",          
                      event = "diabetes",                   
                      exposure = exp_x,                   
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
                      nboot=nb,
                      yvar=1) 
  
  #print(paste0(a,"-",i)) #to monitoring the process of cmest model
  
  save_singlemed_result<-as.data.frame(as.matrix(summary(model_cmest)$summarydf))
  save_singlemed_result<-save_singlemed_result%>%rownames_to_column()
  
  save_singlemed_result$mediators<-paste0(metab_names[i],collapse = "+")
  singlemed_test_result[[i]]<-save_singlemed_result
}

singlemed_test_result_save <- do.call(rbind.data.frame,singlemed_test_result)

singlemed_test_result_save <- singlemed_test_result_save %>% group_by(rowname) %>% mutate(fdr=p.adjust(`P.val`,"BH"))

singlemed_test_result_save$exposure<-exp_x
singlemed_test_result_save$type<-"Single_mediator"

#save results
edih <- singlemed_test_result_save

#save all resluts
res_all <- rbind(amed,ahei,dash,pdi,hpdi,updi,edip,edih)
write.csv(res_all,file="Figure5.csv")