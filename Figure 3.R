library(data.table)
library(dplyr)
library(readxl)
library(reshape2)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(tidyr)
library(ComplexHeatmap)

#--------------------------------------------------------------------------------------------
#
#           chunk1: estimate association of metabolites included in the signatures with dietary patterns
#
#--------------------------------------------------------------------------------------------
#define inverse normal transformation function
inormal <- function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

#read sample data
load("signature_sample.RData")

#signature information
amed0 <- na.omit(sig_list[,1:2]); amed <- data.frame(amed0[,-1]); rownames(amed) <- amed0[,1]
ahei0 <- na.omit(sig_list[,4:5]); ahei <- data.frame(ahei0[,-1]); rownames(ahei) <- ahei0[,1]
dash0 <- na.omit(sig_list[,7:8]); dash <- data.frame(dash0[,-1]); rownames(dash) <- dash0[,1]
opdi0 <- na.omit(sig_list[,10:11]); opdi <- data.frame(opdi0[,-1]); rownames(opdi) <- opdi0[,1]
hpdi0 <- na.omit(sig_list[,13:14]); hpdi <- data.frame(hpdi0[,-1]); rownames(hpdi) <- hpdi0[,1]
updi0 <- na.omit(sig_list[,16:17]); updi <- data.frame(updi0[,-1]); rownames(updi) <- updi0[,1]
edip0 <- na.omit(sig_list[,19:20]); edip <- data.frame(edip0[,-1]); rownames(edip) <- edip0[,1]
edih0 <- na.omit(sig_list[,22:23]); edih <- data.frame(edih0[,-1]); rownames(edih) <- edih0[,1]

names(amed0)[1] <- names(ahei0)[1] <- names(dash0)[1] <- names(opdi0)[1] <- names(hpdi0)[1] <- names(updi0)[1] <- names(edip0)[1] <- names(edih0)[1] <- "HMDB"

#define variable for use
dp <- c("amed1","ahei1","dash1","opdi1","hpdi1","updi1","edip1","edih1")  #list of dietary patterns in the sample data
var <- unique(c(amed0[-1,]$HMDB,ahei0[-1,]$HMDB,dash0[-1,]$HMDB,opdi0[-1,]$HMDB,hpdi0[-1,]$HMDB,updi0[-1,]$HMDB,edip0[-1,]$HMDB,edih0[-1,]$HMDB)) #list of metabolites included in the signatures

#calculate the association between metabolites included in the signatures and dietary patterns
res = data.frame(Diet=NA,HMDB=NA,Est=NA,SE=NA,P=NA)

x = 0

for (i in var){
  for (j in dp){
    data_use <- as.data.frame(train_sample)
    data_use[c(dp,var)] <- apply(data_use[c(dp,var)],2,inormal)
    data_use$exposure <- data_use[,j]
    data_use$outcome <- data_use[,i]
    fit <- coef(summary(glm(outcome~exposure+ageyr+sex+race+smoke+alco+phxhbp+phxchol+fhxdb+antihluse+act+bmi+energy,data=data_use)))
    x = x+1
    res[x,"Diet"] = j
    res[x,"HMDB"] = i
    res[x,3:5] = fit[2,c(1:2,4:4)]
  }
}

res <- right_join(anno[,c("HMDB","Name")],res,by="HMDB")
res_dp <- res

#rename dietary pattern scores
res_dp$Diet <- gsub("amed1","AMED",res_dp$Diet)
res_dp$Diet <- gsub("ahei1","AHEI",res_dp$Diet)
res_dp$Diet <- gsub("dash1","DASH",res_dp$Diet)
res_dp$Diet <- gsub("opdi1","PDI",res_dp$Diet)
res_dp$Diet <- gsub("hpdi1","hPDI",res_dp$Diet)
res_dp$Diet <- gsub("updi1","uPDI",res_dp$Diet)
res_dp$Diet <- gsub("edip1","EDIP",res_dp$Diet)
res_dp$Diet <- gsub("edih1","EDIH",res_dp$Diet)

#FDR correction across each dietary pattern
res_dp_r1 <- res_dp[which(res_dp$Diet=="AMED"&res_dp$HMDB %in% amed0[-1,]$HMDB),]; res_dp_r1$FDR <- p.adjust(res_dp_r1$P,method = "fdr")
res_dp_r2 <- res_dp[which(res_dp$Diet=="AHEI"&res_dp$HMDB %in% ahei0[-1,]$HMDB),]; res_dp_r2$FDR <- p.adjust(res_dp_r2$P,method = "fdr")
res_dp_r3 <- res_dp[which(res_dp$Diet=="DASH"&res_dp$HMDB %in% dash0[-1,]$HMDB),]; res_dp_r3$FDR <- p.adjust(res_dp_r3$P,method = "fdr")
res_dp_r4 <- res_dp[which(res_dp$Diet=="PDI"&res_dp$HMDB %in% opdi0[-1,]$HMDB),];  res_dp_r4$FDR <- p.adjust(res_dp_r4$P,method = "fdr")
res_dp_r5 <- res_dp[which(res_dp$Diet=="hPDI"&res_dp$HMDB %in% hpdi0[-1,]$HMDB),]; res_dp_r5$FDR <- p.adjust(res_dp_r5$P,method = "fdr")
res_dp_r6 <- res_dp[which(res_dp$Diet=="uPDI"&res_dp$HMDB %in% updi0[-1,]$HMDB),]; res_dp_r6$FDR <- p.adjust(res_dp_r6$P,method = "fdr")
res_dp_r7 <- res_dp[which(res_dp$Diet=="EDIP"&res_dp$HMDB %in% edip0[-1,]$HMDB),]; res_dp_r7$FDR <- p.adjust(res_dp_r7$P,method = "fdr")
res_dp_r8 <- res_dp[which(res_dp$Diet=="EDIH"&res_dp$HMDB %in% edih0[-1,]$HMDB),]; res_dp_r8$FDR <- p.adjust(res_dp_r8$P,method = "fdr")

res_dp_r_use <- rbind(res_dp_r1,res_dp_r2,res_dp_r3,res_dp_r4,res_dp_r5,res_dp_r6,res_dp_r7,res_dp_r8)
res_dp_r_use$sig <- ifelse(res_dp_r_use$FDR <= 0.05 & res_dp_r_use$P < 0.05, "**",ifelse(res_dp_r_use$FDR > 0.05 & res_dp_r_use$P < 0.05, "*", ""))

#transform long data to wide data
coef_dp <- res_dp_r_use[,c(2:4)] %>%
  pivot_wider(
    names_from = Diet,   # column to spread into new coZZlumns
    values_from = Est      # column containing values
  ) %>% as.data.frame()

rownames(coef_dp) <- coef_dp$Name
coef_dp <- coef_dp[,-1]

pval_dp <- res_dp_r_use[,c(2:3,8:8)] %>%
  pivot_wider(
    names_from = Diet,   # column to spread into new columns
    values_from = sig      # column containing values
  ) %>% as.data.frame()

rownames(pval_dp) <- pval_dp$Name
pval_dp <- pval_dp[,-1]

rm(list = ls())

#--------------------------------------------------------------------------------------------
#
#           chunk2: plotting
#
#--------------------------------------------------------------------------------------------
#read sample data
load("signature_sample.RData")

#plot elastic net coefficient
coef <- elas
coef <- arrange(coef,desc(AMED))
rownames(coef) <-  coef$name
coef <- coef[,-1]

for (i in 1:8){
  coef[,i] <- as.numeric(coef[,i])
}

coef <- round(coef,3)

use1 <- coef

col_fun = colorRamp2(
  breaks = seq(1, -1, length.out = 9),  # Adjust range based on your data
  colors = brewer.pal(9, "RdBu")
)

ht1 <- Heatmap(use1, 
        col=col_fun,
        rect_gp = gpar(col = "grey50", lwd = 1),
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 6,family="calibri"),  
        #column_split = c(rep('A',each=5),rep('B',each=ncol(coef)-5)),
        show_row_dend = FALSE,
        row_order = order(as.numeric(gsub("row", "", rownames(use1)))),
        #row_split = rep(LETTERS[1:length(table(annot_sp_sig$phylum))],table(annot_sp_sig$phylum)[1:length(table(annot_sp_sig$phylum))]),
        row_title = NULL,
        na_col = "white",
        column_names_side = "top", 
        #left_annotation = row_ha,
        column_names_rot = 90,
        column_names_gp = gpar(fontsize = 7,family="calibri"),
        show_column_dend = FALSE,
        column_order = order(as.numeric(gsub("column", "", colnames(use1)))),
        #column_split = LETTERS[1:dim(htmap_coef)[2]],
        column_title = NULL,
        heatmap_legend_param = list(title = "Scaled EN coefficient"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          v <- coef[i, j]
          if (!is.na(v)) {
            grid.text(sprintf("%.3f", v), x, y,gp = gpar(fontsize = 6, family="calibri",color="grey30"))
          }
        }
)

#plot association of individual metabolites in the signatures and each dietary pattern 
use2 <- lm[match(rownames(use1),lm$Metabolite),]
rownames(use2) <- use2$Metabolite
use2 <- use2[,2:9]

col_fun2 = colorRamp2(
  breaks = seq(0.2, -0.2, length.out = 9),  # Adjust range based on your data
  colors = brewer.pal(9, "RdBu")
)

ht2 <- Heatmap(use2, 
        col=col_fun2,
        rect_gp = gpar(col = "grey80", lwd = 1),
        row_names_side = "left",
        show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 8,fontfamily = "Calibri"),  
        #column_split = c(rep('A',each=8),rep('B',each=1)),
        show_row_dend = FALSE,
        row_order = order(as.numeric(gsub("row", "", rownames(use2)))),
        #row_split = rep(LETTERS[1:length(table(annot_sp_sig$phylum))],table(annot_sp_sig$phylum)[1:length(table(annot_sp_sig$phylum))]),
        row_title = NULL,
        na_col = "white",
        column_names_side = "top", 
        column_gap = unit(c(2), "mm"),
        show_heatmap_legend = TRUE,
        column_names_rot = 90,
        column_names_gp = gpar(fontsize = 8,fontfamily = "Calibri"),
        show_column_dend = FALSE,
        column_order = order(as.numeric(gsub("column", "", colnames(use2)))),
        #column_split = LETTERS[1:dim(htmap_coef)[2]],
        column_title = NULL,
        heatmap_legend_param = list(title = "LM coefficient")
)

#plot association with incident T2D
use3 <- lm[match(rownames(use1),lm$Metabolite),]
rownames(use3) <- use3$Metabolite
use3 <- use3[,10:10] %>% as.data.frame()
rownames(use3) <- rownames(use2)

col_fun3 = colorRamp2(
  breaks = seq(1.6,0.4, length.out = 9),  # Adjust range based on your data
  colors = brewer.pal(9, "RdBu")

ht3 <- Heatmap(use3, 
        col=col_fun3,
        rect_gp = gpar(col = "grey80", lwd = 1),
        row_names_side = "left",
        show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 8,fontfamily = "Calibri"),  
        #column_split = c(rep('A',each=8),rep('B',each=1)),
        show_row_dend = FALSE,
        row_order = order(as.numeric(gsub("row", "", rownames(use3)))),
        #row_split = rep(LETTERS[1:length(table(annot_sp_sig$phylum))],table(annot_sp_sig$phylum)[1:length(table(annot_sp_sig$phylum))]),
        row_title = NULL,
        na_col = "white",
        column_names_side = "top", 
        column_gap = unit(c(2), "mm"),
        column_names_rot = 90,
        column_names_gp = gpar(fontsize = 8,fontfamily = "Calibri"),
        show_column_dend = FALSE,
        show_heatmap_legend = TRUE,
        column_order = order(as.numeric(gsub("column", "", colnames(use3)))),
        #column_split = LETTERS[1:dim(htmap_coef)[2]],
        column_title = NULL,
        heatmap_legend_param = list(title = "With T2D")
)

ht_list <- ht1 + ht2 + ht3

#save the figure
png("Figure3.png",width = 2600, height = 3000, res = 300)
draw(ht_list,gap = unit(4, "mm"),
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

rm(list = ls())
