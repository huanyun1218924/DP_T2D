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
#           chunk1: estimate association with dietary patterns
#
#--------------------------------------------------------------------------------------------
#define function for use
inormal <- function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

#read annotation file
anno <- fread("Annotation_final.csv")
anno[218:227,3] <- "Other lipids"

met_t2d <- fread("Met_T2D.csv") %>% as.data.frame()

#read LVS data
load("Processed_LVS_Metabolome.RData")

dp <- c("amed_av","ahei_av","dash_av","pdi_av","hpdi_av","updi_av","edip_av","edih_av")
var <- met_t2d$HMDB

res = data.frame(Diet=NA,HMDB=NA,Est=NA,SE=NA,P=NA)

x = 0

for (i in var){
  for (j in dp){
    data_use <- lvs
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

res <- right_join(met_t2d[,c("HMDB","name")],res,by="HMDB")
res_dp <- res

#unify the results on dp-met
res_dp$name <- gsub("2-Hydroxy-3-methylbutyric acid or 3-Hydroxyisovaleric acid","2-Hydroxy-3-methylbutyric acid",res_dp$name)
res_dp$name <- gsub("3-Carboxy-4-methyl-5-propyl-2-furanpropanoic acid","CMPF",res_dp$name)
res_dp$name <- gsub("Erythronic acid or Threonic acid","Erythronic acid",res_dp$name)
res_dp$name <- gsub("Glycodeoxycholic acid or Glycochenodeoxycholic acid","Glycodeoxycholic acid",res_dp$name)
res_dp$name <- gsub("2-Hydroxy-3-methyl-pentanoic acid or 2-Hydroxyisocaproic acid","2-Hydroxy-3-methyl-pentanoic acid",res_dp$name)
res_dp$name <- gsub("2-Hydroxybutyric acid or 3-Hydroxybutyric acid","3-Hydroxybutyric acid",res_dp$name)

res_dp$Diet <- gsub("amed_av","AMED",res_dp$Diet)
res_dp$Diet <- gsub("ahei_av","AHEI",res_dp$Diet)
res_dp$Diet <- gsub("dash_av","DASH",res_dp$Diet)
res_dp$Diet <- gsub("pdi_av","PDI",res_dp$Diet)
res_dp$Diet <- gsub("hpdi_av","hPDI",res_dp$Diet)
res_dp$Diet <- gsub("updi_av","uPDI",res_dp$Diet)
res_dp$Diet <- gsub("edip_av","EDIP",res_dp$Diet)
res_dp$Diet <- gsub("edih_av","EDIH",res_dp$Diet)

res_dp_r1 <- res_dp[which(res_dp$Diet=="AMED"&res_dp$HMDB %in% dat2[which(dat2$Diet=="AMED"),]$HMDB),]; res_dp_r1$FDR <- p.adjust(res_dp_r1$P,method = "fdr")
res_dp_r2 <- res_dp[which(res_dp$Diet=="AHEI"&res_dp$HMDB %in% dat2[which(dat2$Diet=="AHEI"),]$HMDB),]; res_dp_r2$FDR <- p.adjust(res_dp_r2$P,method = "fdr")
res_dp_r3 <- res_dp[which(res_dp$Diet=="DASH"&res_dp$HMDB %in% dat2[which(dat2$Diet=="DASH"),]$HMDB),]; res_dp_r3$FDR <- p.adjust(res_dp_r3$P,method = "fdr")
res_dp_r4 <- res_dp[which(res_dp$Diet=="PDI"&res_dp$HMDB %in% dat2[which(dat2$Diet=="PDI"),]$HMDB),];  res_dp_r4$FDR <- p.adjust(res_dp_r4$P,method = "fdr")
res_dp_r5 <- res_dp[which(res_dp$Diet=="hPDI"&res_dp$HMDB %in% dat2[which(dat2$Diet=="hPDI"),]$HMDB),]; res_dp_r5$FDR <- p.adjust(res_dp_r5$P,method = "fdr")
res_dp_r6 <- res_dp[which(res_dp$Diet=="uPDI"&res_dp$HMDB %in% dat2[which(dat2$Diet=="uPDI"),]$HMDB),]; res_dp_r6$FDR <- p.adjust(res_dp_r6$P,method = "fdr")
res_dp_r7 <- res_dp[which(res_dp$Diet=="EDIP"&res_dp$HMDB %in% dat2[which(dat2$Diet=="EDIP"),]$HMDB),]; res_dp_r7$FDR <- p.adjust(res_dp_r7$P,method = "fdr")
res_dp_r8 <- res_dp[which(res_dp$Diet=="EDIH"&res_dp$HMDB %in% dat2[which(dat2$Diet=="EDIH"),]$HMDB),]; res_dp_r8$FDR <- p.adjust(res_dp_r8$P,method = "fdr")

res_dp_r_use <- rbind(res_dp_r1,res_dp_r2,res_dp_r3,res_dp_r4,res_dp_r5,res_dp_r6,res_dp_r7,res_dp_r8)
res_dp_r_use$sig <- ifelse(res_dp_r_use$FDR <= 0.05 & res_dp_r_use$P < 0.05, "**",ifelse(res_dp_r_use$FDR > 0.05 & res_dp_r_use$P < 0.05, "*", ""))

coef_dp <- res_dp_r_use[,c(2:4)] %>%
  pivot_wider(
    names_from = Diet,   # column to spread into new coZZlumns
    values_from = Est      # column containing values
  ) %>% as.data.frame()

rownames(coef_dp) <- coef_dp$name
coef_dp <- coef_dp[,-1]

pval_dp <- res_dp_r_use[,c(2:3,8:8)] %>%
  pivot_wider(
    names_from = Diet,   # column to spread into new columns
    values_from = sig      # column containing values
  ) %>% as.data.frame()

rownames(pval_dp) <- pval_dp$name
pval_dp <- pval_dp[,-1]

#distribution of coefficients
elas <- fread("Coeff_EN.csv") %>% as.data.frame()
lm <- fread("Coeff_LM.csv") %>% as.data.frame()

elas_long <- elas %>%
  pivot_longer(
    cols = names(elas)[2:9],
    names_to = "Trait",
    values_to = "coefficient"
  )

#re-name metabolites
met_t2d$name <- gsub("2-Hydroxy-3-methylbutyric acid or 3-Hydroxyisovaleric acid","2-Hydroxy-3-methylbutyric acid",met_t2d$name)
met_t2d$name <- gsub("3-Carboxy-4-methyl-5-propyl-2-furanpropanoic acid","CMPF",met_t2d$name)
met_t2d$name <- gsub("Erythronic acid or Threonic acid","Erythronic acid",met_t2d$name)
met_t2d$name <- gsub("Glycodeoxycholic acid or Glycochenodeoxycholic acid","Glycodeoxycholic acid",met_t2d$name)
met_t2d$name <- gsub("2-Hydroxy-3-methyl-pentanoic acid or 2-Hydroxyisocaproic acid","2-Hydroxy-3-methyl-pentanoic acid",met_t2d$name)
met_t2d$name <- gsub("2-Hydroxybutyric acid or 3-Hydroxybutyric acid","3-Hydroxybutyric acid",met_t2d$name)

anno_use <- left_join(met_t2d[,c("HMDB","name")],anno[,c("HMDB","Subclass")],by="HMDB")

elas$name <- gsub("2-Hydroxy-3-methylbutyric acid or 3-Hydroxyisovaleric acid","2-Hydroxy-3-methylbutyric acid",elas$name)
elas$name <- gsub("3-Carboxy-4-methyl-5-propyl-2-furanpropanoic acid","CMPF",elas$name)
elas$name <- gsub("Erythronic acid or Threonic acid","Erythronic acid",elas$name)
elas$name <- gsub("Glycodeoxycholic acid or Glycochenodeoxycholic acid","Glycodeoxycholic acid",elas$name)
elas$name <- gsub("2-Hydroxy-3-methyl-pentanoic acid or 2-Hydroxyisocaproic acid","2-Hydroxy-3-methyl-pentanoic acid",elas$name)
elas$name <- gsub("2-Hydroxybutyric acid or 3-Hydroxybutyric acid","3-Hydroxybutyric acid",elas$name)

lm$Metabolite <- gsub("2-Hydroxy-3-methylbutyric acid or 3-Hydroxyisovaleric acid","2-Hydroxy-3-methylbutyric acid",lm$Metabolite)
lm$Metabolite <- gsub("3-Carboxy-4-methyl-5-propyl-2-furanpropanoic acid","CMPF",lm$Metabolite)
lm$Metabolite <- gsub("Erythronic acid or Threonic acid","Erythronic acid",lm$Metabolite)
lm$Metabolite <- gsub("Glycodeoxycholic acid or Glycochenodeoxycholic acid","Glycodeoxycholic acid",lm$Metabolite)
lm$Metabolite <- gsub("2-Hydroxy-3-methyl-pentanoic acid or 2-Hydroxyisocaproic acid","2-Hydroxy-3-methyl-pentanoic acid",lm$Metabolite)
lm$Metabolite <- gsub("2-Hydroxybutyric acid or 3-Hydroxybutyric acid","3-Hydroxybutyric acid",lm$Metabolite)

#--------------------------------------------------------------------------------------------
#
#           chunk2: plotting
#
#--------------------------------------------------------------------------------------------
coef <- elas
coef <- arrange(coef,desc(AMED))
rownames(coef) <-  coef$name
coef <- coef[,-1]
col_fun = colorRamp2(
  breaks = seq(1, -1, length.out = 9),  # Adjust range based on your data
  colors = brewer.pal(9, "RdBu")
)

anno_use <- anno_use[match(rownames(use),anno_use$name),]

row_ha<-rowAnnotation(foo = anno_use$Subclass,
                      col = list(foo = c("Amino acids" = "#4F81BD", "Carbohydrates" = "#8cd2c8", "Glycerophospholipids" = "#fa826e", "Glycerolipids" = "#fab464",
                                         "Sphingolipids" = "#b4dc64", "Fatty acids" = "#facde6", "Other lipids" = "#95a2ff", "Acylcarnitines" = "#beb9dc",
                                         "Nucleotides" = "#B15928", "Cofactors and vitamins" = "#be82be", "Xenobiotics" = "#5F7530"),gp = gpar(col = "black")),
                      annotation_name_side = "top",annotation_name_rot=90,
                      annotation_legend_param = list(foo = list(title = "Subclass",labels = c("Amino acids","Carbohydrates","Glycerophospholipids","Glycerolipids","Sphingolipids",
                                                                                              "Fatty acids","Other lipids","Acylcarnitines","Nucleotides","Cofactors and vitamins","Xenobiotics"))),
                      annotation_label = "Subclass")


for (i in 1:8){
  coef[,i] <- as.numeric(coef[,i])
}

coef <- round(coef,3)

use1 <- coef

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

use2 <- coef_dp[match(rownames(use1),rownames(coef_dp)),]
use2_2 <- pval_dp[match(rownames(use1),rownames(pval_dp)),]

col_fun2 = colorRamp2(
  breaks = seq(0.2, -0.2, length.out = 9),  # Adjust range based on your data
  colors = brewer.pal(9, "RdBu")
)

use2 <- coef_dp[match(rownames(use_ordered),rownames(coef_dp)),]
use2_2 <- pval_dp[match(rownames(use_ordered),rownames(pval_dp)),]

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
        heatmap_legend_param = list(title = "LM coefficient"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          v <- use2_2[i, j]
          if (!is.na(v)) {
            grid.text(sprintf(v), x, y,gp = gpar(fontsize = 8, fontfamily="calibri",color="grey30"))
          }
        }
)

lm_filtered <- met_t2d[match(rownames(use_ordered),met_t2d$name),]

use3 <- lm_filtered[,c("name","beta_t2d")]
rownames(use3) <- use3$name
use3 <- use3[,-1] %>% as.data.frame()
names(use3) <- "T2D"
use3$T2D <- round(exp(use3$T2D),2)
rownames(use3) <- lm_filtered$name

use3_2 <- lm_filtered[,c("name","t2d_sig")]
rownames(use3_2) <- use3_2$name
use3_2 <- use3_2[,-1] %>% as.data.frame()
names(use3_2) <- "T2D"
rownames(use3_2) <- lm_filtered$name

col_fun3 = colorRamp2(
  breaks = seq(1.6,0.4, length.out = 9),  # Adjust range based on your data
  colors = brewer.pal(9, "RdBu")
)

use3 <- use3[match(rownames(use_ordered),rownames(use3)),]
use3_2 <- use3_2[match(rownames(use_ordered),rownames(use3_2)),]

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
        heatmap_legend_param = list(title = "With T2D"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          v <- use3_2[i, j]
          if (!is.na(v)) {
            grid.text(sprintf(v), x, y,gp = gpar(fontsize = 8, fontfamily="calibri",color="grey30"))
          }
        }
)

ht_list <- ht1 + ht2 + ht3

png("Figure3.png",width = 2600, height = 3000, res = 300)
draw(ht_list,gap = unit(4, "mm"),
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()