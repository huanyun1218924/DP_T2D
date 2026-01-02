# --------------------------------------------------------------------------------
#Title: Dietary quality scores, metabolic signature, and cardiometabolic diseases
#Purpose: Compare metabolites among guideline-recommended group, plant-based group, and mechanism-driven diet
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
library(VennDiagram)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(patchwork)
library(venn)

#set working dir
setwd("/udd/nhhyu/DP_T2D/FinalOutput")

#read results
load("Figure3.RData")

#--------------------------------------------------------------------------------------------
#
#           step1: Venn plot
#
#--------------------------------------------------------------------------------------------
#read cross-platform signature
sig_list <- read_excel("Metabolic signature of DQSs_Feb 13.xlsx", sheet = "Cross-platform_SOL (n=122)") %>% as.data.frame()

amed1 <- na.omit(sig_list[,1:2]); amed <- data.frame(amed1[,-1]); rownames(amed) <- amed1[,1]
ahei1 <- na.omit(sig_list[,4:5]); ahei <- data.frame(ahei1[,-1]); rownames(ahei) <- ahei1[,1]
dash1 <- na.omit(sig_list[,7:8]); dash <- data.frame(dash1[,-1]); rownames(dash) <- dash1[,1]
opdi1 <- na.omit(sig_list[,10:11]); opdi <- data.frame(opdi1[,-1]); rownames(opdi) <- opdi1[,1]
hpdi1 <- na.omit(sig_list[,13:14]); hpdi <- data.frame(hpdi1[,-1]); rownames(hpdi) <- hpdi1[,1]
updi1 <- na.omit(sig_list[,16:17]); updi <- data.frame(updi1[,-1]); rownames(updi) <- updi1[,1]
edip1 <- na.omit(sig_list[,19:20]); edip <- data.frame(edip1[,-1]); rownames(edip) <- edip1[,1]
edih1 <- na.omit(sig_list[,22:23]); edih <- data.frame(edih1[,-1]); rownames(edih) <- edih1[,1]

names(amed1) <- c("HMDB","coeff_amed")
names(ahei1) <- c("HMDB","coeff_ahei")
names(dash1) <- c("HMDB","coeff_dash")
names(opdi1) <- c("HMDB","coeff_pdi")
names(hpdi1) <- c("HMDB","coeff_hpdi")
names(updi1) <- c("HMDB","coeff_updi")
names(edip1) <- c("HMDB","coeff_edip")
names(edih1) <- c("HMDB","coeff_edih")

a1 <- amed1[-1,]$HMDB
b1 <- ahei1[-1,]$HMDB
c1 <- dash1[-1,]$HMDB
d1 <- opdi1[-1,]$HMDB
e1 <- hpdi1[-1,]$HMDB
f1 <- updi1[-1,]$HMDB
g1 <- edip1[-1,]$HMDB
h1 <- edih1[-1,]$HMDB

dataForVennDiagram1 <- list(a1,b1,c1)
dataForVennDiagram2 <- list(g1,h1)
dataForVennDiagram3 <- list(d1,e1,f1)

#remove intercept and scale coeff
amed2 <- amed1[-1,]

#get all weight
weight <- full_join(amed1,ahei1,by="HMDB") %>% full_join(dash1,by="HMDB") %>% full_join(opdi1,by="HMDB") %>% full_join(hpdi1,by="HMDB") %>% full_join(updi1,by="HMDB") %>% full_join(edip1,by="HMDB") %>% full_join(edih1,by="HMDB")

#venn plot
pdf("Figure3A-2.pdf",height=3.5, width = 3.5)
venn(dataForVennDiagram1, snames = "", opacity = 0.55, ellipse = F, borders = T, box = F, zcolor = c("#00b0f0","#00b050","#ffc000"))
venn(dataForVennDiagram3, snames = "", opacity = 0.55, ellipse = T, borders = T, box = F, zcolor = c("#00b0f0","#00b050","#ffc000"))
venn(dataForVennDiagram2, snames = "", opacity = 0.55, ellipse = T, borders = T, box = F, zcolor = c("#00b0f0","#00b050","#ffc000"))
dev.off()

#--------------------------------------------------------------------------------------------
#
#           step2: Heatmap to show elastic net model coeffiecient
#
#--------------------------------------------------------------------------------------------
#read results
load("Figure3.RData")

#remove long name
data_h[23,2] <- "2-Hydroxy-3-methylbutyric acid"
data_h[22,2] <- "2-Hydroxy-3-methyl-pentanoic acid"
data_h[27,2] <- "Erythronic acid"
data_h[28,2] <- "Glycodeoxycholic acid"
data_h$name <- gsub("3-Carboxy-4-methyl-5-propyl-2-furanpropanoic acid","CMPF",data_h$name)

data_u[17,2] <- "2-Hydroxy-3-methylbutyric acid"
data_u[21,2] <- "Glycodeoxycholic acid"
data_u$name <- gsub("3-Carboxy-4-methyl-5-propyl-2-furanpropanoic acid","CMPF",data_u$name)

data_p[1,2] <- "3-Hydroxybutyric acid"
data_p[17,2] <- "2-Hydroxy-3-methylbutyric acid"
data_p[16,2] <- "2-Hydroxy-3-methyl-pentanoic acid"
data_p[20,2] <- "Erythronic acid"
data_p[21,2] <- "Glycodeoxycholic acid"
data_p$name <- gsub("3-Carboxy-4-methyl-5-propyl-2-furanpropanoic acid","CMPF",data_p$name)

#re-annotate results to show coefficients
coeff1 <- data_h[,c("name","coeff_amed_std","order_h")]; coeff1$DP <- "AMED"
coeff2 <- data_h[,c("name","coeff_ahei_std","order_h")]; coeff2$DP <- "AHEI-2010"
coeff3 <- data_h[,c("name","coeff_dash_std","order_h")]; coeff3$DP <- "DASH"
names(coeff1) <- names(coeff2) <- names(coeff3) <- c("name","Coefficient","Order","DP")
coeff <- rbind(coeff1,coeff2,coeff3) %>% arrange(desc(Order),Coefficient)

#define direction for signature coefficients
coeff$Direction[which(coeff$Coefficient>0)] = -1
coeff$Direction[which(coeff$Coefficient<0)] = 1
coeff$Direction <- as.factor(coeff$Direction)

#arrange the metabolites by order and by coefficient
coeff <- coeff[with(coeff, order(Order, DP, -Direction, -Coefficient)), ]
coeff$Order <- factor(coeff$Order,levels=c("6","5","4","3","2","1","0"))
coeff$Direction <- factor(coeff$Direction,levels=c("1","-1"))
coeff$DP <- factor(coeff$DP,levels=c("AMED","AHEI-2010","DASH"))
coeff$name <- factor(coeff$name,levels=rev(arrange(coeff1, Order, desc(Coefficient))$name))

#annotation plot
g1 <- ggplot(coeff, aes(fill=Order, x=1, y=name)) + 
  geom_bar(position = "stack",stat='identity',width = 1) +
  scale_fill_manual(values=c("#5F7530","#4F81BD","#B15928","#fa826e","#95a2ff","#b4dc64","#be82be"),guide = guide_legend(reverse = TRUE))+
  theme_bw() + 
  theme(strip.text.x=element_blank(),
        panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_text(color = "black",family = "Calibri"),
        axis.ticks.x=element_blank(),
        legend.position='none') 

#bi-directional bar plot
g2 <- ggplot(coeff, aes(x=name, y=Coefficient,fill=Direction)) +
  geom_bar(position="dodge",stat='identity')+
  scale_fill_manual(values=c("#006EBE","#FA5555"))+
  geom_hline(aes(yintercept = 0),colour="black", linetype="dashed", size=0.5) +
  facet_grid(~DP, scales = "free_x") + 
  theme_bw() + 
  theme(strip.text.x=element_blank(),
        panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black",family = "Calibri"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position='none',
        panel.spacing.x=unit(0.1,"cm"),
        panel.spacing.y=unit(0.05,"cm")) +
  coord_flip()

#heatmap to show association with dietary score
xx <- data_h[,c(12:14)]; xx2 <- data_h[,c("amed_sig2","ahei_sig2","dash_sig2")]
yy <- data_h[,c(62:62)]; yy2 <- data_h[,c("t2d_sig")] %>% as.data.frame(); yy2[75,1] <- ""

g3 <- pheatmap(xx, color=colorRampPalette(c("#006EBE", "white", "#FA5555"))(90), cluster_rows = F,cluster_cols = F,display_numbers = xx2, fontsize = 12, border_color = "black",show_rownames = F,show_colnames = F, legend = F)
g4 <- pheatmap(yy, color=colorRampPalette(c("#006EBE", "white", "#FA5555"))(90), cluster_rows = F,cluster_cols = F,display_numbers = yy2, fontsize = 12, border_color = "black",show_rownames = F,show_colnames = F, legend = F)

#plot EDIP and EDIH
#re-annotate results to show coefficients
coeff1 <- data_u[,c("name","coeff_edip_std","order_u")]; coeff1$DP <- "EDIP"
coeff2 <- data_u[,c("name","coeff_edih_std","order_u")]; coeff2$DP <- "EDIH"
names(coeff1) <- names(coeff2) <- c("name","Coefficient","Order","DP")
coeff <- rbind(coeff1,coeff2) %>% arrange(desc(Order))

#define direction for signature coefficients
coeff$Direction[which(coeff$Coefficient>0)] = -1
coeff$Direction[which(coeff$Coefficient<0)] = 1
coeff$Direction <- as.factor(coeff$Direction)

#define levels
coeff <- coeff[with(coeff, order(Order, DP, -Direction, -Coefficient)), ]
coeff$Direction <- factor(coeff$Direction,levels=c("1","-1"))
coeff$DP <- factor(coeff$DP, levels=c("EDIP","EDIH"))
coeff$Order <- factor(coeff$Order, levels=c("2","1","0"))
coeff$name <- factor(coeff$name,levels=rev(arrange(coeff1, Order, desc(Coefficient))$name))

#annotation plot
g5 <- ggplot(coeff, aes(fill=Order, x=1, y=name)) + 
  geom_bar(position = "stack",stat='identity',width = 1) +
  scale_fill_manual(values=c("#5F7530","#4F81BD","#B15928","#fa826e","#95a2ff","#b4dc64","#be82be"),guide = guide_legend(reverse = TRUE))+
  theme_bw() + 
  theme(strip.text.x=element_blank(),
        panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_text(color = "black",family = "Calibri"),
        axis.ticks.x=element_blank(),
        legend.position='none') 

#bi-directional bar plot
g6 <- ggplot(coeff, aes(x=name, y=Coefficient,fill=Direction)) +
  geom_bar(position="dodge",stat='identity')+
  scale_fill_manual(values=c("#006EBE","#FA5555"))+
  geom_hline(aes(yintercept = 0),colour="black", linetype="dashed", size=0.5) +
  facet_grid(~DP, scales = "free_x") + 
  theme_bw() + 
  theme(strip.text.x=element_blank(),
        panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black",family = "Calibri"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position='none',
        panel.spacing.x=unit(0.1,"cm"),
        panel.spacing.y=unit(0.05,"cm")) +
  coord_flip()

#heatmap to show association with dietary score
xx <- data_u[,c(18:19)]; xx2 <- data_u[,c("edip_sig2","edih_psig")]
yy <- data_u[,c(62:62)]; yy2 <- data_u[,c("t2d_sig")] %>% as.data.frame(); yy2[75,1] <- ""

g7 <- pheatmap(xx, color=colorRampPalette(c("#006EBE", "white", "#FA5555"))(90), cluster_rows = F,cluster_cols = F,display_numbers = xx2, fontsize = 12, border_color = "black",show_rownames = F,show_colnames = F, legend = F)
g8 <- pheatmap(yy, color=colorRampPalette(c("#006EBE", "white", "#FA5555"))(90), cluster_rows = F,cluster_cols = F,display_numbers = yy2, fontsize = 12, border_color = "black",show_rownames = F,show_colnames = F, legend = F)

#plot PDIs
#re-annotate results to show coefficients
coeff1 <- data_p[,c("name","coeff_pdi_std","order_p")]; coeff1$DP <- "PDI"
coeff2 <- data_p[,c("name","coeff_hpdi_std","order_p")]; coeff2$DP <- "hPDI"
coeff3 <- data_p[,c("name","coeff_updi_std","order_p")]; coeff3$DP <- "uPDI"
names(coeff1) <- names(coeff2) <- names(coeff3) <- c("name","Coefficient","Order","DP")
coeff <- rbind(coeff1,coeff2,coeff3) %>% arrange(desc(Order))

#define direction for signature coefficients
coeff[which(coeff$Coefficient>0),"Direction"] <- "-1"
coeff$Direction[which(coeff$Coefficient<=0)] = 1
coeff$Direction <- as.factor(coeff$Direction)

#define levels
coeff <- coeff[with(coeff, order(Order, DP, -Direction, -Coefficient)), ]
coeff$Order <- factor(coeff$Order,levels=c("6","5","4","3","2","1","0"))
coeff$Direction <- factor(coeff$Direction,levels=c("1","-1"))
coeff$DP <- factor(coeff$DP,levels=c("PDI","hPDI","uPDI"))
coeff$name <- factor(coeff$name,levels=rev(arrange(coeff1, Order, desc(Coefficient))$name))

#annotation plot
g9 <- ggplot(coeff, aes(fill=Order, x=1, y=name)) + 
  geom_bar(position = "stack",stat='identity',width = 1) +
  scale_fill_manual(values=c("#5F7530","#4F81BD","#B15928","#fa826e","#95a2ff","#b4dc64","#be82be"),guide = guide_legend(reverse = TRUE))+
  theme_bw() + 
  theme(strip.text.x=element_blank(),
        panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_text(color = "black",family = "Calibri"),
        axis.ticks.x=element_blank(),
        legend.position='none') 

#bi-directional bar plot
g10 <- ggplot(coeff, aes(x=name, y=Coefficient,fill=Direction)) +
  geom_bar(position="dodge",stat='identity')+
  scale_fill_manual(values=c("#006EBE","#FA5555"))+
  geom_hline(aes(yintercept = 0),colour="black", linetype="dashed", size=0.5) +
  facet_grid(~DP, scales = "free_x") + 
  theme_bw() + 
  theme(strip.text.x=element_blank(),
        panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black",family = "Calibri"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position='none',
        panel.spacing.x=unit(0.1,"cm"),
        panel.spacing.y=unit(0.05,"cm")) +
  coord_flip()

#heatmap to show association with dietary score
xx <- data_p[,c(15:17)]; xx2 <- data_p[,c("pdi_sig2","hpdi_sig2","updi_sig2")]
yy <- data_p[,c(62:62)]; yy2 <- data_p[,c("t2d_sig")] %>% as.data.frame(); yy2[72,1] <- ""

g11 <- pheatmap(xx, color=colorRampPalette(c("#006EBE", "white", "#FA5555"))(90), cluster_rows = F,cluster_cols = F,display_numbers = xx2, fontsize = 13, border_color = "black",show_rownames = F,show_colnames = F, legend = F)
g12 <- pheatmap(yy, color=colorRampPalette(c("#006EBE", "white", "#FA5555"))(90), cluster_rows = F,cluster_cols = F,display_numbers = yy2, fontsize = 13, border_color = "black",show_rownames = F,show_colnames = F, legend = F)

p <- plot_grid(g1,g4$gtable,g2,g9,g12$gtable,g10,g5,g8$gtable,g6,nrow = 1,rel_widths=c(3.5,0.6,4.0,3.5,0.6,4.0,3.18,0.6,2.5))

p1 <- plot_grid(g1,g4$gtable,g2,nrow = 1,rel_widths=c(4.0,0.6,4.0))
p2 <- plot_grid(g9,g12$gtable,g10,nrow = 1,rel_widths=c(4.0,0.6,4.0))
p3 <- plot_grid(g5,g8$gtable,g6,nrow = 1,rel_widths=c(4.0,0.6,4.0))

ggsave("Figure3B1.png", plot = p1, width = 5.3, height = 9, dpi = 1000)
ggsave("Figure3B2.png", plot = p2, width = 5.3, height = 9, dpi = 1000)
ggsave("Figure3B3.png", plot = p3, width = 4.8, height = 9, dpi = 1000)