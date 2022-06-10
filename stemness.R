# 干性指数自动计算：提供分组信息即可

rm(list = ls())
pacman::p_load(tidyverse,export,ggpubr)
allrisk <- data.table::fread("lasso500/allRisk.txt")

sample_group <- allrisk %>% select(id,risk) %>% 
  rename(sample = id, group = risk) %>% 
  mutate(sample = paste0(sample,"-01A")) %>% arrange(group)
head(sample_group)
# load("../data/IMMUNE_data/TIDE.Rdata")


# 导入文件 干性指数ZMN和cell--------------------------------------------------------------------
# 记得先备好sample_group,检查一下代码
savefile <- "figure/stemness.pptx"

Plot_Stemness(sample_group = sample_group,savefile = "figure/stemness.pptx")


Plot_Stemness=function(sample_group = null, savefile = "stemness.pptx"){

load("../mitoGENE/outdata/step11_stemness.rds")

# 1.赵队stemness plot ------------------------------------------------------------
  allid <- sample_group$sample
  sameid <- intersect(allid, stem_ZMN$Sample)
  length(sameid)
  
  stem1 <- stem_ZMN[match(sameid, stem_ZMN$Sample),]
  stem1$group <- sample_group$group[match(stem1$Sample,table = sample_group$sample,nomatch = 0)]
  stem1$group <- factor(stem1$group,levels = c("low","high")) # 设置factor
  

library(ggpubr)
p1 <- ggviolin(x="group",y="stemnessScore", fill = "group",data = stem1, size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(group = "group",fill = "white"))+
  xlab("")+
  ylab("")+
  # scale_y_continuous(limits =c(0,1.5),breaks = seq(0,1.5,0.3))+ 
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y.npc = "top",size = 8)+
  theme_classic2()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "stemnessScore")+
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))+
  theme(legend.text = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 14, face = "bold"))

export::graph2ppt(x= p1,file = savefile,width = 6,height = 5, append = T)


# 2.1 mRNAsi plot----------------------------------------------------------------
allid <- sample_group$sample
sameid <- intersect(allid, stem_mRNAsi$Sample)
length(sameid)

stem2 <- stem_mRNAsi[match(sameid, stem_mRNAsi$Sample),]
stem2$group <- sample_group$group[match(stem2$Sample,table = sample_group$sample,nomatch = 0)] 
stem2$group <- factor(stem2$group,levels = c("low","high")) # 设置factor
stem2 %>% dplyr::group_by(group) %>% dplyr::summarise(mean =mean(mRNAsi))
# A tibble: 2 x 2
# group  mean
# <fct> <dbl>
#   1 low   0.302
# 2 high  0.350
library(ggpubr)
p2 <- ggviolin(x="group",y="mRNAsi", fill = "group",data = stem2, size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(group = "group",fill = "white"))+
  xlab("")+
  ylab("")+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y.npc = "top",size = 8)+
  theme_classic2()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "mRNAsi Plot")+
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))+
  theme(legend.text = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 14, face = "bold"))

export::graph2ppt(x=p2,file = savefile,width = 6,height = 5, append = T)

# 2.2 EREG.mRNAsi plot ---------------------------------------------------------------------
library(ggpubr)
p3 <- ggviolin(x="group",y="EREG.mRNAsi", fill = "group",data = stem2, size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(group = "group",fill = "white"))+
  xlab("")+
  ylab("")+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y.npc = "top",size = 8)+
  theme_classic2()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "EREG.mRNAsi Plot")+
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))+
  theme(legend.text = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 14, face = "bold"))
export::graph2ppt(x=p3,file = savefile,width = 6,height = 5, append = T)
# 2.1 mDNAsi plot------------------------------------------------------------------
allid <- sample_group$sample
sameid <- intersect(allid, stem_mDNAsi$Sample)
length(sameid)

stem3 <- stem_mDNAsi[match(sameid, stem_mDNAsi$Sample),]
stem3$group <- sample_group$group[match(stem3$Sample,table = sample_group$sample,nomatch = 0)] 
stem3$group <- factor(stem3$group,levels = c("low","high")) # 设置factor

# 2.1 mDNAsi plot ---------------------------------------------------------------
p4 <- ggviolin(x="group",y="mDNAsi", fill = "group",data = stem3, size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(group = "group",fill = "white"))+
  xlab("")+
  ylab("")+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y.npc = "top",size = 8)+
  theme_classic2()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "mDNAsi Plot")+
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))+
  theme(legend.text = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 14, face = "bold"))

export::graph2ppt(x=p4,file = savefile,width = 6,height = 5, append = T)
}



