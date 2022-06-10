
rm(list = ls())
library(dplyr)
library(ggplot2)


load(file = "../mitoGENE/outdata/LUAD_TCGA_fpkm_symbol_maxgene.Rdata")
fpkm1 <- apply(fpkm, 2, FUN = function(x){2^x - 1}) # 去log化
str(fpkm1)
fpkm1 <- as.data.frame(fpkm1)

# data.table::fwrite(fpkm1, file = "cibersort/DATA.txt",sep = "\t",col.names = T,row.names = T ) # 写出cibersort需要的DATA文件
{
  # data.table::fwrite(fpkm1[,1:150],file = "../maf分组/cibersort/DATA1.txt",sep = "\t",col.names = T,row.names = T ) # 写出cibersort需要的DATA文件
  # data.table::fwrite(fpkm1[,151:300],file = "../maf分组/cibersort/DATA2.txt",sep = "\t",col.names = T,row.names = T ) # 写出cibersort需要的DATA文件
  # data.table::fwrite(fpkm1[,301:450],file = "../maf分组/cibersort/DATA3.txt",sep = "\t",col.names = T,row.names = T ) # 写出cibersort需要的DATA文件
  # data.table::fwrite(fpkm1[,451:585],file = "../maf分组/cibersort/DATA4.txt",sep = "\t",col.names = T,row.names = T ) # 写出cibersort需要的DATA文件
  # tmp1 <- readxl::read_excel("C://Users/JinXing/Desktop/markergenes/EPI-markers.xlsx")
  # tmp2 <- readxl::read_excel("C://Users/JinXing/Desktop/markergenes/immu-marker.xlsx")
  # tmp3 <- readxl::read_excel("C://Users/JinXing/Desktop/markergenes/stomal-marker.xlsx")
  # tmp4 <- data.table::fread("C://Users/JinXing/Desktop/markergenes/Fig1_DE_Lung_atlas_epithelial.csv")
  # tmp5 <- data.table::fread("C://Users/JinXing/Desktop/markergenes/Fig2_DE_Lung_atlas_immune.csv")
  # tmp6 <- data.table::fread("C://Users/JinXing/Desktop/markergenes/immuTmarker.csv")
  # 
  # library(dplyr)
  # genes <- c(tmp1$Gene,tmp2$Gene,tmp3$Gene,tmp4$gene,tmp5$gene,tmp6$gene) %>% unique()
  # genes
  # save(genes,file = "../maf分组/cibersort/marker_genes_all.Rdata")
  # fpkm1 <- as.data.frame(fpkm1) 
  # fpkm1[1:4,1:4]
  # fpkm1 <- fpkm1[intersect(genes,rownames(fpkm1)),]
  # dim(fpkm1)
  # fpkm1[1:4,1:4]
}

# Cibersort
source("cibersort/Cibersort.R")

# Define LM22 file
exp.file <- "../maf分组/cibersort/DATA.txt"
LM22.file <- "../maf分组/cibersort/LM20.txt"
LM10.file <- "../maf分组/cibersort/LM20_czc.txt"
LM50.file <- "cibersort/LM20_czc_origcluster.txt"
TME.results <-  CIBERSORT(LM22.file, exp.file, perm = 1000, QN = TRUE) 
TME.results2 <- CIBERSORT(LM10.file,exp.file,perm=1000,QN = T)
TME.results3 <- CIBERSORT(LM50.file,exp.file,perm=1000,QN = T)

re_hyw <- TME.results[,1:6] %>% t()
re_czc <- TME.results2[,1:11] %>% t()
re_clu_czc <- TME.results3[,1:23] %>% t()

# tme结果作热图 -----------------------------------------------------------------
# TME.results <- read.table("TME.results.output.txt",
#                           sep = "\t",row.names = 1,
#                           header = T,check.names = F) %>% t()

# re <- TME.results[,-(23:25)]

re <- re_clu_czc %>% t()

library(pheatmap)
pheatmap(re,show_rownames = F,cluster_rows = F,show_colnames = F) 
# 热图丑得不堪入目，要去掉在大部分样本中都为0的免疫细胞 
k <- apply(re,2,function(x){sum(x == 0) < nrow(TME.results3)/2})
table(k)

re2 <- as.data.frame(t(re[,k]))

# allid <- read.csv("outdata/allid_493.txt") %>% .[,2]
# phe3 <- read.csv("outdata/phe_finish.csv",row.names = 1)
# rownames(phe3) <-  NULL
# phe3 <- phe3[match(x = allid, table = phe3$sample), ]
# rownames(phe3) <- phe3$sample
# phe3_anno <- phe3 %>% dplyr::select(group, age, gender, smoke_group,ajcc_T, ajcc_N, radiotherapy)
# 
# an <- phe3_anno
# 
# # an = data.frame(group = Group,
# #                 row.names = colnames(exp))
# pheatmap(re2,scale = "row",
#          show_colnames = F,
#          annotation_col = an,
#          cluster_cols = F,
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(50))


# 直方图每个患者的免疫细胞比例 ---------------------------------------------------------------------
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(11,"RdBu"))
# mypalette <- brewer.pal(6,"RdBu")
colnames(re) <- str_replace(colnames(re),pattern = "RNA\\.","Cluster")
save(re,file = "cibersort/cibersort_results_czc_cluster.Rdata")
dat <- re %>% as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>% 
  tidyr::gather(key = Cell_type,value = Proportion,-Sample) %>% 
  mutate(Cell_type = factor(Cell_type,levels = paste0('Cluster',c(0:21,24))))
# dat$group <- an$group[match(dat$Sample, rownames(an))]


library(ggplot2)
ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell_type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(ncol(re)))
library(export)
graph2ppt(file = "figure/cibersort_cells.pptx",append = T,aspect = 1)


# 箱线图展示免疫细胞之间的比较 -------------------------------------------------------------------
ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(ncol(re)))


# 加上顺序
a = dat %>% 
  group_by(Cell_type) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)

dat$Cell_type = factor(dat$Cell_type,levels = a)

ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(ncol(re)))
graph2pdf(file = "figure/f3a.pdf")
######比较两组,肿瘤和非肿瘤##########
# dat$group <- ifelse(endsWith(dat$Sample,suffix = "11A"), "")
library(stringr)
dat$group = ifelse(as.numeric(str_sub(dat$Sample,14,15))<10,"tumor","normal")
library(ggpubr)
ggplot(dat,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,vjust = 0.5))+
  scale_fill_manual(values = c("blue","red"))+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),method = "kruskal.test")

graph2ppt(file = "figure/cibersort_cells.pptx",append = T)
graph2pdf(file = "figure/f3a.pdf")

# 2.用现成的cibersort -----------------------------------------------------------
data <- readr::read_csv("../maf分组/cibersort/cibersortscore.csv")
data <- data %>% 
  filter(CancerType == "LUAD") %>% 
  select(-c(25:27))
data <- data %>% 
  mutate(SampleID = substring(data$SampleID,first = 1,last = 16)) 
data$SampleID <- str_replace_all(string = data$SampleID,pattern = "\\.",replacement = "-")

data1 <- data
names <- data1[,1] %>% as.data.frame() %>% .[,1]
data1 <- data1[, -c(1:2)] %>% as.matrix()
rownames(data1) <- names
data1 <- limma::avereps(data1)
# data1 <- data1[allid, ]
k <- apply(data1,2,function(x) {sum(x == 0) < nrow(data1)/2})
table(k)

data2 <- as.data.frame(t(data1[,k]))

pheatmap(data2,scale = "row",
         show_colnames = F,
         # annotation_col = an,
         cluster_cols = T,cluster_rows = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

graph2ppt(file = "figure/cibersort_cells.pptx",append = T)


dat2 <- data1 %>% as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>% 
  tidyr::gather(key = Cell_type,value = Proportion,-Sample)
dat2$group <- ifelse(as.numeric(str_sub(dat2$Sample,14,15))<10,"tumor","normal")

ggplot(dat2,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell_type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))


ggplot(dat2,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))
# graph2ppt(file = "output/plots/cibersorts22cells",append = T,)

# 加上顺序
# a2 = dat2 %>% 
#   group_by(Cell_type) %>% 
#   summarise(m = median(Proportion)) %>% 
#   arrange(desc(m)) %>% 
#   pull(Cell_type)
# 
# dat2$Cell_type = factor(dat2$Cell_type,levels = a2)

ggplot(dat2,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))
library(ggpubr)
ggplot(dat2,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=70,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(20,1)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),hide.ns = T,method = "kruskal.test")

graph2ppt(file = "figure/cibersort_cells.pptx",append = T)

compare_means(Proportion ~ group,group.by = "Cell_type", data = dat2)


# 只选有意义的作图
dat2sig <- dat2 %>% 
  filter(Cell_type %in% c("Dendritic.cells.resting", 
                          "Mast.cells.resting","Monocytes",
                          "T.cells.CD4.memory.resting",
                          "T.cells.follicular.helper",
                          "T.cells.regulatory..Tregs."))

ggplot(dat2sig,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=0,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(20,1)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),hide.ns = T,method = "kruskal.test")

graph2ppt(file = "figure/cibersort_cells.pptx",append = T)



# 3.czc数据 -----------------------------------------------------------------

re <- re_czc %>% t()

library(pheatmap)
pheatmap(re,show_rownames = F,cluster_rows = F) 
# 热图丑得不堪入目，要去掉在大部分样本中都为0的免疫细胞 
k <- apply(re,2,function(x) {sum(x == 0) < nrow(TME.results)/2})
table(k)

re2 <- as.data.frame(t(re[,k]))

# allid <- read.csv("outdata/allid_493.txt") %>% .[,2]
# phe3 <- read.csv("outdata/phe_finish.csv",row.names = 1)
# rownames(phe3) <-  NULL
# phe3 <- phe3[match(x = allid, table = phe3$sample), ]
# rownames(phe3) <- phe3$sample
# phe3_anno <- phe3 %>% dplyr::select(group, age, gender, smoke_group,ajcc_T, ajcc_N, radiotherapy)
# 
# an <- phe3_anno
# 
# # an = data.frame(group = Group,
# #                 row.names = colnames(exp))
# pheatmap(re2,scale = "row",
#          show_colnames = F,
#          annotation_col = an,
#          cluster_cols = F,
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(50))


# 直方图每个患者的免疫细胞比例 ---------------------------------------------------------------------
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(3,"RdBu"))
# mypalette <- brewer.pal(6,"RdBu")

dat <- re %>% as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>% 
  tidyr::gather(key = Cell_type,value = Proportion,-Sample)
# dat$group <- an$group[match(dat$Sample, rownames(an))]


library(ggplot2)
ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell_type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(11))
library(export)
graph2ppt(file = "figure/cibersort_cells.pptx",append = T)


# 箱线图展示免疫细胞之间的比较 -------------------------------------------------------------------
ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(11))


# 加上顺序
a = dat %>% 
  group_by(Cell_type) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)

dat$Cell_type = factor(dat$Cell_type,levels = a)

ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(11))

######比较两组,肿瘤和非肿瘤##########
# dat$group <- ifelse(endsWith(dat$Sample,suffix = "11A"), "")
library(stringr)
dat$group = ifelse(as.numeric(str_sub(dat$Sample,14,15))<10,"tumor","normal")
library(ggpubr)
ggplot(dat,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,vjust = 0.5))+
  scale_fill_manual(values = mypalette(11)[c(11,1)])+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),method = "kruskal.test")

graph2ppt(file = "figure/cibersort_cells.pptx",append = T)
graph2pdf(file = "figure/f3a.pdf")


