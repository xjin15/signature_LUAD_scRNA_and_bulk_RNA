rm(list = ls())
pacman::p_load(tidyverse,export,ggpubr)
# tide_immune_data <- tide_immune_data %>% 
#   mutate(cell_type = paste0(cell_type,"A")) %>% 
#   rename(sample = cell_type)
# colnames(tide_immune_data) <- colnames(tide_immune_data) %>% str_replace_all(pattern = " ",replacement = ".")
# tide_immune_data
# save(tide_immune_data,file = "../data/IMMUNE_data/TIDE.Rdata")


# 1.1stormal_score+immune_scorebyextimate --------------------------------------------
library(limma)
library(estimate)
load('../mitoGENE/outdata/LUAD_TCGA_fpkm_symbol_maxgene.Rdata')

# data=fpkm
# #删除正常样品
# data <- data[,endsWith(colnames(data),suffix = "01A")]
# out=data[rowMeans(data)>0,]
# 
# #输出整理后的矩阵文件
# out=rbind(ID=colnames(out), out)
# write.table(out, file="estimate/uniq.symbol.txt", sep="\t", quote=F, col.names=F)
# 
# #运行estimate包
# estimate::filterCommonGenes(input.f="estimate/uniq.symbol.txt", 
#                   output.f="estimate/commonGenes.gct", 
#                   id="GeneSymbol")
# 
# estimate::estimateScore(input.ds="estimate/commonGenes.gct",
#               output.ds="estimate/estimateScore.gct")
# 
# #输出每个样品的打分
# scores=read.table("estimate/estimateScore.gct", skip=2, header=T, check.names=F)
# rownames(scores)=scores[,1]
# scores=t(scores[,3:ncol(scores)])
# rownames(scores)=gsub("\\.", "\\-", rownames(scores))
# # scores=scores[,1:3]
# out=rbind(ID=colnames(scores), scores)
# write.table(out,file="scores.txt",sep="\t",quote=F,col.names=F)

# 1.2读取estimate文件并作图  -----------------------------------------------------------
allrisk <- data.table::fread("lasso500/allRisk.txt")
sample_group <- allrisk %>% select(id,risk) %>% 
  rename(sample = id, group = risk) %>% 
  mutate(sample = paste0(sample,"-01A"))
head(sample_group)

score <- data.table::fread("estimate/scores.txt") %>% rename(sample = ID)
library(limma)
library(ggpubr)


# score1=score[sameSample,,drop=F]
# cluster1=cluster[sameSample,,drop=F]
# colnames(cluster1)=c("Cluster")
# data=cbind(score1, cluster1)
# data$Cluster=paste0("Cluster", data$Cluster)

# #设置颜色
# bioCol=c("#0072B5","#BC3C29","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
# bioCol=bioCol[1:length(levels(factor(data$Cluster)))]



#样品取交集
sameid=intersect(score$sample, sample_group$sample)
data <- score[match(sameid,score$sample),]
data$group <- sample_group$group[match(data$sample,table = sample_group$sample,nomatch = 0)]
data$group <-factor(data$group,levels = c("low","high")) 

type=levels(factor(data$group))
data$group=factor(data$group, levels=type)
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
data$TumorPurity
data %>% dplyr::group_by(group) %>% dplyr::summarise(mean =mean(ImmuneScore))
# group  mean
# <fct> <dbl>
#   1 low   1812.
# 2 high  1204.
data %>% dplyr::group_by(group) %>% dplyr::summarise(mean =mean(TumorPurity))
# A tibble: 2 x 2
# group  mean
# <fct> <dbl>
#   1 low   0.599
# 2 high  0.692


#肿瘤微环境差异分析
for(i in colnames(data)[2:(ncol(data)-1)]){
  #绘制箱线图
  print(i)
  boxplot=ggboxplot(data, x="group", y=i, fill="group",
                    xlab="",
                    ylab=i,
                    legend.title="Cluster",
                    palette=c('#00ccff','#ff3333')
  )+ 
    stat_compare_means(comparisons=my_comparisons)
  print(boxplot)
export::graph2ppt(file= "figure/score.pptx",width = 5, height =7 ,append=T)
}


# 2.2TIDE浸润分数 -------------------------------------------------------------
load("../data/IMMUNE_data/TIDE.Rdata")
allrisk <- data.table::fread("lasso500/allRisk.txt")
sample_group <- allrisk %>% select(id,risk) %>% 
  rename(sample = id, group = risk) %>% 
  mutate(sample = paste0(sample,"-01A"))
head(sample_group)

sameid <- intersect(sample_group$sample,tide_immune_data$sample)
allrisk$id <- paste0(allrisk$id,"-01A")

tide <- tide_immune_data[match(sameid,table = tide_immune_data$sample),]
tide$group <- sample_group$group[match(tide$sample,table = sample_group$sample)]
tide$riskscore <- allrisk$riskScore[match(tide$sample,table = allrisk$id)]
tide <- tide %>% select(sample,group,riskscore,everything())

# colnames(tide) %>% grep(pattern = "NK",value = T) -> nkids
# for reproducibility
set.seed(123)
# library(ggcorrplot)
# cormat <- round(cor(tide[,3:ncol(tide)]),digits =3) # 计算相关系数，默认是pearson方法
# p.mat <- cor_pmat(tide[,3:ncol(tide)])
# head(p.mat[, 1:4])
##               mpg          cyl         disp           hp
## mpg  0.000000e+00 6.112687e-10 9.380327e-10 1.787835e-07
## cyl  6.112687e-10 0.000000e+00 1.803002e-12 3.477861e-09
## disp 9.380327e-10 1.803002e-12 0.000000e+00 7.142679e-08
## hp   1.787835e-07 3.477861e-09 7.142679e-08 0.000000e+00
## drat 1.776240e-05 8.244636e-06 5.282022e-06 9.988772e-03
## wt   1.293959e-10 1.217567e-07 1.222311e-11 4.145827e-05
# cor
# ggcorrplot(cormat, hc.order = TRUE,type = "lower", p.mat = p.mat) # 除非没有显著系数
corFilter= 0.2    #相关系数过滤标准
pvalueFilter=0.05       #p值过滤标准
tide <- tide %>% as.data.frame()

cells <- data.frame()
for(i in 3:ncol(tide)){
    x=as.numeric(tide$riskscore)
    y=as.numeric(tide[,i])
    corT=cor.test(x,y)
    cor=corT$estimate
    pvalue=corT$p.value
    if((cor>corFilter) & (pvalue<pvalueFilter)){
      cells <- rbind(cells,cbind(cell = colnames(tide)[i],cor,pvalue,regulation = "posi"))
    }
    if((cor< -corFilter) & (pvalue<pvalueFilter)){
      cells <- rbind(cells,cbind(cell = colnames(tide)[i],cor,pvalue,regulation = "nega"))
    }
}

cellids <- cells$cell
cellids
cells$short <- cellids %>% str_split(pattern = "_",n = 2,simplify = T) %>% .[,1]

data <- tide %>% select(1:3,all_of(cellids))
data[1:5,1:6]

library(ggpubr)

for (i in 4:ncol(data)) {
  
p <- ggplot(data=data, aes(x=riskscore, y=data[,i]))+
     geom_point(size=4,alpha=0.8,color="#6baed6")+
    stat_smooth(method="lm",formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  theme_bw()+ylab(label = NULL)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_cor(data=data, method = "pearson",size=7,color= "red") + 
  theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15, face = "bold",colour = "gray0")) +
  labs(x = "RiskScore",title = str_split(colnames(data)[i],pattern = "_",simplify = T)[1]) +
  theme(plot.title = element_text(family = "Arial",size = 20, hjust = 0.5))

export::graph2ppt(x = p,file = "figure/cor_analysis.pptx",append = T,width = 6,height =6)

}














