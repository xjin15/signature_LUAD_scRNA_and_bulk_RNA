library(dplyr)
library(ggpubr)
library(ggplot2)

# 1.Tcell dysfunction ------------------------------------------------------
allscore <- data.table::fread("testdata/Results/Tumor_Dysf_Excl_scores/TCGA.LUAD.RNASeq.norm_subtract.OS_base") %>% 
  rename(sample = V1) %>% 
  mutate(sample = paste0(sample,"-01A"))

allrisk <- data.table::fread("lasso500/allRisk.txt") 

sample_group <- allrisk %>% select(id,risk) %>% 
  rename(sample = id, group = risk) %>% 
  mutate(sample = paste0(sample,"-01A"),
         group = factor(group,levels = c("low","high"))) %>% 
  arrange(group)
head(sample_group)

sameid <- intersect(sample_group$sample,allscore$sample)
length(sameid)
data <- allscore[allscore$sample %in% sameid,]
data$group <- sample_group$group[match(x=data$sample,table = sample_group$sample)]
data$riskscore <- allrisk$riskScore[match(x=data$sample,table = paste0(allrisk$id,"-01A"))]

data <- data %>% arrange(group)

ggstatsplot::ggbetweenstats(data,x=group,y=Dysfunction,centrality.plotting = FALSE,results.subtitle = FALSE,
               bf.message = FALSE,pairwise.comparisons = FALSE,)+stat_compare_means()+theme_classic2()+theme(legend.position = "none") 
export::graph2ppt(file="figure/TIDE.pptx",width = 5,height = 5,append=T)

ggstatsplot::ggbetweenstats(data,x=group,y=Exclusion,centrality.plotting = FALSE,results.subtitle = FALSE,
               bf.message = FALSE,pairwise.comparisons = FALSE,)+stat_compare_means()+theme_classic2()+theme(legend.position = "none") 
export::graph2ppt(file="figure/TIDE.pptx",width = 5,height = 5,append=T)

# 2.tide网站预测tide评分 --------------------------------------------------------


tmp2 <- data.table::fread("testdata/tide2.csv")
sameid <-intersect(tmp2$Patient,sample_group$sample)
data2 <- tmp2[tmp2$Patient %in% sameid,]

data2 <- merge(data2,data,by.x="Patient",by.y="sample")

ggstatsplot::ggbetweenstats(data2,x=group,y=TIDE,centrality.plotting = FALSE,results.subtitle = FALSE,
                            bf.message = FALSE,pairwise.comparisons = FALSE,)+stat_compare_means()+theme_classic2()+theme(legend.position = "none") 
export::graph2ppt(file="figure/TIDE.pptx",width = 5,height = 5,append=T)

ggstatsplot::ggbarstats(data2,x=Responder,y=group,centrality.plotting = FALSE,results.subtitle = F,
                        bf.message = FALSE,pairwise.comparisons = FALSE,)+theme_classic2()
fisher.test(table(data2$Responder,data2$group)) # p-value = 0.2361

ggstatsplot::ggbetweenstats(data2,x=group,y=Dysfunction.x,centrality.plotting = FALSE,results.subtitle = FALSE,
                            bf.message = FALSE,pairwise.comparisons = FALSE,)+stat_compare_means()+theme_classic2()+theme(legend.position = "none") 
export::graph2ppt(file="figure/TIDE.pptx",width = 5,height = 5,append=T)


# 3.immuneAI预测免疫治疗 --------------------------------------------------------


imm_ai <- data.table::fread("../maf分组/cibersort/immucellAI/ImmuCellAI_icb_result.txt") %>% rename(sample =V1)
imm_ai <- imm_ai[imm_ai$sample %in% sample_group$sample]
imm_ai$group <- sample_group$group[match(imm_ai$sample,sample_group$sample)]
table(imm_ai$Response,imm_ai$group)

imm_ai %>% group_by(group) %>% summarise(mean_ct=mean(InfiltrationScore))

ggstatsplot::ggbarstats(imm_ai,x=Response,y=group,centrality.plotting = FALSE,results.subtitle = F,
                            bf.message = FALSE,pairwise.comparisons = FALSE,)+theme_classic2()
export::graph2ppt(file="figure/TIDE.pptx",width = 5,height = 5,append=T)

fisher.test(table(imm_ai$Response,imm_ai$group)) # p-value = 0.001229


