rm(list = ls())
library(tinyarray)
library(pacman)
p_load(limma,ggpubr,tidyverse,ggthemes)
# 1.肺癌和癌旁的差异基因 ------------------------------------------------------------
load("../mitoGENE/outdata/LUAD_TCGA_fpkm_symbol_maxgene.Rdata")
sample_group <- data.frame(pat = substring(colnames(fpkm),1,12))
sample_group$sample <- colnames(fpkm)
sample_group <- sample_group %>% 
  as.data.frame() %>% 
  mutate(group = ifelse(as.numeric(substring(sample,14,15))<10,"Tumor","Normal"))
sample_group$group %>% table()

sample_group[sample_group$group=="Normal" & endsWith(sample_group$sample,"A"),'pat'] -> np
np %>% anyDuplicated();np %>% length()

tp <- sample_group$pat[endsWith(sample_group$sample,"01A")];length(tp)


intersect(np,tp) -> peidui
anyDuplicated(peidui)

# 选择癌和癌旁的配对做差异分析
peidui %>% length();tp %>% length()

sample_group1 <- sample_group
sample_group <- sample_group %>% 
  filter( pat %in% peidui & endsWith(sample,"A")) %>% 
  arrange(group)
sample_group$group %>% table()
# Normal  Tumor 
# 56     56 
mrna <- fpkm[,sample_group$sample]
########### 设计limma比较分组矩阵以及PCA################
{
  p_load(FactoMineR,factoextra) # PCA用到的包
  sample_group$group |> table()
  
  group_list <- c(rep(0,56), rep(1, 56))
  group <- factor(group_list, 
                  levels = c(1,0), 
                  labels = c("Tumor","Normal"))
  group
  design <- model.matrix(~0+group)
  rownames(design) <- colnames(mrna)
  colnames(design) <- levels(group)
  design
  contrast.matrix <- makeContrasts(Tumor - Normal,levels = design)
  contrast.matrix
  # voom(mrna, design, plot = T) 过滤掉低表达基因以后就会有一条很平滑的曲线
  ##### PCA图
  {
    dat <-  as.data.frame(t(mrna)) # 画PCA图时要求是行名是样本名，列名时探针名，。
    # 格式要求data.frame
    dat.pca <- PCA(dat, graph = F)
    # fviz_pca_ind按样本  fviz_pca_var按基因
    fviz_pca_ind(dat.pca,
                 geom.ind = "point", # c("point", "text)2选1
                 col.ind = group, # color by groups
                 palette = c("#00AFBB", "#E7B800"),# 自定义颜色
                 addEllipses = T, # 加圆圈
                 legend.title = "Groups"# 图例名称
    ) + coord_fixed(1)
    export::graph2ppt(file="figure/PCA.pptx", append = T)
    
    # plotMDS(mrna, col = as.numeric(group))
  }
}
#######进行limma 差异分析######
{
  fit <- lmFit(mrna,design)
  fit1 <- contrasts.fit(fit, contrast.matrix)
  fit1 <- eBayes(fit1)
  qqt(fit1$t, df = fit1$df.prior+fit1$df.residual, pch = 16, cex = 0.2)
  abline(0,1) #QQ图看正态分布?
  all_diff <- topTable(fit1, 
                       adjust.method = 'fdr',
                       coef=1,
                       p.value = 1,
                       lfc <- log(1,2),
                       number = Inf,
                       sort.by = 'logFC')
  head(all_diff)
  # 所有差异加一列ID改成gene名字
  all_diff$ID <- rownames(all_diff)
  # 加一列表示logP
  all_diff$logP<- -log10(all_diff$adj.P.Val)
  #保存为csv
  write.csv(all_diff,file="outdata/fpkm_alldiff.csv")
}


######简单火山图#########

{
  
  FCvalue <- 1
  padj <- 0.05
  #将adj.P.Val<0.05,logFC>0.5的基因设置为显著上调基因
  #将adj.P.Val<0.05,logFC<0.5的基因设置为显著下调基因
  
  all_diff$group <-  "not-significant"
  all_diff$group[which((all_diff$adj.P.Val< padj) & (all_diff$logFC> FCvalue))] = "up-regulated"
  all_diff$group[which((all_diff$adj.P.Val< padj) & (all_diff$logFC< -FCvalue))] = "down-regulated"
  #查看上调和下调基因的数目
  table(all_diff$group)
  all_diff$group <- as.factor(all_diff$group)
  all_diff$logP <- -log10(all_diff$adj.P.Val)
  #火山图，加上颜色
  ggscatter(all_diff,x="logFC",y="logP",
            color = "group",
            palette = c("green","gray","red"),
            size = 1.5)+theme_base()
  #再加上辅助线
  p <- ggscatter(all_diff,x="logFC",y="logP",
                 color = "group",
                 palette = c("green","gray","red"),
                 size = 1.5)+theme_base() + 
    geom_hline(yintercept = 1.30,linetype = "dashed") + 
    geom_vline(xintercept = c(-FCvalue,FCvalue),linetype = "dashed")
  p
  # dev.size("px")
  # ggsave(p,filename = "output/KEAP1_diffgene_number_volcanoplot.pdf")      ##,width = ,height = )
}
#加上排名前十的基因
{
  all_diff$label <-  ""
  #对差异基因的p值由小到大排序，学会一个order算法！
  all_diff <- all_diff[order(all_diff$adj.P.Val),]
  #
  all_diff$X <- rownames(all_diff)
  all_diff$X <- NULL
  #高表达的基因中，选择adj.P.val最小的10个
  up.genes <- head(all_diff$ID[which(all_diff$group == "up-regulated")],10)
  #低表达的基因中，选择adj.P.val最小的10个
  down.genes <- head(all_diff$ID[which(all_diff$group == "down-regulated")],10)
  #以上两步合并加入label中
  all_diff.top10.genes <- c(as.character(up.genes),as.character(down.genes))
  all_diff$label[match(all_diff.top10.genes,all_diff$ID)] <- all_diff.top10.genes
}

#### 火山图最终成图########

p <- ggplot(data = all_diff, 
            aes(x = logFC, y = logP)) +
  geom_point(alpha=0.7, size=2, aes(color=group) ) +
  scale_color_manual(values = c("#00468BCC","gray","#ED0000CC"))+
  geom_vline(xintercept=c(-FCvalue,FCvalue),lty="dashed",col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty="dashed",col="black",lwd=0.8) +
  # ggrepel::geom_text_repel(aes(label = label), box.padding = unit(1, "lines"),
  #                          point.padding = unit(1, "lines"), show.legend = F, 
  #                          segment.color = 'black', size = 3,max.overlaps = 40)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.position = "right"
  )
p+theme_bw()
library(export)
graph2ppt(file="figure/DEanalysis.pptx", width = 7, aspectr = 1.5, append = T)
graph2pdf(file="figure/DEanalysis.pdf")
save(up.genes,down.genes,all_diff,file = "outdata/diff_analysis.Rdata")



# 3.韦恩图 -------------------------------------------------------------------

## 3.1加载dendtric cells marker -------------------------------------------------------------------
load("outdata/WGCNA_browngenes_andC1markers.Rdata") 
load("outdata/diff_analysis.Rdata") 
library(tinyarray)
library(dplyr)
# Cluster1 是NK/NKT细胞
Cluster1marker <- C1.marker %>% filter(p_val_adj < 0.05) %>% pull(gene) %>% as.character()
upgenes <- all_diff$ID[which(all_diff$group == "up-regulated")]
downgenes <- all_diff$ID[which(all_diff$group == "down-regulated")]



x=list(
  upgenes = upgenes,
  C1_marker = Cluster1marker,
  brown_module =brown.genes
)
tinyarray::draw_venn(x,name = '',color = c("#cec2f1","#89cae9","#84acaf"))
# graph2ppt(file = "output/plots/czc9象限图.pptx",append = T)

y=list(
  downgenes = downgenes,
  C1_marker = Cluster1marker,
  brown_module =brown.genes
)
tinyarray::draw_venn(y,name = '',color = c("#cec2f1","#89cae9","#84acaf"))
length(Cluster1marker)

venn_data <- data.frame(gene = c(upgenes,downgenes,brown.genes,Cluster1marker)) %>% 
  mutate(group = c(rep("upgene",924), rep("downgene",1048), rep("brown_module",1477),rep("C1_marker",500) )       )
write.table(venn_data,file = "outdata/venndata.txt",sep = "\t",row.names = F)


# 3.2clustermarker 选top500的 ---------------------------------------------------------------------
Cluster1marker <- C1.marker %>% filter(abs(avg_log2FC) > 1 & p_val_adj<0.05) %>% pull(gene) %>% as.character()
# 范围太小了，筛不出什么基因，所以不用这个策略 185个
Cluster1marker <- C1.marker %>% top_n(500,wt = avg_log2FC) %>% pull(gene) %>% as.character()

lassogene <- intersect(downgenes,brown.genes) %>% intersect(Cluster1marker);length(lassogene)
# "ADRB2"  "SPOCK2" "LAMP3"  "SFTPC"  "ITM2A"  "CCDC69" "ICAM2"  "CCND2"  "CHI3L2"
lassogene2 <- intersect(upgenes,brown.genes) %>% intersect(Cluster1marker);length(lassogene2)
# "IDH2"
lassogene_500 <- c(lassogene2,lassogene) %>% unique()
length(lassogene_500)

lassogene_185 <- c(lassogene2,lassogene) %>% unique()
length(lassogene_185)

save(lassogene_500,lassogene_185,file = "outdata/lasso_next_gene.Rdata")

# 4.LASSO ---------------------------------------------------------------------
# 代码来自day6
# 先整理数据为rt 行为样本名，列为基因名 最前面2列是生存时间和生存状态
rm(list = ls())


# 4.1 lasso建模 --------------------------------------------------------------
# install.packages("survival")
# install.packages("caret")
# install.packages("glmnet")
# install.packages("survminer")
# install.packages("timeROC")

#引用包
library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)
library(RColorBrewer)
library(ggplot2)
library(ggsci)
load("outdata/lasso_next_gene.Rdata")
load("../mitoGENE/outdata/LUAD_TCGA_fpkm_symbol_maxgene.Rdata")
load("../mitoGENE/outdata/LUAD_cln_clean.rds")
rt=read.table("outdata/uniSigExp.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
# 先整理数据为rt 行为样本名，列为基因名 最前面2列是生存时间和生存状态
rt[1:4,1:4]
fpkm[1:4,1:4]
dt <- t(fpkm)
dt <- dt[,union(lassogene_185,lassogene_500)]
index <- rownames(dt) %>% endsWith(suffix = "01A")

dt <- dt[index,] 
dt <- as.data.frame(dt)
dt$futime <- phe_f$OS.time[match(x=rownames(dt),table = phe_f$sample)]
dt$fustat <- phe_f$OS[match(x=rownames(dt),table = phe_f$sample)]
dt <- dt %>% dplyr::select(futime,fustat,everything())
dt[1:4,1:4]
rt1 <- rt
rt <- dt
library(stringr)
rownames(rt)=stringr::str_sub(rownames(rt),1,12)
rt$futime=rt$futime/12
rt$futime[rt$futime<=0]=0.003
rt %>% is.na() %>% table()
test <- na.omit(rt);dim(test);str(test)
rt <- test


# 有9个人没有生存时间，去除。

#对分组进行循环跑，找出train和test都显著的分组
# 为什么这步要设置2000循环，不给大家固定结果呢，
# 是为了让大家都能跑出略不一样的结果！！避免雷同
# 等待！
# 这一步跑完不建议重复跑

# rt_185 <- rt[,c('futime','fustat',lassogene_185)];dim(rt_185)
rt_500 <- rt[,c('futime','fustat',lassogene_500)];dim(rt_500)
# 先试试用rt185
# rt <- rt_185

rt <- rt_500
# for(i in 1:5000){
#   print(i)
#   set.seed(i)
# set.seed(59) 
 set.seed(4489) 
  #############对数据进行分组#############
  inTrain<-createDataPartition(y=rt[,3], p=0.5, list=F)
  train<-rt[inTrain,]
  test<-rt[-inTrain,]
  trainOut=cbind(id=row.names(train), train)
  testOut=cbind(id=row.names(test), test)
  
  #lasso回归
  x=as.matrix(train[,c(3:ncol(train))])
  y=data.matrix(Surv(train$futime,train$fustat))
  fit <- glmnet(x, y, family = "cox", maxit = 1000)
  cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
  # plot(fit,xvar = 'lambda')
  # plot(cvfit)

  coef <- coef(fit, s = cvfit$lambda.min)
  index <- which(coef != 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
  if(nrow(geneCoef)<2){next}
  
  #输出train组风险值
   trainFinalGeneExp=train[,lassoGene]
  # myFun=function(x){crossprod(as.numeric(x), actCoef)}
  # trainScore=apply(trainFinalGeneExp, 1, myFun)
  trainScore=predict(cvfit, newx=as.matrix(train[,c(3:ncol(train))]), s="lambda.min", type="response")
  outCol=c("futime", "fustat", lassoGene)
  risk=as.vector(ifelse(trainScore> median(trainScore), "high", "low"))
  train=cbind(train[,outCol], riskScore=as.vector(trainScore), risk)
  trainRiskOut=cbind(id=rownames(train), train)
  
  #输出test组风险值
  testFinalGeneExp=test[,lassoGene]
  #testScore=apply(testFinalGeneExp, 1, myFun)
  testScore=predict(cvfit, newx=as.matrix(test[,c(3:ncol(test))]), s="lambda.min", type="response")
  outCol=c("futime", "fustat", lassoGene)
  risk=as.vector(ifelse(testScore>median(trainScore), "high", "low"))
  test=cbind(test[,outCol], riskScore=as.vector(testScore), risk)
  testRiskOut=cbind(id=rownames(test), test)
  
  #生存差异pvalue	
  diff=survdiff(Surv(futime, fustat) ~risk, data=train)
  pValue=1-pchisq(diff$chisq, df=1)
  diffTest=survdiff(Surv(futime, fustat) ~risk, data=test)
  pValueTest=1-pchisq(diffTest$chisq, df=1)
  
  #ROC曲线下面积
  predictTime=1    #预测时间
  roc=timeROC(T=train$futime, delta=train$fustat,
              marker=trainScore, cause=1,
              weighting='aalen',
              times=c(predictTime), ROC=TRUE)
  rocTest=timeROC(T=test$futime, delta=test$fustat,
                  marker=testScore, cause=1,
                  weighting='aalen',
                  times=c(predictTime), ROC=TRUE)	
  # 训练集AUC在0.7以上，验证集0.6以上
  if((pValue<0.05) & (roc$AUC[2]>0.7) & 
     (pValueTest<0.05) & (rocTest$AUC[2]>0.65) & 
     (length(lassoGene) <= 7 ) & ( length(lassoGene) >= 3 )){
    #输出分组结果
    # write.table(trainOut,file="train.data.txt",sep="\t",quote=F,row.names=F)
    # write.table(testOut,file="test.data.txt",sep="\t",quote=F,row.names=F)
    # 4.2自定义lasso图形 -----------------------------------------------------------
    plot(fit, xvar = "lambda", label = TRUE)
    plot(cvfit)
    abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
    tidy_df <- broom::tidy(fit)
    tidy_cvdf <- broom::tidy(cvfit)
    mypalette <- c(brewer.pal(11,"BrBG"),brewer.pal(11,"Spectral"),brewer.pal(5,"Accent"))
    
p1 <- ggplot(tidy_df, aes(step, estimate, group = term,color=term)) +
      geom_line(size=1.2)+
      geom_hline(yintercept = 0)+
      ylab("Coefficients")+
      scale_color_manual(name="variable",values = mypalette)+
      theme_bw()
p1
export::graph2ppt(file = "lasso.pptx",width = 5,height = 5,append = T)

p2 <- ggplot(tidy_df, aes(lambda, estimate, group = term, color = term)) +
  geom_line(size=1.2)+
  geom_hline(yintercept = 0)+
  scale_x_log10(name = "Log Lambda")+
  ylab("Coefficients")+
  scale_color_manual(name="variable",values = mypalette)+
  theme_bw()
p2
p3 <- ggplot()+
  geom_point(data=tidy_cvdf, aes(lambda,estimate))+
  geom_errorbar(data = tidy_cvdf, aes(x=lambda,ymin=conf.low,ymax=conf.high))+
  scale_x_log10(name = "Log Lambda")+
  ylab("Coefficients")+ geom_vline(xintercept = cvfit$lambda.min,linetype = "dashed")+
  theme_bw()
p3
library(patchwork)

p2 + p3

export::graph2ppt(file = "lasso.pptx",width = 10,height = 5,append = T)
输出模型公式和风险文件
write.table(geneCoef, file="geneCoef.txt", sep="\t", quote=F, row.names=F)
write.table(trainRiskOut,file="trainRisk.txt",sep="\t",quote=F,row.names=F)
write.table(testRiskOut,file="testRisk.txt",sep="\t",quote=F,row.names=F)
#所有样品的风险值
allRiskOut=rbind(trainRiskOut, testRiskOut)
write.table(allRiskOut,file="allRisk.txt",sep="\t",quote=F,row.names=F)
    break
  }
# }




#引用包
library(survival)
library(survminer)

#自定义提取高低风险并做生存曲线的函数
bioSurvival=function(inputFile=null,outFile=null){
  #读取输入文件
  rt=read.table(inputFile, header=T, sep="\t")
  #比较高低风险组生存差异，得到显著性p值
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("#bc3c29", "#0072b5"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  print(surPlot)
  export::graph2ppt(file = "D://R language/R_Projects/cibersort/lasso.pptx",append =T)
  export::graph2pdf(file = outFile)
}
bioSurvival(inputFile="lasso185/trainRisk.txt", outFile="lasso185/trainSurv.pdf")
bioSurvival(inputFile="lasso500/trainRisk.txt", outFile="lasso500/trainSurv.pdf")
bioSurvival(inputFile="lasso185/testRisk.txt", outFile="lasso185/testSurv.pdf")
bioSurvival(inputFile="lasso500/testRisk.txt", outFile="lasso500/testSurv.pdf")



#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")


#引用包
library(survival)
library(survminer)
library(timeROC)

#自定义绘制ROC曲线函数
bioROC=function(inputFile=null, rocFile=null){
  predictTime=1    #预测时间
  #读取输入文件
  rt=read.table(inputFile, header=T, sep="\t")
  #ROC曲线
  ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
                 marker=rt$riskScore, cause=1,
                 weighting='aalen',
                 times=c(predictTime), ROC=TRUE)
  
  plot(ROC_rt, time=predictTime, col='#bc3c29', title=FALSE, lwd=2)
  legend('bottomright', cex=1.3,
         paste0('AUC=',sprintf("%.03f",ROC_rt$AUC[2])),
         col="white", lwd=1, bty = 'n')
  export::graph2pdf(file = rocFile, width=5, height=5)
  
}

bioROC(inputFile="lasso500/trainRisk.txt",rocFile="lasso500/train.ROC.pdf")
bioROC(inputFile="lasso500/testRisk.txt",rocFile="lasso500/test.ROC.pdf")



#install.packages("pheatmap")


library(pheatmap)         #引用包

bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null,heatmapFile=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
  rt=rt[order(rt$riskScore),]      #按照riskScore对样品排序
  
  #绘制风险曲线
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  lowMax=max(rt$riskScore[riskClass=="low"])
  line=rt[,"riskScore"]
  line[line>10]=10
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("#0072b5",lowLength),rep("#bc3c29",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk", "Low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  export::graph2pdf(file = riskScoreFile,width=7, height=4)
  
  #绘制生存状态图
  color=as.vector(rt$fustat)
  color[color==1]="#bc3c29"
  color[color==0]="#0072b5"
  plot(rt$futime, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("#bc3c29","#0072b5"),cex=1.2)
  abline(v=lowLength,lty=2)
  export::graph2pdf(file = survStatFile,width=7, height=4)
  #绘制风险热图
  rt1=rt[c(3:(ncol(rt)-2))]
  rt1=t(rt1)
  annotation=data.frame(type=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  p <- pheatmap(rt1, 
           annotation=annotation, 
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           show_colnames = F,
           scale="row",
           color = colorRampPalette(c(rep("#0072b5",3.5), "white", rep("#bc3c29",3.5)))(50),
           fontsize_col=3,
           fontsize=7,
           fontsize_row=8)
  export::graph2pdf(x=p, file = heatmapFile,width= 7, height=4)
}
#tarin组风险曲线
bioRiskPlot(inputFile="lasso500/trainRisk.txt",riskScoreFile="lasso500/train.riskScore.pdf",survStatFile="lasso500/train.survStat.pdf",heatmapFile="lasso500/train.heatmap.pdf")
#test组风险曲线
bioRiskPlot(inputFile="lasso500/testRisk.txt",riskScoreFile="lasso500/test.riskScore.pdf",survStatFile="lasso500/test.survStat.pdf",heatmapFile="lasso500/test.heatmap.pdf")

