
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width=6.6, height=4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #绘制右边的森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
  dev.off()
}
############绘制森林图函数############

#定义独立预后分析函数
indep=function(riskFile=null,cliFile=null,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
  # cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #读取临床文件
  cli=test2
  #数据合并
  sameSample=intersect(row.names(cli),row.names(risk))
  risk=risk[sameSample,]
  cli=cli[sameSample,]
  rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
  dim(rt)
  if(anyNA(rt) == T){rt <- na.omit(rt)}
  dim(rt)
  #单因素独立预后分析
  uniTab=data.frame()
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    uniTab=rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="#0072b5")
  
  
  #多因素独立预后分析
  uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
  rt1=rt[,c("futime","fustat",as.vector(uniTab[,"id"]))]
  multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
  multiCoxSum=summary(multiCox)
  multiTab=data.frame()
  multiTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  multiTab=cbind(id=row.names(multiTab),multiTab)
  write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="#bc3c29")
}

# 注意为何用clinical2,因为独立预后可以是连续变量上

#独立预后分析
setwd("lasso500/")
indep(riskFile="trainRisk.txt",
      cliFile="clinical2.txt",
      uniOutFile="train.uniCox.txt",
      multiOutFile="train.multiCox.txt",
      uniForest="train.uniForest.pdf",
      multiForest="train.multiForest.pdf")
indep(riskFile="testRisk.txt", 
      cliFile="clinical2.txt",
      uniOutFile="test.uniCox.txt",
      multiOutFile="test.multiCox.txt",
      uniForest="test.uniForest.pdf",
      multiForest="test.multiForest.pdf")

##############临床表型亚组生存分析,不做##########
#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(survminer)
riskFile="allRisk.txt"      #风险输入文件
cliFile="clinical.txt"      #临床输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #读取临床文件

#数据合并
sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
data=cbind(futime=risk[,1],fustat=risk[,2],cli,risk=risk[,"risk"])

#对临床信息进行循环
for(i in colnames(data[,3:(ncol(data)-1)])){
  rt=data[,c("futime","fustat",i,"risk")]
  rt=rt[(rt[,i]!="unknow"),]
  colnames(rt)=c("futime","fustat","clinical","risk")
  tab=table(rt[,"clinical"])
  tab=tab[tab!=0]
  #对每个临床信息里面的分类进行循环
  for(j in names(tab)){
    rt1=rt[(rt[,"clinical"]==j),]
    tab1=table(rt1[,"risk"])
    tab1=tab1[tab1!=0]
    labels=names(tab1)
    if(length(labels)==2){
      titleName=j
      if((i=="age") | (i=="Age") | (i=="AGE")){
        titleName=paste0("age",j)
      }
      diff=survdiff(Surv(futime, fustat) ~risk,data = rt1)
      pValue=1-pchisq(diff$chisq,df=1)
      if(pValue<0.001){
        pValue="p<0.001"
      }else{
        pValue=paste0("p=",sprintf("%.03f",pValue))
      }
      fit <- survfit(Surv(futime, fustat) ~ risk, data = rt1)
      #绘制生存曲线
      surPlot=ggsurvplot(fit, 
                         data=rt1,
                         conf.int=F,
                         pval=pValue,
                         pval.size=6,
                         title=paste0("Patients with ",titleName),
                         legend.title="Risk",
                         legend.labs=labels,
                         font.legend=12,
                         xlab="Time(years)",
                         break.time.by = 1,
                         palette=c("#bc3c29", "#0072b5"),
                         risk.table=TRUE,
                         risk.table.title="",
                         risk.table.col = "strata",
                         risk.table.height=.25)
      #输出图片
      j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
      pdf(file=paste0("survival.",i,"_",j,".pdf"),onefile = FALSE,
          width = 4,        #图片的宽度
          height =5.5)        #图片的高度
      print(surPlot)
      dev.off()
    }
  }
}

#############风险图######

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")


#引用包
library(limma)
library(pheatmap)
ClusterFile="cluster.txt"     #分型结果文件
cliFile="clinical.txt"        #临床数据文件
riskFile="allRisk.txt"        #风险文件
scoreFile="scores.txt"        #肿瘤微环境文件

#读取分型的结果文件
Cluster=read.table(ClusterFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(Cluster)=c("Cluster")
rownames(Cluster)=stringr::str_sub(rownames(Cluster),1,12)
#读取临床数据文件
cli <-  data.table::fread("clinical.txt") 
cli <- as.data.frame(cli)
row.names(cli) <- cli[,1]

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取肿瘤微环境打分文件，并对数据进行整理
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
score=avereps(score)
rownames(score)=stringr::str_sub(rownames(score),1,12)

#合并数据
samSample=intersect(row.names(Cluster), row.names(cli))
Cluster=Cluster[samSample,"Cluster",drop=F]
cli=cli[samSample,,drop=F]
risk=risk[samSample,,drop=F]
score=score[samSample,,drop=F]
score[,"ImmuneScore"]=ifelse(score[,"ImmuneScore"]>median(score[,"ImmuneScore"]), "High", "Low")
data=cbind(risk, Cluster, score[,"ImmuneScore",drop=F], cli)
data=data[order(data$riskScore),,drop=F]      #根据风险打分对样品排序
Type=data[,(ncol(risk):ncol(data))]      #提取临床信息，作为热图注释文件
exp=data[,(3:(ncol(risk)-2))]       #提取lncRNA表达量
Type$Cluster=paste0("Cluster", Type$Cluster)

#对临床性状进行循环，得到显著性标记
sigVec=c("risk")
for(clinical in colnames(Type[,2:ncol(Type)])){
  data=Type[c("risk", clinical)]
  colnames(data)=c("risk", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  tableStat=table(data)
  stat=chisq.test(tableStat)
  pvalue=stat$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  sigVec=c(sigVec, paste0(clinical, Sig))
  #print(paste(clinical, pvalue, Sig, sep="\t"))
}
colnames(Type)=sigVec

#定义热图注释的颜色
colorList=list()
#Type=Type[apply(Type,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
bioCol=c("#FF0000","#0066FF","#0066FF","#FF0000","#FF9900","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
         "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
j=0
for(cli in colnames(Type[,1:ncol(Type)])){
  cliLength=length(levels(factor(Type[,cli])))
  cliCol=bioCol[(j+1):(j+cliLength)]
  j=j+cliLength
  names(cliCol)=levels(factor(Type[,cli]))
  if("unknow" %in% levels(factor(Type[,cli]))){
    cliCol["unknow"]="grey75"}
  colorList[[cli]]=cliCol
}

#热图可视化
pdf("heatmap_risk——clinical.pdf", height=5, width=9)
pheatmap(t(exp),
         annotation=Type,
         annotation_colors = colorList,
         color = colorRampPalette(c(rep("#0072b5",5), "white", rep("#bc3c29",5)))(100),
         cluster_cols =F,
         cluster_rows =F,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6)
dev.off()




######不同表型比较得分
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#引用包
library(limma)
library(ggpubr)
ClusterFile="cluster.txt"     #分型结果文件
cliFile="clinical.txt"        #临床数据文件
riskFile="allRisk.txt"        #风险文件
scoreFile="scores.txt"        #肿瘤微环境打分文件

#读取分型结果文件
Cluster=read.table(ClusterFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(Cluster)=c("Cluster")
rownames(Cluster)=stringr::str_sub(rownames(Cluster),1,12)
#读取临床数据文件
cli <-  data.table::fread("clinical.txt") 
cli <- as.data.frame(cli)
row.names(cli) <- cli[,1]

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取免疫结果文件，并对数据进行整理
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
score=avereps(score)
rownames(score)=stringr::str_sub(rownames(score),1,12)
#合并数据
samSample=intersect(row.names(Cluster), row.names(cli))
Cluster=Cluster[samSample,"Cluster",drop=F]
cli=cli[samSample,,drop=F]
risk=risk[samSample,,drop=F]

score=score[samSample,,drop=F]
score[,"ImmuneScore"]=ifelse(score[,"ImmuneScore"]>median(score[,"ImmuneScore"]), "High", "Low")
data=cbind(risk, Cluster, score[,"ImmuneScore",drop=F], cli)
rt=data[order(data$riskScore),,drop=F] 
rt=rt[,((ncol(risk)-1):ncol(rt))]
rt=rt[,-2]
rt$Cluster=paste0("Cluster", rt$Cluster)

#临床相关性分析，输出图形结果
for(clinical in colnames(rt[,2:ncol(rt)])){
  data=rt[c("riskScore", clinical)]
  colnames(data)=c("riskScore", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  #设置比较组
  group=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  #绘制箱线图
  boxplot=ggboxplot(data, x="clinical", y="riskScore", color="clinical",
                    xlab=clinical,
                    ylab="Risk score",
                    legend.title=clinical,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  #输出图片
  pdf(file=paste0(clinical, ".pdf"), width=5.5, height=5)
  print(boxplot)
  dev.off()
}



####检查点
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#引用包
library(limma)
library(ggpubr)
gene="PDCD1LG2"         #基因的标准名称
showName="PD-L2"     #图形里面显示的基因名称

#定义绘制图形的函数
load('LUAD_tpm.Rdata')
data=exprSet_tcga_mRNA
#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2", "1", group)
data=data[,group==0]

#提取目标基因表达量
data=rbind(data, gene=data[gene,])
exp=t(data[c("gene",gene),])
row.names(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\",  row.names(exp))
exp=avereps(exp)
rownames(exp)=stringr::str_sub(rownames(exp),1,12)
#读取风险数据文件
riskFile="allRisk.txt"
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(exp), row.names(risk))
exp=exp[sameSample,]
#exp[exp>quantile(exp,0.99)]=quantile(exp,0.99)
risk=risk[sameSample,]
data=cbind(as.data.frame(exp), as.data.frame(risk))

#设置比较组
data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制boxplot
boxplot=ggboxplot(data, x="risk", y="gene", color="risk",
                  xlab="",
                  ylab=paste0(showName, " expression"),
                  legend.title="",
                  palette = c("#0072b5", "#bc3c29"),
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons)

#输出图片
pdf(file='PDL2-risk.pdf', width=4, height=5.5)
print(boxplot)
dev.off()



#### 免疫相关性

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")


#引用包
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
immFile="CIBERSORT-Results.txt"       #免疫细胞浸润的结果文件
riskFile="allRisk.txt"                #风险输入文件
pFilter=0.05       #免疫细胞的过滤条件

#读取免疫结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

#删掉正常样品
group=sapply(strsplit(row.names(immune),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
immune=immune[group==0,]
row.names(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(immune))
immune=avereps(immune)

rownames(immune)=stringr::str_sub(rownames(immune),1,12)
#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#交集样品数据
sameSample=intersect(row.names(immune),row.names(risk))
immune1=immune[sameSample,]
risk1=risk[sameSample,]

#相关性检验
outTab=data.frame()
x=as.numeric(risk1[,"riskScore"])
#按免疫细胞循环
for(j in colnames(immune1)){
  y=as.numeric(immune1[,j])
  if(sd(y)>0.001){
    df1=as.data.frame(cbind(x,y))
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pValue=corT$p.value
    p1=ggplot(df1, aes(x, y)) + 
      xlab("Risk score")+ ylab(j)+
      geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
      stat_cor(method = 'spearman', aes(x =x, y =y))
    if(pValue<pFilter){
      pdf(file=paste0(j,".pdf"), width=5, height=4.6)
      print(p1)
      dev.off()
      outTab=rbind(outTab,cbind(Cell=j, pValue))
    }
  }
}
write.table(outTab,file="immuneCor.result.txt",sep="\t",row.names=F,quote=F)




