# p_load(do,tidyverse,ggplot2)
rm(list=ls())
library(dplyr)
# 载入数据 --------------------------------------------------------------------
load("../mitoGENE/outdata/LUAD_cln_clean.rds")
allrisk <- data.table::fread("lasso500/allRisk.txt")
sample_group <- allrisk %>% select(id,risk) %>% 
  rename(sample = id, group = risk) %>% 
  mutate(sample = paste0(sample,"-01A")) %>% 
  arrange(group)
head(sample_group)

# table one 制作 ------------------------------------------------------------------
library(compareGroups)
table(sample_group$group)
str(phe_f)
dt <- phe_f |> filter(sample %in% sample_group$sample) |> 
  dplyr::select( sample, age, gender,race, starts_with("smoke"),
                 starts_with("ajcc"), starts_with("stage"),
                 radiotherapy, resection_site) |> 
  dplyr::mutate(group = sample_group$group[match(sample,sample_group$sample)]) 

str(dt)

dt_factor <- dt[,-1] |> 
  transmute(
    group = factor(group,levels = c("low","high")),
    Age = as.numeric(age),
    Age_group  = as.factor(ifelse(Age > 60,">60","<=60")),
    Age_median  = as.factor(ifelse(Age > 66,">66","<=66")),
    Gender = as.factor(gender),
    Race = factor(race,levels = c("white","black","other"),labels = c("white","black","other")),
    Smoke = as.factor(smoke_group),
    Tstage = as.factor(ajcc_T),
    Nstage = as.factor(ajcc_N),
    Mstage = as.factor(ajcc_M),
    Stage_group = as.factor(stage_group),
    Stage = factor(stage2, levels = c("i","ii","iii","iv"), labels = c("I","II","III","IV")),
    Radiotherapy = as.factor(radiotherapy),
    Tumor_site = as.factor(resection_site),
  ) 



phetable2 <- descrTable(formula = group~. , # 公式，左边为分组，右边为变量
                        data = dt_factor,  # 数据集
                        method = 4, # 根据数据实际情况，自动选择统计方法
                        alpha = 0.05,  # 显著性水平
                        Q1 = 0.25, Q3 = 0.75,  # 默认输出p25和p75的分位数结果
                        show.n = T, # 显示样本量
                        show.ci = F, # 显示置信区间，默认是F
                        conf.level = 0.95, # 置信区间范围
                        type = 2, # 分类变量会显示频数和百分比
                        show.p.overall = TRUE, # 显示P值
                        digits.p = 3, # p值小数点位数
                        sd.type = 2, # 1位mean（sd），2位mean±sd
)#加载包
phetable1 <- descrTable(group ~ ., 
                        data = dt_factor, 
                        method = 2, #对数值变量。 1-强制分析为 "正态分布"；2-强制分析为 "连续非正态"；3-强制分析为 "分类"；4-正态检验后决定
                        type = 2 )
print(phetable1)
print(phetable2)

export2word(phetable2,file = "table1.docx",) #导出时注意 修改word名字







# OS和PFS ------------------------------------------------------------------
pacman::p_load(survival, survminer, export)


# 创建一个生存分析专用的dataframe 包含进去OS时间和PFS时间
dt_sur <- dt_factor |> 
  mutate(
    sample = dt$sample
  ) |> 
  inner_join(phe_f[,1:5],by = c("sample"="sample")) 
# |> 
# mutate(OS = as.factor(OS),
#        PFS = as.factor(PFS))
# 以下代码直接运行即可
# 筛选生存时间小于100 的病人做OS分析
data1 <- dt_sur %>% filter(OS.time <= 84) 

fit<-survfit(Surv(OS.time,OS)~group, data = data1) 

surv_summary(fit) #查看生存率及其标准误
surv_median(fit = fit,combine = F) # 查看中位生存时间
# strata   median    lower    upper
# 1  group=low 50.53333 40.96667       NA
# 2 group=high 37.83333 32.53333 48.46667
# 注意因为生存数据非正态，所以中位生存时间的95% 置信区间不是 对称的。
pp <- ggsurvplot(fit, data=data1, linetype = 1, palette = c("jco"),size=1,
                 surv.scale = c("percent"),pval = TRUE,
                 legend.title = 'OS', 
                 break.time.by =12,
                 xlim = c(0,84),
                 risk.table = TRUE,
                 risk.table.title = "Patients at risk",
                 ylab = "Overall Survival, %",xlab = "Months",
                 font.x = c(20,"plain","black"),font.y = c(20,"plain","black"),
                 font.tickslab = c(20,"plain","black"),
                 risk.table.fontsize = 6,font.legend =c(20,"plain","black"),
                 font.main = c(20,"plain","black"),pval.size = 8)

######risk table 的修改
pp$table <- ggpar(
  pp$table,
  font.title    = c(13, "bold.italic", "green"),  ###risk table 标题的修改
  font.subtitle = c(15, "bold", "pink"),  ###risk table小标题的修改
  font.caption  = c(11, "plain", "darkgreen"), ####插入字的修改
  font.x        = c(18, "plain", "black"), ### risk table x的修改
  font.y        = c(18, "plain", "black"),### risk table y的修改
  font.xtickslab = c(18, "plain", "black"),### risk table x 坐标轴的修改
  font.ytickslab = c(18),  ### risk table y 坐标轴的修改
  legend=c(0.8,0.88),
  censor.size=3
)
print(pp)
graph2ppt(file="figure/OS.pptx", append = T,width=7,heigh=6)

#### PFS

data1 <- dt_sur %>% filter(PFS.time <= 84) 

fit<-survfit(Surv(PFS.time,PFS)~group, data = data1) 

surv_summary(fit) #查看生存率及其标准误
surv_median(fit = fit,combine = F) # 查看中位生存时间
# strata   median    lower    upper
# 1  group=low 40.30000 34.20000 52.26667
# 2 group=high 25.73333 21.13333 32.90000
pp <- ggsurvplot(fit, data=data1, linetype = 1, palette = c("jco"),size=1,
                 surv.scale = c("percent"),pval = TRUE,
                 legend.title = 'PFS', 
                 break.time.by =12,
                 xlim = c(0,84),
                 risk.table = TRUE,
                 risk.table.title = "Patients at risk",
                 ylab = "Overall Survival, %",xlab = "Months",
                 font.x = c(20,"plain","black"),font.y = c(20,"plain","black"),
                 font.tickslab = c(20,"plain","black"),
                 risk.table.fontsize = 6,font.legend =c(20,"plain","black"),
                 font.main = c(20,"plain","black"),pval.size = 8)

######risk table 的修改
pp$table <- ggpar(
  pp$table,
  font.title    = c(13, "bold.italic", "green"),  ###risk table 标题的修改
  font.subtitle = c(15, "bold", "pink"),  ###risk table小标题的修改
  font.caption  = c(11, "plain", "darkgreen"), ####插入字的修改
  font.x        = c(18, "plain", "black"), ### risk table x的修改
  font.y        = c(18, "plain", "black"),### risk table y的修改
  font.xtickslab = c(18, "plain", "black"),### risk table x 坐标轴的修改
  font.ytickslab = c(18),  ### risk table y 坐标轴的修改
  legend=c(0.8,0.88),
  censor.size=3
)
print(pp)

graph2ppt(file="figure/OS.pptx", append = T,width=7,heigh=6)

save(dt_sur, sample_group,data1, file = "outdata/step5_sur+vars.rds")
