pkg = c('kableExtra','dplyr','tidyr','reshape','RColorBrewer','ggplot2','ggpubr','stringi','scales','tibble')
suppressPackageStartupMessages(sapply(pkg, require, character = T))
setwd('/Users/m216453/Documents/Mayo_Research/SemiSimulation/0925/sim/new/')
source('~/Documents/Mayo_project/Code/DailyCode.R')

filename = 'allC0';dir = 'sim/';type = 'allC';lineplot =F; boxplot =F;barplot =T;summary.plot = T;na.plot = F; nas = NULL;fdr = 0.05;
output = 'plot/';covariate.type ='binary'
load(paste0(dir,filename,'_res.Rdata'))
res0 <- reshape2::melt(res)
colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
res1 <- res0 %>% dplyr::select(c('depth.conf.factors','covariate.types',name0, 'nOTUs', 'measures', 'methods', 'value','iters')) %>% 
  filter(covariate.types ==covariate.type & nSams %in% c('nSam_L1','nSam_L2','nSam_L4') & nOTUs %in% c('nOTU_L1','nOTU_L3','nOTU_L5'))%>% droplevels()
levels(res1$nOTUs) <- c('OTU=50','OTU=200','OTU=500')
levels(res1$nSams) <- c('sample=50','sample=100','sample=200')
levels(res1$depth.conf.factors) <- c('None','Depth confounding')
formula = paste0('value ~ covariate.types + depth.conf.factors + nOTUs + nSams + legend + label +methods + measures')
formula.na = paste0('value ~ covariate.types + depth.conf.factors + nOTUs + nSams +methods + measures')


res1$methods = as.character(res1$methods)
res1$methods[res1$methods =='Wilcox' & res1$covariate.types =='binary'] = 'TSS+Wilcoxon'
res1$methods[res1$methods =='Wilcox' & res1$covariate.types =='continuous'] = 'TSS+Spearman'
res1$methods[res1$methods =='Rarefy' & res1$covariate.types =='binary'] = 'Rarefy+Wilcoxon'
res1$methods[res1$methods =='Rarefy' & res1$covariate.types =='continuous'] = 'Rarefy+Spearman'
res1$methods = gsub('ANCOMBC','ANCOM-BC',res1$methods)
res1$methods = gsub('glmquassi','GLM(quasipoisson)',res1$methods)
res1$methods = gsub('DESeq2.gmpr','GMPR+DESeq2',res1$methods)
res1$methods = gsub('DESeq2.Wrench','Wrench+DESeq2',res1$methods)
res1$methods = gsub('edgeR.gmpr','GMPR+edgeR',res1$methods)
res1$methods = gsub('edgeR.Wrench','Wrench+edgeR',res1$methods)
res1$methods = gsub('MSeq2.Wrench','Wrench+metagenomeSeq',res1$methods)
res1$methods = gsub('MSeq2','metagenomeSeq',res1$methods)
res1$methods = gsub('eBayW','eBay(Wilcoxon)',res1$methods)
res1$methods = gsub('eBayt','eBay(t-test)',res1$methods)
res1$methods = gsub('BBinomial','Beta-binomial',res1$methods)
res1$methods = gsub('Aldex2we','Aldex2(t-test)',res1$methods)
res1$methods = gsub('^Aldex2$','Aldex2(Wilcoxon)',res1$methods)

if(covariate.type == 'continuous'){
  sub = c('DESeq2','Wrench+DESeq2', 'edgeR', 'Wrench+edgeR')
  na.delete = c('eBay(Wilcoxon)','eBay(t-test)','RAIDA','GMPR+DESeq2','Wrench+DESeq2','GMPR+edgeR','Wrench+edgeR','Wrench+metagenomeSeq','metagenomeSeq','Aldex2(t-test)','CLRBC')
  delete = c('eBay(Wilcoxon)','eBay(t-test)','RAIDA','DESeq2','Wrench+DESeq2','edgeR','Wrench+edgeR','Wrench+metagenomeSeq','metagenomeSeq','Aldex2(t-test)','Aldex2(Wilcoxon)','CLRBC')
}else{
  sub = c('DESeq2','Wrench+DESeq2', 'edgeR', 'Wrench+edgeR','Wrench+metagenomeSeq','eBay(Wilcoxon)','Aldex2(t-test)')
  na.delete = c('eBay(t-test)','GMPR+DESeq2','Wrench+DESeq2','GMPR+edgeR','Wrench+edgeR','Wrench+metagenomeSeq','Aldex2(t-test)','CLRBC')
  delete = c('DESeq2','Wrench+DESeq2','edgeR','Wrench+edgeR','Wrench+metagenomeSeq','eBay(t-test)','Aldex2(t-test)','CLRBC')
}

res.df2 <- data_summary(data= res2, formula = formula)
measure = 'FP'
grid.formula = '. ~ nSams'
none = res.df2 %>% filter(measures == 'FP')

methods = unique(none$methods)
x1 = none %>% dplyr::select(depth.conf.factors,nOTUs, nSams, methods, value) #%>% filter(nOTUs =='OTU=50') %>% dplyr::select(-nOTUs)
x1$nOTUs = factor(x1$nOTUs, levels = c('OTU=50','OTU=200','OTU=500'))
x1$nSams = factor(x1$nSams, levels = c('sample=50','sample=100','sample=200'))
x1$depth.conf.factors = factor(x1$depth.conf.factors, levels = c('None','Depth confounding'))
x2 = x1 %>% unite('grp',c('depth.conf.factors','nOTUs','nSams'))
x3 = x2 %>% spread(c('grp'), value)
x4 = x3[,c('methods',
           "None_OTU=50_sample=50","None_OTU=50_sample=100","None_OTU=50_sample=200",
           "None_OTU=200_sample=50","None_OTU=200_sample=100","None_OTU=200_sample=200",
           "None_OTU=500_sample=50","None_OTU=500_sample=100","None_OTU=500_sample=200",
           "Depth confounding_OTU=50_sample=50","Depth confounding_OTU=50_sample=100","Depth confounding_OTU=50_sample=200",
           "Depth confounding_OTU=200_sample=50","Depth confounding_OTU=200_sample=100","Depth confounding_OTU=200_sample=200",
           "Depth confounding_OTU=500_sample=50","Depth confounding_OTU=500_sample=100","Depth confounding_OTU=500_sample=200")]
colnames(x4)[1] = 'Sample size'
colnames(x4) = gsub('.*sample=','',colnames(x4))
head(x4)
x4 = x4 %>% column_to_rownames('Sample size')
x4[x4 < 1] = '***'
x4[x4 >= 1 & x4 <2] = '**'
x4[x4 >= 2 & x4 <5] = '*'
x4[x4 >= 5] = 'x'
x41 = x4[,c(1:9)] 
x41$score = apply(x41, 1, function(x) str_count(x, "\\*") %>% sum())
x42 = x4[,c(10:18)]
x42$score = apply(x42, 1, function(x) str_count(x, "\\*") %>% sum())
x4 = merge(x41, x42, by = 0)
colnames(x4) = gsub('\\..*','',colnames(x4))
colnames(x4)[1] = 'Sample size'
x4$Overall = x41$score + x42$score 
x4 = x4[order(x4$Overall, decreasing = T),]
colnames(x4)[c(11,21:22)] = c(" "," "," ")
rownames(x4) = NULL

count.ct = function(data = x4,j = 2){
  ct = NULL
  for(i in 1:22){
    ct=c(ct,str_count(x4[i,j], "\\*") %>% sum())
  }
  ct = gsub('^3$',brewer.pal(11,'PiYG')[9],ct)
  ct = gsub('^2$',brewer.pal(8,'Set2')[6],ct)
  ct = gsub('^1$',brewer.pal(8,'RdGy')[2],ct)
  ct = gsub('^0$',brewer.pal(12,'Set3')[9],ct)
  return(ct)
} 

count.score = function(data = x4,j = 11, color1 = brewer.pal(8,'Paired')[1], color2 = 'black', value = 27){
  ct = x4[,j]
  ct[ct==value] = color1
  ct[-grep('^\\#',ct)] =color2
  return(ct)
}

kbl(x4, align = 'c') %>%
  kable_classic_2(full_width = F)%>% 
  add_header_above(c("Taxa number" = 1, "50" =3,"200" = 3,"500" =3,"Score" = 1,"50" =3,"200" = 3,"500" =3,"Score" = 1,"Overall" = 1), bold = T) %>%
  add_header_above(c("Stool-binary" = 1, "None" =9," " = 1, "Depth confounding" = 9," " = 1," " = 1), bold = T) %>%
  row_spec(0,bold=TRUE)  %>%
  column_spec(2, background= count.ct(data = x4,2),bold = T)%>%
  column_spec(3, background= count.ct(data = x4,3),bold = T)%>%
  column_spec(4, background= count.ct(data = x4,4),bold = T)%>%
  column_spec(5, background= count.ct(data = x4,5),bold = T)%>% 
  column_spec(6, background= count.ct(data = x4,6),bold = T)%>% 
  column_spec(7, background= count.ct(data = x4,7),bold = T)%>% 
  column_spec(8, background= count.ct(data = x4,8),bold = T)%>% 
  column_spec(9, background= count.ct(data = x4,9),bold = T)%>% 
  column_spec(10, background= count.ct(data = x4,10),bold = T)%>% 
  column_spec(12, background= count.ct(data = x4,12),bold = T)%>% 
  column_spec(13, background= count.ct(data = x4,13),bold = T)%>% 
  column_spec(14, background= count.ct(data = x4,14),bold = T)%>% 
  column_spec(15, background= count.ct(data = x4,15),bold = T)%>% 
  column_spec(16, background= count.ct(data = x4,16),bold = T)%>% 
  column_spec(17, background= count.ct(data = x4,17),bold = T)%>% 
  column_spec(18, background= count.ct(data = x4,18),bold = T)%>%
  column_spec(19, background= count.ct(data = x4,19),bold = T) %>% 
  column_spec(20, background= count.ct(data = x4,20),bold = T) %>%
  column_spec(11, color = count.score(data = x4, j = 11,color1 = 'white', color2 = 'black'),
              background = count.score(data = x4, j = 11,color1 = brewer.pal(8,'Dark2')[1], color2 = 'white')) %>%
  column_spec(21, color = count.score(data = x4, j = 21,color1 = 'white', color2 = 'black'),
              background = count.score(data = x4, j = 21,color1 = brewer.pal(8,'Dark2')[1], color2 = 'white')) %>%
  column_spec(22, color = count.score(data = x4, j = 22,color1 = 'white', color2 = 'black', value = 54),
              background = count.score(data = x4, j = 22,color1 = brewer.pal(8,'Dark2')[3], color2 = 'white', value = 54),bold = T)


