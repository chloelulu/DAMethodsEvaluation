pkg = c('dplyr','tidyr','reshape','RColorBrewer','ggplot2','ggpubr','stringr','scales','tibble','kableExtra','formattable')
suppressPackageStartupMessages(sapply(pkg, require, character = T))
setwd('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/')
source('~/Documents/Mayo_project/Code/DailyCode.R')

dir = 'SimulationEvaluation/';lineplot =F; boxplot =F;barplot =T;summary.plot = T;na.plot = F; nas = NULL;fdr = 0.05;
output = '../result/SimulationEvaluation/'
RES.DF2 = RES2= list()
for(filename in c('allD5','allD6','allC5','allC6')){
  for(covariate.type in c('continuous','binary')){
    if(length(grep('C',filename))==1){
      type ='allC'
    }
    if(length(grep('D',filename))==1){
      type ='allD'
    }
    
    load(paste0(dir,filename,'_res.Rdata'))
    res0 <- reshape2::melt(res)
    colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
    
    name = as.symbol('covariate.eff.means')
    name0 = 'covariate.eff.means'
    
    res0 <- res0 %>% dplyr::select(c('diff.otu.modes','covariate.types',name0, 'diff.otu.pcts','measures', 'methods', 'value','iters')) %>% filter(diff.otu.modes !='mix' &covariate.types ==covariate.type)
    res1 <- res0 %>% filter(covariate.eff.means %in% c('L2','L3','L4')) %>% droplevels()
    levels(res1$covariate.eff.means) <- c('Weak effect','Moderate effect','Strong effect')
    formula = paste0('value ~ diff.otu.modes+',name,'+ diff.otu.pcts+ methods + measures')
    formula.TPR = paste0('median ~ diff.otu.modes+',name,'+ diff.otu.pcts+  measures')
    formula.na = paste0('value ~ diff.otu.modes+',name,'+ diff.otu.pcts+ methods + measures')
    levels(res1$diff.otu.modes) = c("Abundant", "Rare")
    levels(res1$diff.otu.pcts) = c("Low density", "Medium density","High density")
    
    res1$methods = as.character(res1$methods)
    res1$methods[res1$methods =='Wilcox.gmpr' & res1$covariate.types =='binary'] = 'GMPR+Wilcoxon'
    res1$methods[res1$methods =='Wilcox.Wrench' & res1$covariate.types =='binary'] = 'Wrench+Wilcoxon'
    res1$methods[res1$methods =='Wilcox' & res1$covariate.types =='binary'] = 'TSS+Wilcoxon'
    res1$methods[res1$methods =='Wilcox.gmpr' & res1$covariate.types =='continuous'] = 'GMPR+Spearman'
    res1$methods[res1$methods =='Wilcox' & res1$covariate.types =='continuous'] = 'TSS+Spearman'
    res1$methods[res1$methods =='ttest.gmpr' & res1$covariate.types =='binary'] = 'GMPR+t-test'
    res1$methods[res1$methods =='ttest.Wrench' & res1$covariate.types =='binary'] = 'Wrench+t-test'
    res1$methods[res1$methods =='Rarefyttest' & res1$covariate.types =='binary'] = 'Rarefy+t-test'
    res1$methods[res1$methods =='ttest' & res1$covariate.types =='binary'] = 'TSS+t-test'
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
    res1$methods = gsub('^CLRBC$','LinDA',res1$methods)
    
    if(covariate.type == 'continuous'){
      sub = c('DESeq2','GMPR+DESeq2','edgeR', 'GMPR+edgeR','ANCOM-BC','DACOMP','LDM','GLM(quasipoisson)','Beta-binomial','Rarefy+Spearman','TSS+Spearman')# for continuous
      na.delete = c('eBay(Wilcoxon)','eBay(t-test)','RAIDA','GMPR+DESeq2','Wrench+DESeq2','GMPR+edgeR','Wrench+edgeR','Wrench+metagenomeSeq','metagenomeSeq','Aldex2(t-test)','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon')
      delete = c('eBay(Wilcoxon)','eBay(t-test)','RAIDA','DESeq2','Wrench+DESeq2','edgeR','Wrench+edgeR','Wrench+metagenomeSeq','metagenomeSeq','Aldex2(t-test)','Aldex2(Wilcoxon)','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon','RioNorm2','LINDA2','LinDA','GMPR+Spearman')
      sub0 =c('DESeq2','Wrench+DESeq2','GMPR+DESeq2','edgeR', 'Wrench+edgeR','GMPR+edgeR','ANCOM-BC','DACOMP','LDM','GLM(quasipoisson)','Beta-binomial','Rarefy+Spearman','TSS+Spearman')
    }else{
      na.delete = c('eBay(t-test)','GMPR+DESeq2','Wrench+DESeq2','GMPR+edgeR','Wrench+edgeR','Wrench+metagenomeSeq','Aldex2(t-test)')
      sub = c("Aldex2(Wilcoxon)","ANCOM-BC","Beta-binomial","DACOMP","eBay(Wilcoxon)",
              "GLM(quasipoisson)","GMPR+DESeq2","GMPR+edgeR","LDM","mbzinb",
              "RAIDA","Rarefy+t-test","Rarefy+Wilcoxon","TSS+t-test","TSS+Wilcoxon","metagenomeSeq",
              'DESeq2','Wrench+DESeq2','edgeR', 'Wrench+edgeR','Wrench+metagenomeSeq','eBay(t-test)','Aldex2(t-test)','RioNorm2') 
      delete = c('DESeq2','Wrench+DESeq2','edgeR','Wrench+edgeR','Wrench+metagenomeSeq','eBay(t-test)','Aldex2(t-test)','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon','LINDA2','LinDA','RioNorm2')
    }
    
    res1 = res1 %>% filter(methods %in% sub)
    
    summary <- data_summary(data= res1, formula ='value ~ methods + measures')
    
    res1_FDR  = res1[res1$measures =='FDR',]
    res1_FDR$value = res1_FDR$value/0.05
    res1_noFDR  = res1[res1$measures !='FDR',]
    res1 = rbind(res1_FDR, res1_noFDR)
    res2 <- res1[!is.na(res1$value),];dim(res2)
    res.na <- res1[(is.na(res1$value) & res1$measures =='FDR'),] %>% filter(!(methods %in% na.delete))
    res.df2 <- data_summary(data= res2, formula = formula)
    RES.DF2[[paste0(filename, covariate.type)]] = res.df2
    
    x = res.df2 %>% filter(!!as.name(name0) ==names(table(res.df2[name0]))[1] & diff.otu.pcts =='Low density' & measures =='TPR'& diff.otu.modes =='Rare')
    ord = x[order(x$value),]$methods
    res.df2$methods = factor(res.df2$methods, levels = ord)
    res2$methods = factor(res2$methods, levels = ord)
    
    letters1 = c(letters,rev(toupper(letters)))
    for(i in 1:length(unique(res.df2$methods))){
      res.df2[(res.df2$methods ==ord[i]),'label'] = letters1[i]
      res2[(res2$methods ==ord[i]),'label'] = letters1[i]
    }
    
    res.df2$legend = paste0(res.df2$label,':',res.df2$methods)
    res2$legend = paste0(res2$label,':',res2$methods)
    
    for(i in 1:length(unique(res.df2$methods))){
      res.df2[is.na(res.df2$value),'label'] = NA
    }
    
    
    if(filename %in% c(paste0(type,c(1:7,10,130,140,120, 100,60,50)))){
      res.df2 <- tidyr::complete(res.df2, diff.otu.modes,!!as.name(name),diff.otu.pcts, legend, measures)
    }
    
    for(i in 1:length(unique(res.df2$methods))){
      res.df2[is.na(res.df2$value),'label'] = NA
    }
    
    res.df2 = res.df2 %>% filter(measures %in% c('TPR','FDR','MCC','F1'))
    
    RES2[[paste0(filename, covariate.type)]] = res2
    
  }
}


# d56_continuous = rbind(RES.DF2$allD6continuous %>% mutate(include5 = 'Yes'),RES.DF2$allD5continuous %>% mutate(include5 = 'No'))
# save(d56_continuous,file = 'SimulationEvaluation/d56_continuous.Rdata')
# d56_binary = rbind(RES.DF2$allD6binary %>% mutate(include5 = 'Yes'),RES.DF2$allD5binary %>% mutate(include5 = 'No'))
# save(d56_binary,file = 'SimulationEvaluation/d56_binary.Rdata')
# 
# c56_continuous = rbind(RES.DF2$allC6continuous %>% mutate(include5 = 'Yes'),RES.DF2$allC5continuous %>% mutate(include5 = 'No'))
# save(c56_continuous,file = 'SimulationEvaluation/c56_continuous.Rdata')
# c56_binary = rbind(RES.DF2$allC6binary %>% mutate(include5 = 'Yes'),RES.DF2$allC5binary %>% mutate(include5 = 'No'))
# save(c56_binary,file = 'SimulationEvaluation/c56_binary.Rdata')
# 
# 
# 
# d56_continuous = rbind(RES2$allD6continuous %>% mutate(include5 = 'Yes'),RES2$allD5continuous %>% mutate(include5 = 'No'))
# save(d56_continuous,file = 'SimulationEvaluation/d56_continuous_res2.Rdata')
# d56_binary = rbind(RES2$allD6binary %>% mutate(include5 = 'Yes'),RES2$allD5binary %>% mutate(include5 = 'No'))
# save(d56_binary,file = 'SimulationEvaluation/d56_binary_res2.Rdata')
# 
# c56_continuous = rbind(RES2$allC6continuous %>% mutate(include5 = 'Yes'),RES2$allC5continuous %>% mutate(include5 = 'No'))
# save(c56_continuous,file = 'SimulationEvaluation/c56_continuous_res2.Rdata')
# c56_binary = rbind(RES2$allC6binary %>% mutate(include5 = 'Yes'),RES2$allC5binary %>% mutate(include5 = 'No'))
# save(c56_binary,file = 'SimulationEvaluation/c56_binary_res2.Rdata')



files = list.files('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/')
files = files[grep('_res2',files)]
files = gsub('_res2','',files)
file = files[1]
for(file in files){
  load(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/',file))
  data = get(gsub('.Rdata','',file))

  res.df2 = data %>% filter(diff.otu.modes =='Abundant') %>% dplyr::select(-diff.otu.modes)
  res.df2$include5 = as.factor(res.df2$include5)
  levels(res.df2$include5) = c('No','Yes')
  x = res.df2 %>% filter(covariate.eff.means =="Weak effect" & diff.otu.pcts =='Low density' & measures =='TPR' & include5 =='No')
  ord = x[order(x$value),]$methods
  res.df2$methods = factor(res.df2$methods, levels = ord)
  
  letters1 = c(letters,rev(toupper(letters)))
  for(i in 1:length(unique(res.df2$methods))){
    res.df2[(res.df2$methods ==ord[i]),'label'] = letters1[i]
  }
  
  res.df2$legend = paste0(res.df2$label,':',res.df2$methods)
  
  for(i in 1:length(unique(res.df2$methods))){
    res.df2[is.na(res.df2$value),'label'] = NA
  }
  
  res.df2 <- tidyr::complete(res.df2, include5,covariate.eff.means,diff.otu.pcts, legend, measures)
  
  for(i in 1:length(unique(res.df2$methods))){
    res.df2[is.na(res.df2$value),'label'] = NA
  }
  cols <- c('ZicoSeq'=brewer.pal(9,'Paired')[5], 'Omnibus' = brewer.pal(9,'PuBuGn')[7], 'mbzinb' = brewer.pal(9,'PuBuGn')[8], 
            'GLM(quasipoisson)' = brewer.pal(9,'PuRd')[3],'LinDA'=brewer.pal(9,'Set1')[8],'LINDA2'=brewer.pal(9,'Set1')[7],
            'RAIDA'=brewer.pal(11,'BrBG')[7],'RioNorm2' = brewer.pal(11,'BrBG')[8],
            'eBay(Wilcoxon)' = brewer.pal(8,'RdPu')[4],'eBay(t-test)' = brewer.pal(8,'RdPu')[6], 'ANCOM-BC' = brewer.pal(8,'Set2')[7],
            'LDM'=brewer.pal(9,'Set1')[1], 'DACOMP'=brewer.pal(8,'Dark2')[6],'Aldex2(Wilcoxon)'=brewer.pal(8,'YlGn')[3],'Aldex2(t-test)'=brewer.pal(8,'YlGn')[2],
            'Beta-binomial'=brewer.pal(11,'BrBG')[2],'DESeq2'=brewer.pal(9,'Greys')[7], 'Wrench+DESeq2'= brewer.pal(9,'Greys')[3], 'GMPR+DESeq2'=brewer.pal(9,'Greys')[5],
            'TSS+Wilcoxon'=brewer.pal(9,'Purples')[7],'Rarefy+Wilcoxon'= brewer.pal(9,'Set1')[5], 'GMPR+Wilcoxon'= brewer.pal(9,'Purples')[5], 'Wrench+Wilcoxon'= brewer.pal(9,'Purples')[3], 
            'TSS+Spearman'=brewer.pal(9,'Purples')[7],'Rarefy+Spearman'=brewer.pal(9,'Set1')[5], 'GMPR+Spearman'= brewer.pal(9,'Purples')[5], 
            'TSS+t-test'=brewer.pal(9,'GnBu')[7],'Rarefy+t-test'= brewer.pal(9,'GnBu')[5], 'GMPR+t-test'= brewer.pal(9,'GnBu')[5], 'Wrench+t-test'= brewer.pal(9,'GnBu')[3], 
            'edgeR'=brewer.pal(9,'Greens')[7], 'Wrench+edgeR'=brewer.pal(9,'Greens')[3], 'GMPR+edgeR'=brewer.pal(9,'Greens')[5], 
            'metagenomeSeq'=brewer.pal(9,'Blues')[7], 'Wrench+metagenomeSeq'=brewer.pal(9,'Blues')[5])
  
  cols1 = cols
  
  s = unique(res.df2$legend)
  for(i in 1:length(cols1)){
    if(names(cols1)[i] %in% gsub('.*:','',s)){
      names(cols1)[i] = s[which((names(cols1)[i] == gsub('.*:','',s)))]
    }else{
      cat('This legend does not exist! \n')
    }
  }
  
  
  measure = 'TPR'
  measure1 = 'FDR'
  measure2 = 'na'
  ylab = 'FDR inflation rate'
  grid.formula = 'diff.otu.pcts ~ diff.otu.modes'
  res.df2 = res.df2 %>% filter(measures %in% c('TPR','FDR'))
  
  thw =theme(axis.text.x = element_text(color="black", size =26),
             axis.text.y = element_text(color="black", size = 26),
             axis.title = element_text(color="black", size = 30),
             strip.text = element_text(size = 30),
             strip.background = element_rect(fill="white",color = "black", size = 1),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_rect(colour = "black",size= 1),
             legend.title=element_text(size=30),
             legend.text = element_text(size=30),
             plot.title = element_text(size=22))
  # dat <- data_summary(res2, formula = 'value ~ measures + methods')
  # save(dat, file = paste0(filename , '_', covariate.type, '_summary.Rdata'))
  if(length(grep('^c',file))==1){
    type.name = 'Stool'
  }else{
    type.name = 'Vaginal'
  }
  
  ## Add kable table
  head(res.df2)
  for(i in 1:3){
    cat(names(table(res.df2[,i])),'\n')
  }
  
  if(length(grep('continuous',file)) ==1){
    covariate.type ='continuous'
  }else{
    covariate.type ='binary'
  }
  
  if(covariate.type == 'continuous'){
    sub = c('DESeq2','GMPR+DESeq2','edgeR', 'GMPR+edgeR','ANCOM-BC','DACOMP','LDM','GLM(quasipoisson)','Beta-binomial','Rarefy+Spearman','TSS+Spearman')# for continuous
    na.delete = c('eBay(Wilcoxon)','eBay(t-test)','RAIDA','GMPR+DESeq2','Wrench+DESeq2','GMPR+edgeR','Wrench+edgeR','Wrench+metagenomeSeq','metagenomeSeq','Aldex2(t-test)','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon')
    delete = c('eBay(Wilcoxon)','eBay(t-test)','RAIDA','DESeq2','Wrench+DESeq2','edgeR','Wrench+edgeR','Wrench+metagenomeSeq','metagenomeSeq','Aldex2(t-test)','Aldex2(Wilcoxon)','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon','RioNorm2','LINDA2','LinDA','GMPR+Spearman')
    sub0 =c('DESeq2','GMPR+DESeq2','edgeR','GMPR+edgeR','ANCOM-BC','DACOMP','LDM','GLM(quasipoisson)','Beta-binomial','Rarefy+Spearman','TSS+Spearman')
  }else{
    na.delete = c('eBay(Wilcoxon)','GMPR+DESeq2','Wrench+DESeq2','GMPR+edgeR','Wrench+edgeR','Wrench+metagenomeSeq','Aldex2(Wilcoxon)')
    sub = c("Aldex2(Wilcoxon)",'Aldex2(t-test)',"ANCOM-BC","Beta-binomial","DACOMP",
            "GLM(quasipoisson)","LDM","mbzinb","eBay(Wilcoxon)",'eBay(t-test)',
            "RAIDA","Rarefy+t-test","Rarefy+Wilcoxon","TSS+t-test","TSS+Wilcoxon",
            "metagenomeSeq",'Wrench+metagenomeSeq','DESeq2','Wrench+DESeq2',"GMPR+DESeq2",
            'edgeR', 'Wrench+edgeR',"GMPR+edgeR",'RioNorm2') 
    delete = c('GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon','LINDA2','LinDA','RioNorm2')
  }
  
  
  name0 = 'covariate.eff.means'
  name = as.symbol(name0)
  fdr = aggregate(value~., data= res.df2 %>% filter(measures =='FDR')%>% dplyr::select(covariate.eff.means,diff.otu.pcts, include5, methods, value) , function(x)  mean(x[!is.na(x)])) %>%
    filter(!(methods %in% delete)) %>% 
    unite('grp',c('include5','diff.otu.pcts','covariate.eff.means')) %>% 
    spread(c('grp'), value)
  sum(is.na(fdr))
  fdr[is.na(fdr)] = 10
  
  
  fdr = fdr[,c('methods',
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[1]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[2]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[3]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[1]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[2]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[3]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[1]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[2]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[3]),
               
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[1]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[2]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[3]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[1]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[2]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[3]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[1]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[2]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[3])
  )]
  
  
  colnames(fdr)[1] = name0
  colnames(fdr) = gsub('.*\\_','',colnames(fdr))
  fdr = fdr %>% column_to_rownames(name0)
  fdr = apply(fdr, 2, function(x) ifelse(x<=1,"***", ifelse(x<=2, "**", ifelse(x<=3, "*", 'x'))))
  fdr1 = as.data.frame(fdr[,c(1:9)])
  fdr1$score = apply(fdr1, 1, function(x) str_count(x, "\\*") %>% sum())
  fdr2 = as.data.frame(fdr[,c(10:18)])
  fdr2$score = apply(fdr2, 1, function(x) str_count(x, "\\*") %>% sum())
  fdr = merge(fdr1, fdr2, by = 0)
  fdr$Overall = rowSums(fdr[c(11,21)])
  colnames(fdr) = gsub('\\..*','',colnames(fdr))
  colnames(fdr)[1] = name0
  fdr = fdr[order(fdr$Overall, decreasing = T),]
  colnames(fdr)[c(11,21:22)] = c(" "," "," ")
  rownames(fdr) = NULL
  
  kbl(fdr, align = 'c') %>%
    kable_classic_2(full_width = F) %>%
    add_header_above(c("signal density" = 1, "5%" =3,"10%" = 3,"20%" =3,"Score" = 1,"5%" =3,"10%" = 3,"20%" =3,"Score" = 1,"Overall" = 1), bold = T) %>%
    add_header_above(c(" " = 1, "Random" =9," " = 1, "Include top 5 taxa" = 9," " = 1," " = 1), bold = T) %>%
    row_spec(0,bold=TRUE)  %>%
    column_spec(2, background= count.ct(data = fdr,2))%>%
    column_spec(3, background= count.ct(data = fdr,3))%>%
    column_spec(4, background= count.ct(data = fdr,4))%>%
    column_spec(5, background= count.ct(data = fdr,5))%>%
    column_spec(6, background= count.ct(data = fdr,6))%>%
    column_spec(7, background= count.ct(data = fdr,7))%>%
    column_spec(8, background= count.ct(data = fdr,8))%>%
    column_spec(9, background= count.ct(data = fdr,9))%>%
    column_spec(10, background= count.ct(data = fdr,10))%>%
    column_spec(12, background= count.ct(data = fdr,12))%>%
    column_spec(13, background= count.ct(data = fdr,13))%>%
    column_spec(14, background= count.ct(data = fdr,14))%>%
    column_spec(15, background= count.ct(data = fdr,15))%>%
    column_spec(16, background= count.ct(data = fdr,16))%>%
    column_spec(17, background= count.ct(data = fdr,17))%>%
    column_spec(18, background= count.ct(data = fdr,18))%>%
    column_spec(19, background= count.ct(data = fdr,19)) %>%
    column_spec(20, background= count.ct(data = fdr,20)) %>%
    column_spec(11, color = count.score(data = fdr, j = 11,color1 = 'white', color2 = 'black'),
                background = count.score(data = fdr, j = 11,color1 = brewer.pal(8,'Dark2')[1], color2 = 'white')) %>%
    column_spec(21, color = count.score(data = fdr, j = 21,color1 = 'white', color2 = 'black'),
                background = count.score(data = fdr, j = 21,color1 = brewer.pal(8,'Dark2')[1], color2 = 'white')) %>%
    column_spec(22, color = count.score(data = fdr, j = 22,color1 = 'white', color2 = 'black', value = 54),
                background = count.score(data = fdr, j = 22,color1 = brewer.pal(8,'Dark2')[3], color2 = 'white', value = 54),bold = T) %>%
    save_kable(file = paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/SimulationEvaluation/','compositional',covariate.type,'_',type.name,'_FDR_kableTable.pdf'))
  
  
  ## TPR
  tpr = aggregate(value~., data = res.df2 %>% filter(measures =='TPR')%>% dplyr::select(name,diff.otu.pcts, include5, methods, value), function(x)  mean(x[!is.na(x)])) %>% 
    filter(!(methods %in% delete)) %>% 
    unite('grp',c('include5','diff.otu.pcts',name)) %>% 
    spread(c('grp'), value)
  
  tpr = tpr[,c('methods',
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[1]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[2]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[3]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[1]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[2]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[3]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[1]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[2]),
               paste0(levels(res.df2$include5)[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[3]),
               
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[1]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[2]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[3]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[1]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[2]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[3]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[1]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[2]),
               paste0(levels(res.df2$include5)[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[3])
  )]
  
  
  colnames(tpr)[1] = name0
  colnames(tpr) = gsub('.*\\_','',colnames(tpr))
  tpr = tpr %>% column_to_rownames(name0)
  tpr = tpr[fdr[,name0],]
  
  tpr.rank = apply(tpr, 2, rank)
  tpr1 = tpr[,c(1:9)]  
  str(tpr1)
  tpr1$`TPR score` = apply(tpr.rank[,c(1:9)], 1, function(x) sum(x)) 
  tpr1$`FDR score` = fdr[,11]
  
  tpr2 = tpr[,c(10:18)] 
  tpr2$`TPR score` = apply(tpr.rank[,c(10:18)], 1, function(x) sum(x))
  tpr2$`FDR score` = fdr[,21]
  
  tpr = merge(tpr1, tpr2, by = 0) %>% mutate(`Overall TPR` = `TPR score.y` + `TPR score.x`,`Overall FDR` = `FDR score.y` + `FDR score.x`) 
  colnames(tpr) = gsub('\\..*','',colnames(tpr))
  colnames(tpr)[1] = name0
  tpr = tpr[order(tpr$`Overall TPR`, decreasing = T),]
  colnames(tpr)[grep('FDR|TPR', colnames(tpr))] = c(" "," "," ")
  rownames(tpr) = NULL
  for(i in c(2:10,13:21)){
    tpr[,i] = format(round(tpr[,i], digits = 2), nsmall = 2)
  }
  
  
  fdr = fdr %>% column_to_rownames(name0)
  fdr = fdr[tpr[,1],] %>% rownames_to_column(name0)
  
  colnames(tpr)[1] = 'Effect size'
  
  tpr[[11]] <- color_tile("white", brewer.pal(8,'Dark2')[1])(tpr[[11]])
  tpr[[12]] <- color_tile("white", brewer.pal(8,'Dark2')[1])(tpr[[12]])
  
  tpr[[22]] <- color_tile("white", brewer.pal(8,'Dark2')[1])(tpr[[22]])
  tpr[[23]] <- color_tile("white", brewer.pal(8,'Dark2')[1])(tpr[[23]])
  
  tpr[[24]] <- color_tile("white", brewer.pal(8,'Dark2')[3])(tpr[[24]])
  tpr[[25]] <- color_tile("white", brewer.pal(8,'Dark2')[3])(tpr[[25]])
  
  
  kbl(tpr, escape = F, align = 'c') %>%
    kable_classic_2(full_width = F) %>%
    add_header_above(c("Signal density" = 1, "5%" =3,"10%" = 3,"20%" =3,"TPR Score" = 1,"FDR Score" = 1,"5%" =3,"10%" = 3,"20%" =3,"TPR Score" = 1,"FDR Score" = 1,"Overall TPR" = 1,"Overall FDR" = 1), bold = T) %>%
    add_header_above(c("Compositional" = 1, "Random" =9," " = 2, "Include top 5 taxa" = 9," " = 2," " = 2), bold = T) %>%
    row_spec(0,bold=TRUE)  %>% 
    column_spec(2, background= count.ct(data = fdr,2))%>%
    column_spec(3, background= count.ct(data = fdr,3))%>%
    column_spec(4, background= count.ct(data = fdr,4))%>%
    column_spec(5, background= count.ct(data = fdr,5))%>% 
    column_spec(6, background= count.ct(data = fdr,6))%>% 
    column_spec(7, background= count.ct(data = fdr,7))%>% 
    column_spec(8, background= count.ct(data = fdr,8))%>% 
    column_spec(9, background= count.ct(data = fdr,9))%>% 
    column_spec(10, background= count.ct(data = fdr,10))%>% 
    column_spec(13, background= count.ct(data = fdr,13-1))%>% 
    column_spec(14, background= count.ct(data = fdr,14-1))%>% 
    column_spec(15, background= count.ct(data = fdr,15-1))%>% 
    column_spec(16, background= count.ct(data = fdr,16-1))%>% 
    column_spec(17, background= count.ct(data = fdr,17-1))%>% 
    column_spec(18, background= count.ct(data = fdr,18-1))%>%
    column_spec(19, background= count.ct(data = fdr,19-1)) %>% 
    column_spec(20, background= count.ct(data = fdr,20-1))%>%
    column_spec(21, background= count.ct(data = fdr,21-1)) %>%
    column_spec(11, bold = T) %>%
    column_spec(12, bold = T)%>%
    column_spec(22, bold = T)%>%
    column_spec(23, bold = T)%>%
    column_spec(24, bold = T)%>%
    column_spec(25, bold = T) %>%
    # as_image(file = paste0(getwd(),'/plot/',covariate.type,'_',name1,'_TPR_kableTable.pdf'))
    save_kable(file = paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/SimulationEvaluation/','compositional',covariate.type,'_',type.name,'_TPR_kableTable.png'), zoom = 5)
  
  
  
  size = 6
  res.df2$include5 = gsub('Yes','Compositional: include top 5 taxa',res.df2$include5)
  res.df2$include5 = gsub('No','Compositional',res.df2$include5)
  obj1 <- ggplot(res.df2 %>% filter(measures == measure & (methods %in% sub)) %>% droplevels() %>% na.omit(), aes(x = !!as.name(name0), y = value,  fill = legend)) +
    theme_bw() +
    geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7) +
    geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
    scale_fill_manual(values = cols1) +
    scale_y_continuous(limits = c(0,1),expand = c(0.1, 0, 0, 0)) +
    facet_grid(as.formula("diff.otu.pcts ~ include5"), scales = 'free_y')+
    labs(y = 'Power', x = '', color = "", fill = '') +
    thw +
    geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5) +
    labs(title = paste0(type.name,'_',covariate.type))
  
  obj2 <- ggplot(res.df2 %>% filter(measures == measure1 & (methods %in% sub))%>% na.omit(), aes(x = !!as.name(name0), y = value,  fill = legend)) +
    theme_bw() +
    geom_bar(position = position_dodge2(width = 0.7, preserve = "single"),stat="identity", width = 0.7)+
    geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
    scale_fill_manual(values = cols1) +
    scale_y_continuous(expand = c(0.1, 0, 0, 0)) +
    facet_grid(as.formula("diff.otu.pcts ~ include5"), scales = 'free_y') +
    xlab('') +
    ylab(ylab)+
    labs(color = "", fill = '') +
    guides(fill = FALSE)+thw+
    geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5) +
    geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5)  +
    labs(title = paste0(type.name,'_',covariate.type))
  
  p =ggarrange(obj1, obj2, nrow = 2,common.legend = T)
  ggsave(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/SimulationEvaluation/','Compositional',covariate.type,'_',type.name,'.pdf'), width = 30, height = 30, dpi = 100)
  

  
  
  load(paste0("/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/",gsub('.Rdata','_res2.Rdata',file)))
  res2 = get(gsub('.Rdata','',file))
  res2 = res2 %>% filter((methods %in% sub)) %>% dplyr::select(-c(label, legend))
  me = unique(res2$methods)

  letters1 = c(letters,rev(toupper(letters)))
  for(i in 1:length(unique(res2$methods))){
    res2[(res2$methods ==me[i]),'label'] = letters1[i]
  }
  res2$legend = paste0(res2$label,':',res2$methods)
  head(res2)
  
    # dir.create('plot/SummaryPlot')
  for(k in c('TPR','FDR')){
    res3 = res2%>% filter(measures == k & !(methods %in% delete)) 
    head(res3)
    if(k =='FDR'){
      res3[res3$measures=='FDR','value'] = (res3[res3$measures=='FDR','value'])
    }
    
    if(k =='TPR'){
      res3[res3$measures=='TPR','value'] = -(res3[res3$measures=='TPR','value'])
    }
    
    res4 = res3 %>% dplyr::group_by(iters,!!as.name(name)) %>% mutate(rank=dense_rank(-value)) %>%
      dplyr::select(c(name0,'label','methods','measures','rank','value'))
    
    fac <- with(res4, reorder(label, -rank, stats::median, order = TRUE)) # higher order = higher power
    res4$label <- factor(res4$label, levels = levels(fac))
    
    p1 = ggplot(res4,aes(x = label, y = rank, fill= methods)) +
      stat_boxplot(geom = "errorbar", width = 0.7,lwd=0.2) +
      geom_boxplot(outlier.alpha = 0.5, outlier.size = 0.2, outlier.colour = 'grey') +
      theme_bw() +
      scale_fill_manual(values = cols) +
      theme(axis.text.x = element_text(color="black", size =26, vjust = 0.25, hjust = 1),
            axis.text.y = element_text(color="black", size = 26),
            axis.title = element_text(color="black", size = 26),
            strip.text = element_text(size = 26),
            strip.background = element_rect(fill="white",color = "black", size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black",size= 1),
            legend.position = 'none') +
      ylab(paste('Rank of',k)) + xlab('') +
      labs(title = paste0(type.name,'_',covariate.type))
    ggsave(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/SimulationEvaluation/','SummaryPlot/',k,'_compositional_',covariate.type,'_',type.name,'_oneplot_rank.pdf'), width = 7, height = 6, dpi = 100)
  }
}



