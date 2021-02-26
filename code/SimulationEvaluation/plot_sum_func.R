pkg = c('dplyr','tidyr','reshape','RColorBrewer','ggplot2','ggpubr','stringr','scales','tibble','kableExtra','formattable')
suppressPackageStartupMessages(sapply(pkg, require, character = T))
setwd('/Users/m216453/Documents/Mayo_Research/SemiSimulation/0925/sim/new/')
setwd('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/')
source('~/Documents/Mayo_project/Code/DailyCode.R')

filename = 'allD0';dir = 'SimulationEvaluation/';type = 'allD';lineplot =F; boxplot =F;barplot =T;summary.plot = T;na.plot = F; nas = NULL;fdr = 0.05;
covariate.type ='binary';output = '../result/SimulationEvaluation/'


plot_sim <- function(filename = 'allC1',dir = 'SimulationEvaluation/',type = 'allC',output = '../result/SimulationEvaluation/',
                     lineplot =F, boxplot =T,barplot =T,summary.plot =T,na.plot = F, time = T,kable = T, recommendation = F,MCC =T,
                     nas = NULL,covariate.type ='binary'){
  if(filename %in% c('allC0')){
    load(paste0(dir,filename,'_melt.Rdata'))
  }else{
    load(paste0(dir,filename,'_res.Rdata'))
    res0 <- reshape2::melt(res)
    colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
  }
  # save(res0, file = 'allC0_res_melt.Rdata')
  if(filename %in% c('allC1','allC10','allC140','allC130','allC120','allC110','allC5','allC6','allC60','allD1','allD5','allD6','allC50','allD50')){
    name = as.symbol('covariate.eff.means')
    name0 = 'covariate.eff.means'
  }
  
  if(filename %in% c('allC0','allC2','allD0','allD2')){
    name = as.symbol('nSams')
    name0 = 'nSams'
  }
  
  if(filename %in% c('allC3','allD3')){
    name = as.symbol('depth.mus')
    name0 = 'depth.mus'
  }
  
  if(filename %in% c('allC4','allD4')){
    name = as.symbol('nOTUs')
    name0 = 'nOTUs'
  }
  
  if(filename %in% c('allC7','allD7')){
    name = as.symbol('depth.conf.factors')
    name0 = 'depth.conf.factors'
  }
  
  if(filename %in% c('allC0','allD0')){
    res0 <- res0 %>% dplyr::select(c('depth.conf.factors','covariate.types',name0, 'nOTUs', 'measures', 'methods', 'value','iters')) %>% filter(covariate.types ==covariate.type)
    if(filename  %in% c('allD0','allC0')){
      res1 <- res0 %>% filter(nSams %in% c('nSam_L1','nSam_L2','nSam_L4') & nOTUs %in% c('nOTU_L1','nOTU_L3','nOTU_L5')) %>% droplevels()
    }
    levels(res1$nOTUs) <- c('OTU=50','OTU=200','OTU=500')
    levels(res1$nSams) <- c('sample=50','sample=100','sample=200')
    levels(res1$depth.conf.factors) <- c('None','Depth confounding')
    formula = paste0('value ~ covariate.types + depth.conf.factors + nOTUs + nSams + +methods + measures')
    formula.na = paste0('value ~ covariate.types + depth.conf.factors + nOTUs + nSams +methods + measures')
  }else{
    res0 <- res0 %>% dplyr::select(c('diff.otu.modes','covariate.types',name0, 'diff.otu.pcts','measures', 'methods', 'value','iters')) %>% filter(diff.otu.modes !='mix' &covariate.types ==covariate.type)
    if(filename %in% c('allC1','allC10','allC120','allC130','allC140','allD1','allC5','allD5','allC6','allC60','allD6','allC50')){
      if(filename %in% c('allC50')){
        res1 <- res0 %>% filter(covariate.eff.means %in% c('L3','L4','L5')) %>% droplevels()
      }else{
        res1 <- res0 %>% filter(covariate.eff.means %in% c('L2','L3','L4')) %>% droplevels()
      }
      levels(res1$covariate.eff.means) <- c('Weak effect','Moderate effect','Strong effect')
    }
    
    if(filename %in% c('allC2','allD2')){
      res1 <- res0 %>% filter(nSams %in% c('nSam_L1','nSam_L2','nSam_L4')) %>% droplevels()
      levels(res1$nSams) <- c('sample=50','sample=100','sample=200')
    }
    
    if(filename %in% c('allC3','allD3')){
      res1 <- res0 %>% filter(depth.mus %in% c('D2','D3','D4')) %>% droplevels()
      levels(res1$depth.mus) <- c('Low depth','Medium depth','High depth')
    }
    
    if(filename %in% c('allC4','allD4')){
      res1 <- res0 %>% filter(nOTUs %in% c('nOTU_L1','nOTU_L3','nOTU_L5')) %>% droplevels()
      levels(res1$nOTUs) <- c('OTU=50','OTU=200','OTU=500')
    }
    
    if(filename %in% c('allC7','allD7')){
      res1 <- res0 %>% filter(depth.conf.factors %in% c('DL1','DL2','DL3')) %>% droplevels()
      levels(res1$depth.conf.factors) <- c('Depth confounding +','Depth confounding ++','Depth confounding +++')
    }
    
    formula = paste0('value ~ diff.otu.modes+',name,'+ diff.otu.pcts+ methods + measures')
    formula.TPR = paste0('median ~ diff.otu.modes+',name,'+ diff.otu.pcts+  measures')
    formula.na = paste0('value ~ diff.otu.modes+',name,'+ diff.otu.pcts+ methods + measures')
    levels(res1$diff.otu.modes) = c("Abundant", "Rare")
    levels(res1$diff.otu.pcts) = c("Low density", "Medium density","High density")
  }
  
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
  
  res1 = res1 %>% filter(methods %in% sub)
  
  summary <- data_summary(data= res1, formula ='value ~ methods + measures')
  
  res1_FDR  = res1[res1$measures =='FDR',]
  res1_FDR$value = res1_FDR$value/0.05
  res1_noFDR  = res1[res1$measures !='FDR',]
  res1 = rbind(res1_FDR, res1_noFDR)
  res2 <- res1[!is.na(res1$value),];dim(res2)
  res.na <- res1[(is.na(res1$value) & res1$measures =='FDR'),] %>% filter(!(methods %in% na.delete))
  if(nrow(res.na) >0){
    res.na$value[is.na(res.na$value)] = 1
    na.sum <- aggregate(as.formula(formula.na), res.na, function(x) sum(x[!is.na(x)])) %>% mutate(value = value/1000)
    na.sum <- aggregate(value ~ methods, na.sum, function(x) sum(x)/18)
    # save(na.sum, file = paste0(filename , '_', covariate.type, '_na.Rdata'))
  }else{
    na.sum = NULL
    cat('All methods in all iterations can produce result!\n')
  }
  
  cols <- c('ZicoSeq'=brewer.pal(9,'Paired')[5], 'Omnibus' = brewer.pal(9,'PuBuGn')[7], 'mbzinb' = brewer.pal(9,'PuBuGn')[8], 
            'GLM(quasipoisson)' = brewer.pal(9,'PuRd')[3],'LinDA'=brewer.pal(9,'Set1')[8],'LINDA2'=brewer.pal(9,'Set1')[7],
            'RAIDA'=brewer.pal(11,'BrBG')[7],'RioNorm2' = brewer.pal(11,'BrBG')[10],
            'eBay(Wilcoxon)' = brewer.pal(8,'RdPu')[4],'eBay(t-test)' = brewer.pal(8,'RdPu')[6], 'ANCOM-BC' = brewer.pal(8,'Set2')[7],
            'LDM'=brewer.pal(9,'Set1')[1], 'DACOMP'=brewer.pal(8,'Dark2')[6],'Aldex2(Wilcoxon)'=brewer.pal(8,'YlGn')[3],'Aldex2(t-test)'=brewer.pal(8,'YlGn')[2],
            'Beta-binomial'=brewer.pal(11,'BrBG')[2],'DESeq2'=brewer.pal(9,'Greys')[7], 'Wrench+DESeq2'= brewer.pal(9,'Greys')[3], 'GMPR+DESeq2'=brewer.pal(9,'Greys')[5],
            'TSS+Wilcoxon'=brewer.pal(9,'Purples')[7],'Rarefy+Wilcoxon'= brewer.pal(9,'Set1')[5], 'GMPR+Wilcoxon'= brewer.pal(9,'Purples')[5], 'Wrench+Wilcoxon'= brewer.pal(9,'Purples')[3], 
            'TSS+Spearman'=brewer.pal(9,'Purples')[7],'Rarefy+Spearman'=brewer.pal(9,'Set1')[5], 'GMPR+Spearman'= brewer.pal(9,'Purples')[5], 
            'TSS+t-test'=brewer.pal(9,'GnBu')[7],'Rarefy+t-test'= brewer.pal(9,'GnBu')[5], 'GMPR+t-test'= brewer.pal(9,'GnBu')[5], 'Wrench+t-test'= brewer.pal(9,'GnBu')[3], 
            'edgeR'=brewer.pal(9,'Greens')[7], 'Wrench+edgeR'=brewer.pal(9,'Greens')[3], 'GMPR+edgeR'=brewer.pal(9,'Greens')[5], 
            'metagenomeSeq'=brewer.pal(9,'Blues')[7], 'Wrench+metagenomeSeq'=brewer.pal(9,'Blues')[5])
  res.df2 <- data_summary(data= res2, formula = formula)
  # c5 = res.df2
  # c51 = c5 %>% mutate(include5 = 'No')
  
  # c6 = res.df2 %>% mutate(include5 = 'Yes')
  # c56 = rbind(c51, c6) 
  # save(c56,file = paste0('c56',covariate.type,'.Rdata'))
  
  
  # d5 = res.df2
  # d51 = d5 %>% mutate(include5 = 'No')
  
  # d6 = res.df2 %>% mutate(include5 = 'Yes')
  # d56 = rbind(d51, d6) 
  # save(d56,file = paste0('d56',covariate.type,'.Rdata'))
  
  if(filename %in% c('allC0','allD0')){
    x = res.df2 %>% filter(!!as.name(name0) ==names(table(res.df2[name0]))[1]&measures =='FP' & nOTUs =='OTU=50'& depth.conf.factors=='None')
  }else{
    x = res.df2 %>% filter(!!as.name(name0) ==names(table(res.df2[name0]))[1] & diff.otu.pcts =='Low density' & measures =='TPR'& diff.otu.modes =='Rare')
  }
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
  
  if(filename %in% paste0(type,c(0))){
    res.df2 <- tidyr::complete(res.df2, depth.conf.factors, nOTUs, nSams, legend)
  }
  
  for(i in 1:length(unique(res.df2$methods))){
    res.df2[is.na(res.df2$value),'label'] = NA
  }
  
  cols1 = cols
  
  if(barplot|lineplot|boxplot){
    s = unique(res.df2$legend)
    for(i in 1:length(cols1)){
      if(names(cols1)[i] %in% gsub('.*:','',s)){
        names(cols1)[i] = s[which((names(cols1)[i] == gsub('.*:','',s)))]
      }else{
        cat('This legend does not exist! \n')
      }
    }
  }
  
  
  if(filename %in% c(paste0(type,c(1:8,10,140,130,120,100, 60,50)))){
    measure = 'TPR'
    measure1 = 'FDR'
    measure2 = 'na'
    ylab = 'FDR inflation rate'
    grid.formula = 'diff.otu.pcts ~ diff.otu.modes'
    res.df2 = res.df2 %>% filter(measures %in% c('TPR','FDR','MCC','F1'))
  } else {
    measure = 'FP'
    measure2 = 'na'
    grid.formula = '. ~ nSams'
    res.df2 = res.df2 %>% filter(measures == 'FP')
  }
  
  
  if(!(filename %in% c('allC0','allD0'))){
    res000 <- data_summary1(data= res.df2 %>% filter(!(methods %in% delete)), formula = formula.TPR)
    colnames(res000)[ncol(res000)] = 'max_median'
    summary.TPR <- merge(res.df2, res000, by=c('diff.otu.modes',name0, 'diff.otu.pcts','measures'))
  }else{
    summary.TPR <-NULL
  }
  
  head(res.df2)
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
  if(type =='allC'){
    type.name = 'Stool'
  }else{
    type.name = 'Vaginal'
  }
  
  ## Add kable table
  # head(res.df2)
  # res2.5 = res2
  # res2.6 = res2
  # res.df2_5 = res.df2
  # res2 = rbind(res2.5 %>% mutate(include5 = 'No'), res2.6 %>% mutate(include5 = 'Yes'))
  # c56 = rbind(res.df2_5 %>% mutate(include5 = 'No'), res.df2_6 %>% mutate(include5 = 'Yes'))
  # save(c56, file = 'sim/d56_vaginal.Rdata')
  # save(res2, file = paste0(covariate.type,'_d56_vaginal_res2.Rdata'))
  # save(res2, file = paste0(covariate.type,'_c56_stool_res2.Rdata'))
  for(i in 1:3){
    cat(names(table(res.df2[,i])),'\n')
  }
  
    fdr = aggregate(value~., data= res.df2 %>% filter(measures =='FDR'& !(methods %in% delete))%>% dplyr::select(name,diff.otu.pcts, diff.otu.modes, methods, value) , function(x)  mean(x[!is.na(x)])) %>%
      unite('grp',c('diff.otu.modes','diff.otu.pcts',name)) %>% 
      spread(c('grp'), value)
    
    sum(is.na(fdr))
    fdr[is.na(fdr)] = 10
    fdr = fdr[,c('methods',
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[1]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[2]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[3]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[1]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[2]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[3]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[1]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[2]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[3]),
                 
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[1]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[2]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[3]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[1]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[2]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[3]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[1]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[2]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[3])
    )]
    
    names.recommend = colnames(fdr)[-1]
    colnames(fdr)[1] = name0
    fdr = fdr %>% column_to_rownames(name0)
    fdr = apply(fdr, 2, function(x) ifelse(x<=1,"***", ifelse(x<=2, "**", ifelse(x<=3, "*", 'x'))))
    fdr.schema1 = fdr
    
    colnames(fdr) = gsub('.*\\_','',colnames(fdr))
    fdr.schema = fdr
    fdr = as.data.frame(fdr)
    fdr$Overall = apply(fdr, 1, function(x) str_count(x, "\\*") %>% sum())
    fdr = fdr %>% rownames_to_column(name0)
    
    # fdr1 = as.data.frame(fdr[,c(1:9)])
    # fdr1$score = apply(fdr1, 1, function(x) str_count(x, "\\*") %>% sum())
    # fdr2 = as.data.frame(fdr[,c(10:18)])
    # fdr2$score = apply(fdr2, 1, function(x) str_count(x, "\\*") %>% sum())
    # fdr = merge(fdr1, fdr2, by = 0)
    # fdr$Overall = rowSums(fdr[c(11,21)])
   
    colnames(fdr) = gsub('\\..*','',colnames(fdr))
    colnames(fdr)[1] = name0
    fdr = fdr[order(fdr$Overall, decreasing = T),]
    colnames(fdr)[20] = c(" ")
    rownames(fdr) = NULL
    
    # kbl(fdr, align = 'c') %>%
    #   kable_classic_2(full_width = F) %>%
    #   add_header_above(c("signal density" = 1, "5%" =3,"10%" = 3,"20%" =3,"Score" = 1,"5%" =3,"10%" = 3,"20%" =3,"Score" = 1,"Overall" = 1), bold = T) %>%
    #   add_header_above(c(" " = 1, "Abundant" =9," " = 1, "Rare" = 9," " = 1," " = 1), bold = T) %>%
    #   row_spec(0,bold=TRUE)  %>%
    #   column_spec(2, background= count.ct(data = fdr,2))%>%
    #   column_spec(3, background= count.ct(data = fdr,3))%>%
    #   column_spec(4, background= count.ct(data = fdr,4))%>%
    #   column_spec(5, background= count.ct(data = fdr,5))%>%
    #   column_spec(6, background= count.ct(data = fdr,6))%>%
    #   column_spec(7, background= count.ct(data = fdr,7))%>%
    #   column_spec(8, background= count.ct(data = fdr,8))%>%
    #   column_spec(9, background= count.ct(data = fdr,9))%>%
    #   column_spec(10, background= count.ct(data = fdr,10))%>%
    #   column_spec(12, background= count.ct(data = fdr,12))%>%
    #   column_spec(13, background= count.ct(data = fdr,13))%>%
    #   column_spec(14, background= count.ct(data = fdr,14))%>%
    #   column_spec(15, background= count.ct(data = fdr,15))%>%
    #   column_spec(16, background= count.ct(data = fdr,16))%>%
    #   column_spec(17, background= count.ct(data = fdr,17))%>%
    #   column_spec(18, background= count.ct(data = fdr,18))%>%
    #   column_spec(19, background= count.ct(data = fdr,19)) %>%
    #   column_spec(20, background= count.ct(data = fdr,20)) %>%
    #   column_spec(11, color = count.score(data = fdr, j = 11,color1 = 'white', color2 = 'black'),
    #               background = count.score(data = fdr, j = 11,color1 = brewer.pal(8,'Dark2')[1], color2 = 'white')) %>%
    #   column_spec(21, color = count.score(data = fdr, j = 21,color1 = 'white', color2 = 'black'),
    #               background = count.score(data = fdr, j = 21,color1 = brewer.pal(8,'Dark2')[1], color2 = 'white')) %>%
    #   column_spec(22, color = count.score(data = fdr, j = 22,color1 = 'white', color2 = 'black', value = 54),
    #               background = count.score(data = fdr, j = 22,color1 = brewer.pal(8,'Dark2')[3], color2 = 'white', value = 54),bold = T) %>%
    #   save_kable(file = paste0(getwd(),'/plot/',filename,covariate.type,'_',type.name,'_FDR_kableTable.pdf'))
    
    
    ## TPR
    tpr = aggregate(value~., data = res.df2 %>% filter(measures =='TPR'& !(methods %in% delete))%>% dplyr::select(name,diff.otu.pcts, diff.otu.modes, methods, value), function(x)  mean(x[!is.na(x)])) %>% 
      unite('grp',c('diff.otu.modes','diff.otu.pcts',name)) %>% 
      spread(c('grp'), value)
    cat(as.character(tpr[,'methods']),'\n')
    cat(length(tpr[,'methods']),'\n')
    tpr = tpr[,c('methods',
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[1]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[2]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[3]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[1]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[2]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[3]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[1]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[2]),
                 paste0(levels(res.df2$diff.otu.modes)[1],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[3]),
                 
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[1]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[2]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[1],'_',levels(res.df2[[name0]])[3]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[1]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[2]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[2],'_',levels(res.df2[[name0]])[3]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[1]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[2]),
                 paste0(levels(res.df2$diff.otu.modes)[2],'_',levels(res.df2$diff.otu.pcts)[3],'_',levels(res.df2[[name0]])[3])
    )]
    
    
    colnames(tpr)[1] = name0
    colnames(tpr) = gsub('.*\\_','',colnames(tpr))
    colnames(tpr) = gsub(' effect','',colnames(tpr))
    tpr = tpr %>% column_to_rownames(name0)
    tpr = tpr[fdr[,name0],]
  
    tpr.rank = apply(tpr, 2, rank)
    ## for recommendation matrix
    fdr.schema = fdr.schema[rownames(tpr.rank),]
    tpr0 = tpr %>% as.matrix()
    recommend = matrix(NA,nrow(tpr),ncol(tpr))
    rownames(recommend) = rownames(fdr.schema)
    colnames(recommend) = gsub(' effect','',names.recommend)

    for(i in 1:ncol(tpr0)){
      # FDR < 0.05 & TPR > 0.5, we recommend it
      ord = tpr0[,i][which(fdr.schema[,i] %in% c('***'))]
      qt = quantile(ord, probs = c(0, 0.25, 0.5, 0.75, 1))[2]
      names = names(
        tpr0[,i][which(fdr.schema[,i] %in% c('***'))][names(ord[ord>= qt])]
      )
      
      # names = names(tpr0[,i][which(fdr.schema[,i] =='***')][order(-tpr0[,i][which(fdr.schema[,i] =='***')])[1:3]])
      # FDR < 0.1 & TPR > 0.75 of all avaiable methods, we recommend it
      ord = tpr0[,i][which(fdr.schema[,i] %in% c('***','**'))]
      qt = quantile(ord, probs = c(0, 0.25, 0.5, 0.75, 1))[4]
      names2 = names(
        tpr0[,i][which(fdr.schema[,i] %in% c('***','**'))][names(ord[ord>= qt])]
        )
      names2 = names2[!(names2 %in% names)]
      recommend[rownames(recommend) %in% names,i] = '***'
      recommend[rownames(recommend) %in% names2,i] = '**'
    }
    
    recommend[is.na(recommend)] = ' '
    recommend1 = recommend 
    colnames(recommend1) = gsub('.*\\_','',colnames(recommend1))

    if(recommendation){
      kbl(recommend1, escape = F, align = 'c') %>%
        kable_classic_2(full_width = F) %>%
        add_header_above(c("Signal density" = 1, "5%" =3,"10%" = 3,"20%" =3,"5%" =3,"10%" = 3,"20%" =3), bold = T) %>%
        add_header_above(c("Differential mode of taxa" = 1, "Abundant" =9,"Rare" = 9), bold = T) %>%
        column_spec(2, background= count.ct(data = recommend1,1)) %>%
        column_spec(3, background= count.ct(data = recommend1,2)) %>%
        column_spec(4, background= count.ct(data = recommend1,3)) %>%
        column_spec(5, background= count.ct(data = recommend1,4)) %>%
        column_spec(6, background= count.ct(data = recommend1,5)) %>%
        column_spec(7, background= count.ct(data = recommend1,6)) %>%
        column_spec(8, background= count.ct(data = recommend1,7)) %>%
        column_spec(9, background= count.ct(data = recommend1,8)) %>%
        column_spec(10, background= count.ct(data = recommend1,9)) %>%
        column_spec(11, background= count.ct(data = recommend1,10)) %>%
        column_spec(12, background= count.ct(data = recommend1,11)) %>%
        column_spec(13, background= count.ct(data = recommend1,12)) %>%
        column_spec(14, background= count.ct(data = recommend1,13)) %>%
        column_spec(15, background= count.ct(data = recommend1,14)) %>%
        column_spec(16, background= count.ct(data = recommend1,15)) %>%
        column_spec(17, background= count.ct(data = recommend1,16)) %>%
        column_spec(18, background= count.ct(data = recommend1,17)) %>%
        column_spec(19, background= count.ct(data = recommend1,18)) %>%
        save_kable(file = paste0(getwd(),'/plot/',filename,covariate.type,'_',type.name,'_Recommendation_kableTable.png'), zoom = 5)
      
    }

    ## MCC/F1
    if(MCC){
      MCC = aggregate(value~., data= res.df2 %>% filter(measures =='MCC'& !(methods %in% delete))%>% dplyr::select(name,diff.otu.pcts, diff.otu.modes, methods, value) , function(x)  mean(x[!is.na(x)])) %>%
        unite('grp',c('diff.otu.modes','diff.otu.pcts',name)) %>% 
        spread(c('grp'), value) %>% column_to_rownames('methods')
      F1 = aggregate(value~., data= res.df2 %>% filter(measures =='F1'& !(methods %in% delete))%>% dplyr::select(name,diff.otu.pcts, diff.otu.modes, methods, value) , function(x)  mean(x[!is.na(x)])) %>%
        unite('grp',c('diff.otu.modes','diff.otu.pcts',name)) %>% 
        spread(c('grp'), value)%>% column_to_rownames('methods')
      MCC1 = MCC[rownames(fdr.schema1),]
      F11 = F1[rownames(fdr.schema1),]
      # FDR >= 2 stars MCC/F1 top1
      recommend2 = matrix(NA,nrow(MCC),ncol(MCC))
      rownames(recommend2) = rownames(fdr.schema1)
      colnames(recommend2) = colnames(fdr.schema1)
      
      try({
        for(i in 1:ncol(MCC1)){
        # FDR < 1 & TPR > 0.5, we recommend it
          
        f.FDR = which(fdr.schema1[,i] %in% c('***','**'))
          f.MCC = MCC[,i,drop=F][f.FDR,,drop=F]
          f.MCC = f.MCC[order(-f.MCC[,1]),,drop =F]
          if(nrow(f.MCC)>1){
            max.MCC = rownames(f.MCC)[1]#[1:2]
            f.F1 = F1[,i,drop=F][f.FDR,,drop=F]
            f.F1 = f.F1[order(-f.F1[,1]),,drop =F]
            max.F1 = rownames(f.F1)[1]#[1:2]
          }else{
            max.MCC = rownames(f.MCC)[1]
            f.F1 = F1[,i,drop=F][f.FDR,,drop=F]
            f.F1 = f.F1[order(-f.F1[,1]),,drop =F]
            max.F1 = rownames(f.F1)[1]
          }
          if(!is.na(max.F1)){
            recommend2[max.F1,i] = '***'
          }
          if(!is.na(max.MCC)){
          recommend2[max.MCC,i] = '***'
          }
        
      } })
    }

      
      tpr$`TPR score` = apply(tpr.rank, 1, function(x) sum(x))
      tpr$`FDR score` = fdr[,20]
      tpr = tpr %>% rownames_to_column(name0)
      # tpr1 = tpr[,c(1:9)]  
      # str(tpr1)
      # tpr1$`TPR score` = apply(tpr.rank[,c(1:9)], 1, function(x) sum(x))
      # tpr1$`FDR score` = fdr[,11]
      # tpr2 = tpr[,c(10:18)] 
      # tpr2$`TPR score` = apply(tpr.rank[,c(10:18)], 1, function(x) sum(x))
      # tpr2$`FDR score` = fdr[,21]
      # tpr = merge(tpr1, tpr2, by = 0) %>% mutate(`Overall TPR` = `TPR score.y` + `TPR score.x`,`Overall FDR` = `FDR score.y` + `FDR score.x`) 
      
      tprscore = tpr[,c(1,20)] # extract for heatmap prepare
      
      colnames(tpr) = gsub('\\..*','',colnames(tpr))
      colnames(tpr)[1] = name0
      tpr = tpr[order(tpr$`TPR score`, decreasing = T),]
      colnames(tpr)[grep('FDR|TPR', colnames(tpr))] = c(" "," "," ")
      rownames(tpr) = NULL
      for(i in c(2:19)){
        tpr[,i] = format(round(tpr[,i], digits = 2), nsmall = 2)
      }
      
      
      fdr = fdr %>% column_to_rownames(name0)
      fdr = fdr[tpr[,1],] %>% rownames_to_column(name0)
      
      if(filename %in% c('allC1','allD1','allC5','allD5','allC6','allD6','allC10')){
        colnames(tpr)[1] = 'Effect size'
      }
      
      if(filename %in% c('allC2','allD2')){
        colnames(tpr)[1] = 'Sample size'
      }
      
      if(filename %in% c('allC3','allD3')){
        colnames(tpr)[1] = 'Depth'
      }
      
      if(filename %in% c('allC4','allD4')){
        colnames(tpr)[1] = 'Taxa number'
      }
      
      
      if(filename %in% c('allC7','allD7')){
        colnames(tpr)[1] = 'Depth confounding'
      }
      # for(i in c(2:10,12:20)){
      #   tpr[[i]] <- ifelse(
      #     fdr[[i]] =='***',
      #     cell_spec(tpr[[i]], color = 'black', background=brewer.pal(11,'PiYG')[9], align = 'center'),
      #     ifelse(
      #       fdr[[i]] =='**',
      #       cell_spec(tpr[[i]], color = 'black', background= brewer.pal(8,'Set2')[6], align = 'center'),
      #       ifelse(
      #         fdr[[i]] =='*',
      #         cell_spec(tpr[[i]], color = 'black', background= brewer.pal(8,'RdGy')[2], align = 'center'),
      #         cell_spec(tpr[[i]], color = 'black', background= brewer.pal(12,'Set3')[9], align = 'center'))))
      # }
      
      tpr[[20]] <- color_tile("white", brewer.pal(8,'Dark2')[1])(tpr[[20]])
      tpr[[21]] <- color_tile("white", brewer.pal(8,'Dark2')[1])(tpr[[21]])
      
      # tpr[[22]] <- color_tile("white", brewer.pal(8,'Dark2')[1])(tpr[[22]])
      # tpr[[23]] <- color_tile("white", brewer.pal(8,'Dark2')[1])(tpr[[23]])
      # 
      # tpr[[24]] <- color_tile("white", brewer.pal(8,'Dark2')[3])(tpr[[24]])
      # tpr[[25]] <- color_tile("white", brewer.pal(8,'Dark2')[3])(tpr[[25]])
      
      if(kable){
        
        kbl(tpr, escape = F, align = 'c') %>%
          kable_classic_2(full_width = F) %>%
          add_header_above(c("Signal density" = 1, "5%" =3,"10%" = 3,"20%" =3,"5%" =3,"10%" = 3,"20%" =3,"TPR Score" = 1,"FDR Score" = 1), bold = T) %>%
          add_header_above(c("Differential mode of taxa" = 1, "Abundant" =9,"Rare" = 9," " = 2), bold = T) %>%
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
          column_spec(11, background= count.ct(data = fdr,11))%>% 
          column_spec(12, background= count.ct(data = fdr,12))%>% 
          column_spec(13, background= count.ct(data = fdr,13))%>% 
          column_spec(14, background= count.ct(data = fdr,14))%>% 
          column_spec(15, background= count.ct(data = fdr,15))%>% 
          column_spec(16, background= count.ct(data = fdr,16))%>% 
          column_spec(17, background= count.ct(data = fdr,17))%>% 
          column_spec(18, background= count.ct(data = fdr,18))%>%
          column_spec(19, background= count.ct(data = fdr,19)) %>% 
          # column_spec(20, background= count.ct(data = fdr,20-1))%>%
          # column_spec(21, background= count.ct(data = fdr,21-1)) %>%
          column_spec(20, bold = T) %>%
          column_spec(21, bold = T) %>%
          # column_spec(22, bold = T)%>%
          # column_spec(23, bold = T)%>%
          # column_spec(24, bold = T)%>%
          # column_spec(25, bold = T)%>%
          # as_image(file = paste0(getwd(),'/plot/',covariate.type,'_',name1,'_TPR_kableTable.pdf'))
          save_kable(file = paste0(output,filename,covariate.type,'_',type.name,'_TPR_kableTable.png'), zoom = 5)
        
      }
      ## For kable of false discover FPs
      

  
  if(summary.plot ==T & !(filename %in% c('allC0','allD0'))){
    # dir.create('plot/SummaryPlot')
    for(k in c('TPR','FDR')){
      res3 = res2 %>% separate(iters, c('iter','m'),sep = '-') %>% dplyr::select(-m) %>% filter(measures == k & (methods %in% sub)) 
      if(k =='FDR'){
        res3[res3$measures=='FDR','value'] = (res3[res3$measures=='FDR','value'])
      }
      
      if(k =='TPR'){
        res3[res3$measures=='TPR','value'] = -(res3[res3$measures=='TPR','value'])
      }
      
      res4 = res3 %>% dplyr::group_by(iter,!!as.name(name)) %>% mutate(rank=dense_rank(-value)) %>%
        dplyr::select(c(name0,'label','methods','measures','rank','value'))
      
      fac <- with(res4, reorder(label, -rank, stats::median, order = TRUE)) # higher order = higher power
      res4$label <- factor(res4$label, levels = levels(fac))
      
      p1 = ggplot(res4,aes(x = label, y = rank, fill= methods)) +
        stat_boxplot(geom = "errorbar", width = 0.7,lwd=0.2) +
        geom_boxplot(outlier.alpha = 0.5, outlier.size = 0.2, outlier.colour = 'grey') +
        theme_bw() +
        scale_fill_manual(values = cols) +
        theme(axis.text.x = element_text(color="black", size =26),
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
      ggsave(paste0(output,'SummaryPlot/',k,'_',filename,'_',covariate.type,'_oneplot_rank.pdf'), width = 8, height = 3, dpi = 100)
    }
  }
  
  
  if('na' %in% res1$measures){
    na.df <- res1 %>% filter(measures == measure2)
    na.df[is.na(na.df$value),'value'] <- 0
    na.df$value <- na.df$value *100
    na.df1 = data_summary(na.df, formula = 'value~methods')
    na.df1 <- na.df1[order(-na.df1$value),]
    ord <- within(na.df1, methods <- factor(methods, levels=na.df1$methods))
    if(na.plot){
      na.plt <- ggplot(ord, aes(x = methods, y = value, fill = methods)) +
        theme_bw() +
        geom_bar(stat="identity",position = position_dodge2(width = 0.9, preserve = "single"),  width = 0.9) +
        geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.9, size = 0.2,
                      position=position_dodge2(.9, preserve = "single")) +
        scale_fill_manual(values = cols) +
        ylab('# NAs/ # tests, %') + xlab('') +
        theme(axis.text.x = element_text(color="black", size =26, angle = 90, vjust = 0.95, hjust =1),
              axis.text.y = element_text(color="black", size = 26),
              axis.title = element_text(color="black", size = 30),
              strip.text = element_text(size = 30),
              strip.background = element_rect(fill="white",color = "black", size = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black",size= 1),
              legend.title=element_text(size=30),legend.position = 'none',
              legend.text = element_text(size=30)) +
        scale_y_continuous(trans = sqrt_trans(),
                           breaks = trans_breaks("sqrt", function(x) x^2),
                           labels = trans_format("sqrt", math_format(.x^2)))
      ggsave(paste0(output,filename,covariate.type,'_','V35_PercentOfQvalNA.pdf'), width = 12, height = 10, dpi = 100)
    }
  }
  
  
  if(filename == paste0(type,'0')){
    plt1 <- ggplot(res.df2 %>% filter(measures == measure & (methods %in% sub) & covariate.types ==covariate.type) %>% droplevels(), aes(x = nOTUs, y = value,  fill = legend)) +
      theme_bw() +
      geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7) +
      geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
      scale_fill_manual(values = cols1) +
      geom_text(aes(y=0,label=label), position=position_dodge(width=0.7),size=6,vjust = 1) +
      facet_grid(nSams~depth.conf.factors, scales = 'free_y',space = "free_x")+
      ylab('FPs') + xlab('') +
      labs(color = "Methods", title = paste0('Default: no differential taxa (',covariate.type,')'), fill = '')+
      thw +theme(legend.position = 'top') +
      scale_y_continuous(trans = sqrt_trans(),
                         breaks = trans_breaks("sqrt", function(x) x^2),
                         labels = trans_format("sqrt", math_format(.x^2)))+
      geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size = 0.5) +
      labs(title = paste0(type.name,'_',covariate.type))

    ggsave(paste0(output,filename,covariate.type,'_','barplot.pdf'),  width = 26, height = 20, dpi = 100)
    
  }else{
    # Pwer plot: lineplot(only for power), boxplot, barplot
    if(lineplot){
      size = 1.5
      line1 <- ggplot(res.df2 %>% filter(measures == measure & !(methods %in% delete)&covariate.types ==covariate.type), aes(x = !!as.name(name0), y = value, colour = methods, group= methods)) +
        theme_bw() +
        geom_errorbar(aes(ymin=ymin, ymax=ymax,width=1), size = size,  position=position_dodge(0) ) +
        geom_line(position=position_dodge(0), size = size) +
        geom_point(position=position_dodge(0),size=size, shape=21, fill="white")+
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        scale_color_manual(values = cols) +
        scale_y_continuous(limits = c(0,1),expand = c(0.1, 0, 0, 0)) + 
        labs(y = 'TPR', x = stri_trans_totitle(name), color = "Methods") +
        thw
      line11 <- ggplot(res.df2 %>% filter(measures == measure & (methods %in% sub)&covariate.types ==covariate.type), aes(x = !!as.name(name0), y = value, colour = methods,group= methods)) +
        theme_bw() +
        geom_errorbar(aes(ymin=ymin, ymax=ymax,width=1), size = size,  position=position_dodge(0) ) +
        geom_line(position=position_dodge(0), size = size) +
        geom_point(position=position_dodge(0),size=size, shape=21, fill="white")+
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        scale_color_manual(values = cols) +
        scale_y_continuous(limits = c(0,1),expand = c(0.1, 0, 0, 0)) + 
        labs(y = 'TPR', x = stri_trans_totitle(name), color = "Methods") +
        thw
      
      line2 <- ggplot(res.df2 %>% filter(measures == measure1 & !(methods %in% delete)&covariate.types ==covariate.type), aes(x = !!as.name(name0), y = value, colour = methods, group= methods)) +
        theme_bw() +
        geom_errorbar(aes(ymin=ymin, ymax=ymax,width=1), size = size,  position=position_dodge(0) ) +
        geom_line(position=position_dodge(0), size = size) +
        geom_point(position=position_dodge(0),size=size, shape=21, fill="white")+
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        scale_color_manual(values = cols) +
        scale_y_continuous(limits = c(0,1),expand = c(0.1, 0, 0, 0)) + 
        labs(y = 'TPR', x = stri_trans_totitle(name), color = "Methods") +
        thw
      line22 <- ggplot(res.df2 %>% filter(measures == measure1 & (methods %in% sub)&covariate.types ==covariate.type), aes(x = !!as.name(name0), y = value, colour = methods,group= methods)) +
        theme_bw() +
        geom_errorbar(aes(ymin=ymin, ymax=ymax,width=1), size = size,  position=position_dodge(0) ) +
        geom_line(position=position_dodge(0), size = size) +
        geom_point(position=position_dodge(0),size=size, shape=21, fill="white")+
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        scale_color_manual(values = cols) +
        scale_y_continuous(limits = c(0,1),expand = c(0.1, 0, 0, 0)) + 
        labs(y = 'TPR', x = stri_trans_totitle(name), color = "Methods") +
        thw
      p =ggarrange(line1, line2, nrow = 2,common.legend = T)
      ggsave(paste0(output,filename,covariate.type,'_','V35_line1.pdf'), width = 30, height = 40, dpi = 100)
      
      p1 =ggarrange(line11, line22, nrow = 2,common.legend = T)
      ggsave(paste0(output,filename,covariate.type,'_','V35_line2.pdf'), width = 30, height = 40, dpi = 100)
      
    }
    
    size = 6
    if(barplot){
      obj1 <- ggplot(res.df2 %>% filter(measures == measure & !(methods %in% delete)) %>% droplevels(), aes(x = !!as.name(name0), y = value,  fill = legend)) +
        theme_bw() +
        geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7) +
        geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
        scale_fill_manual(values = cols1) +
        scale_y_continuous(limits = c(0,1),expand = c(0.1, 0, 0, 0)) +
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        labs(y = 'TPR', x = '', color = "", fill = '') +
        thw +
        geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5) +
        labs(title = paste0(type.name,'_',covariate.type))

      obj11 <- ggplot(res.df2 %>% filter(measures == measure & (methods %in% sub)) %>% droplevels(), aes(x = !!as.name(name0), y = value,  fill = legend)) +
        theme_bw() +
        geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7)+
        geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
        scale_y_continuous(limits = c(0,1), expand = c(0.1, 0, 0, 0)) +
        scale_fill_manual(values = cols1) +
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        ylab('TPR') +
        xlab('') +
        labs(color = "", fill = '') +
        thw +
        geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5)  +
        labs(title = paste0(type.name,'_',covariate.type))

      obj2 <- ggplot(res.df2 %>% filter(measures == measure1 & !(methods %in% delete)), aes(x = !!as.name(name0), y = value,  fill = legend)) +
        theme_bw() +
        geom_bar(position = position_dodge2(width = 0.7, preserve = "single"),stat="identity", width = 0.7)+
        geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
        scale_fill_manual(values = cols1) +
        scale_y_continuous(expand = c(0.1, 0, 0, 0)) +
        facet_grid(as.formula(grid.formula), scales = 'free_y') +
        xlab('') +
        ylab(ylab)+
        labs(color = "", fill = '') +
        guides(fill = FALSE)+thw+
        geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5) +
        geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5)  +
        labs(title = paste0(type.name,'_',covariate.type))

      obj22 <- ggplot(res.df2 %>% filter(measures == measure1 & (methods %in% sub)), aes(x = !!as.name(name0), y = value,  fill = legend)) +
        theme_bw() +
        geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7)+
        geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
        scale_fill_manual(values = cols1) +
        scale_y_continuous(expand = c(0.1, 0, 0, 0)) +
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        ylab(ylab) +
        xlab('') +
        labs(color = "", fill = '') +
        guides(fill = FALSE)+thw+
        geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5) +
        geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5)  +
        labs(title = paste0(type.name,'_',covariate.type))
      p =ggarrange(obj1, obj2, nrow = 2,common.legend = T)
      ggsave(paste0(output,filename,covariate.type,'_','V35_part1.pdf'), width = 30, height = 30, dpi = 100)
      p =ggarrange(obj11, obj22, nrow = 2,common.legend = T)
      ggsave(paste0(output,filename,covariate.type,'_','V35_part2.pdf'),  width = 30, height = 30, dpi = 100)
      
      tpr.plt <- ggplot(res.df2 %>% filter(measures == measure & (methods %in% sub)) %>% droplevels(), aes(x = !!as.name(name0), y = value,  fill = legend)) +
        theme_bw() +
        geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7) +
        geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
        scale_fill_manual(values = cols1) +
        scale_y_continuous(limits = c(0,1),expand = c(0.1, 0, 0, 0)) +
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        labs(y = 'TPR', x = '', color = "", fill = '') +
        thw + 
        geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5) +
        labs(title = paste0(type.name,'_',covariate.type))

      fdr.plt <- ggplot(res.df2 %>% filter(measures == measure1 & (methods %in% sub)), aes(x = !!as.name(name0), y = value,  fill = legend)) +
        theme_bw() +
        geom_bar(position = position_dodge2(width = 0.7, preserve = "single"),stat="identity", width = 0.7)+
        geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
        scale_fill_manual(values = cols1) +
        scale_y_continuous(expand = c(0.1, 0, 0, 0)) +
        facet_grid(as.formula(grid.formula), scales = 'free_y') +
        xlab('') +
        ylab(ylab)+
        labs(color = "", fill = '') +
        guides(fill = FALSE)+thw+
        geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5) +
        geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5)  +
        labs(title = paste0(type.name,'_',covariate.type))
      
      p =ggarrange(tpr.plt, fdr.plt, nrow = 2,common.legend = T)
      ggsave(paste0(output,filename,covariate.type,'.pdf'),  width = 35, height = 30, dpi = 100)
      
      
    }
    
    if(boxplot){
      box1 <- ggplot(res2 %>% filter(measures == measure & !(methods %in% delete)&covariate.types ==covariate.type) %>% droplevels(), aes(x = !!as.name(name0), y = value,  fill = legend)) +
        theme_bw() +
        stat_boxplot(geom = "errorbar", width = 0.7,lwd=0.1) +
        geom_boxplot(width=0.7,outlier.size = 0, outlier.colour = 'grey30',lwd=0.1) +
        # geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, jitter.height = 0.05, jitter.width = 0.075, errorbar.draw = TRUE)+
        # geom_half_boxplot(position = position_dodge(width = 0.9),nudge = 0.05, outlier.color = NA) +
        # geom_half_point(transformation = position_quasirandom(width = .9, groupOnX = TRUE)) +
        scale_fill_manual(values = cols1) +
        scale_color_manual(values = cols1) +
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        scale_y_continuous(expand = c(0.1, 0, 0, 0)) +
        ylab('TPR') +
        xlab('') +
        labs(fill = '')+thw +
        geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5)

      box11 <- ggplot(res2 %>% filter(measures == measure & (methods %in% sub)&covariate.types ==covariate.type) %>% droplevels(), aes(x = !!as.name(name0), y = value,  fill = legend)) +
        theme_bw() +
        stat_boxplot(geom = "errorbar", width = 0.7,lwd=0.1) +
        geom_boxplot(width=0.7,outlier.size = 0.1, outlier.colour = 'grey30',lwd=0.1) +
        scale_fill_manual(values = cols1) +
        scale_color_manual(values = cols1) +
        scale_y_continuous(expand = c(0.1, 0, 0, 0)) +
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        ylab('TPR') +
        xlab('') +
        labs(fill = '')+thw +
        geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5)

      box2 <- ggplot(res2 %>% filter(measures == measure1 & !(methods %in% delete)&covariate.types ==covariate.type), aes(x = !!as.name(name0), y = value,  fill = legend)) +
        theme_bw() +
        stat_boxplot(geom = "errorbar", width = 0.7,lwd=0.1) +
        geom_boxplot(width=0.7,outlier.size = 0.1, outlier.colour = 'grey30',lwd=0.1) +
        scale_fill_manual(values = cols1) +
        scale_color_manual(values = cols1) +
        scale_y_continuous(expand = c(0.1, 0, 0, 0)) +
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        ylab(ylab) +
        xlab('') +
        labs(fill = '')+thw +
        geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5) +
        geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5)
      box22 <- ggplot(res2 %>% filter(measures == measure1 & (methods %in% sub)&covariate.types ==covariate.type), aes(x = !!as.name(name0), y = value,  fill = legend)) +
        theme_bw() +
        stat_boxplot(geom = "errorbar", width = 0.7,lwd=0.1) +
        geom_boxplot(width=0.7,outlier.size = 0.1, outlier.colour = 'grey30',lwd=0.1) +
        scale_fill_manual(values = cols1) +
        scale_color_manual(values = cols1) +
        scale_y_continuous(expand = c(0.1, 0, 0, 0)) +
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        ylab(ylab) +
        xlab('') +
        labs(fill = '')+thw +
        geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5) +
        geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5)
      p =ggarrange(box1, box2, nrow = 2,common.legend = T)
      ggsave(paste0(output,filename,covariate.type,'_','V35_part1_box.pdf'),  width = 30, height = 30, dpi = 100)
      p =ggarrange(box11, box22, nrow = 2,common.legend = T)
      ggsave(paste0(output,filename,covariate.type,'_','V35_part2_box.pdf'), width = 30, height = 30, dpi = 100)
      
      box1 <- ggplot(res2 %>% filter(measures == measure & (methods %in% sub)&covariate.types ==covariate.type) %>% droplevels(), aes(x = !!as.name(name0), y = value,  fill = legend)) +
        theme_bw() +
        stat_boxplot(geom = "errorbar", width = 0.7,lwd=0.1) +
        geom_boxplot(width=0.7,outlier.size = 0, outlier.colour = 'grey30',lwd=0.1) +
        # geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, jitter.height = 0.05, jitter.width = 0.075, errorbar.draw = TRUE)+
        # geom_half_boxplot(position = position_dodge(width = 0.9),nudge = 0.05, outlier.color = NA) +
        # geom_half_point(transformation = position_quasirandom(width = .9, groupOnX = TRUE)) +
        scale_fill_manual(values = cols1) +
        scale_color_manual(values = cols1) +
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        scale_y_continuous(expand = c(0.1, 0, 0, 0)) + 
        ylab('TPR') +
        xlab('') +
        labs(fill = '')+thw +
        geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5) 
      
      box2 <- ggplot(res2 %>% filter(measures == measure1 & (methods %in% sub)&covariate.types ==covariate.type), aes(x = !!as.name(name0), y = value,  fill = legend)) +
        theme_bw() +
        stat_boxplot(geom = "errorbar", width = 0.7,lwd=0.1) +
        geom_boxplot(width=0.7,outlier.size = 0.1, outlier.colour = 'grey30',lwd=0.1) +
        scale_fill_manual(values = cols1) +
        scale_color_manual(values = cols1) +
        scale_y_continuous(expand = c(0.1, 0, 0, 0)) + 
        facet_grid(as.formula(grid.formula), scales = 'free_y')+
        ylab(ylab) +
        xlab('') +
        labs(fill = '')+thw +
        geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5) +
        geom_text(aes(y=0, label=label), position=position_dodge(width=0.7),size=size, vjust = 1.5) 
      p =ggarrange(box1, box2, nrow = 2,common.legend = T)
      ggsave(paste0(output,filename,covariate.type,'_boxplot.pdf'),  width = 35, height = 30, dpi = 100)
      
    }
  }
  return(list(na.sum = na.df1, 
              summary = summary, # does include the max(median TPR)) for heatmap score FDR and FPR
              summary.TPR = summary.TPR, # include the max(median TPR)) for heatmap score TPR
              res2 = res2,tprscore =tprscore,res.df2 = res.df2, tpr.rank = tpr.rank,
              recommend = recommend,recommend2 = recommend2
              ))
}



c10c = plot_sim(filename = 'allC10',covariate.type = 'continuous',type ='allC',MCC = F)
c10b = plot_sim(filename = 'allC10',covariate.type = 'binary',type ='allC',MCC = F) # I forgot to calculate MCC in sample size =10, and we do not need this result in manuscript

D50c = plot_sim(filename = 'allD50',covariate.type = 'continuous',type ='allD',MCC = T)
D50b = plot_sim(filename = 'allD50',covariate.type = 'binary',type ='allD',MCC = T)

c50c = plot_sim(filename = 'allC50',covariate.type = 'continuous',type ='allC',MCC = F)
c50b = plot_sim(filename = 'allC50',covariate.type = 'binary',type ='allC',MCC = F)


c140c = plot_sim(filename = 'allC140',covariate.type = 'continuous',type ='allC')
c140b = plot_sim(filename = 'allC140',covariate.type = 'binary',type ='allC')

c130c = plot_sim(filename = 'allC130',covariate.type = 'continuous',type ='allC')
c130b = plot_sim(filename = 'allC130',covariate.type = 'binary',type ='allC')

c120c = plot_sim(filename = 'allC120',covariate.type = 'continuous',type ='allC')
c120b = plot_sim(filename = 'allC120',covariate.type = 'binary',type ='allC',)

c60c = plot_sim(filename = 'allC60',covariate.type = 'continuous',type ='allC')
c60b = plot_sim(filename = 'allC60',covariate.type = 'binary',type ='allC')

## continuous covariate
c0c = plot_sim(filename = 'allC0',covariate.type = 'continuous',type ='allC')
c0b = plot_sim(filename = 'allC0',covariate.type = 'binary',type ='allC')
c1c = plot_sim(filename = 'allC1',covariate.type = 'continuous',type ='allC')
c2c = plot_sim(filename = 'allC2',covariate.type = 'continuous',type ='allC')
c3c = plot_sim(filename = 'allC3',covariate.type = 'continuous',type ='allC')
c4c = plot_sim(filename = 'allC4',covariate.type = 'continuous',type ='allC')
c5c = plot_sim(filename = 'allC5',covariate.type = 'continuous',type ='allC')
c6c = plot_sim(filename = 'allC6',covariate.type = 'continuous',type ='allC')
c7c = plot_sim(filename = 'allC7',covariate.type = 'continuous',type ='allC')

d0c = plot_sim(filename = 'allD0',covariate.type = 'continuous',type ='allD')
d0b = plot_sim(filename = 'allD0',covariate.type = 'binary',type ='allD')
d1c = plot_sim(filename = 'allD1',covariate.type = 'continuous',type ='allD')
d2c = plot_sim(filename = 'allD2',covariate.type = 'continuous',type ='allD')
d3c = plot_sim(filename = 'allD3',covariate.type = 'continuous',type ='allD')
d4c = plot_sim(filename = 'allD4',covariate.type = 'continuous',type ='allD')
d5c = plot_sim(filename = 'allD5',covariate.type = 'continuous',type ='allD')
d6c = plot_sim(filename = 'allD6',covariate.type = 'continuous',type ='allD')
d7c = plot_sim(filename = 'allD7',covariate.type = 'continuous',type ='allD')

# binary covariate
c1b = plot_sim(filename = 'allC1',covariate.type = 'binary',type ='allC')
c2b = plot_sim(filename = 'allC2',covariate.type = 'binary',type ='allC')
c3b = plot_sim(filename = 'allC3',covariate.type = 'binary',type ='allC')
c4b = plot_sim(filename = 'allC4',covariate.type = 'binary',type ='allC')
c5b = plot_sim(filename = 'allC5',covariate.type = 'binary',type ='allC')
c6b = plot_sim(filename = 'allC6',covariate.type = 'binary',type ='allC')
c7b = plot_sim(filename = 'allC7',covariate.type = 'binary',type ='allC')
d1b = plot_sim(filename = 'allD1',covariate.type = 'binary',type ='allD')
d2b = plot_sim(filename = 'allD2',covariate.type = 'binary',type ='allD')
d3b = plot_sim(filename = 'allD3',covariate.type = 'binary',type ='allD')
d4b = plot_sim(filename = 'allD4',covariate.type = 'binary',type ='allD')
d5b = plot_sim(filename = 'allD5',covariate.type = 'binary',type ='allD')
d6b = plot_sim(filename = 'allD6',covariate.type = 'binary',type ='allD')
d7b = plot_sim(filename = 'allD7',covariate.type = 'binary',type ='allD')

getwd()
save(c1b,file = 'SimulationEvaluation/c1b.Rdata')
save(c2b,file = 'SimulationEvaluation/c2b.Rdata')
save(c3b,file = 'SimulationEvaluation/c3b.Rdata')
save(c4b,file = 'SimulationEvaluation/c4b.Rdata')
save(c5b,file = 'SimulationEvaluation/c5b.Rdata')
save(c6b,file = 'SimulationEvaluation/c6b.Rdata')
save(c7b,file = 'SimulationEvaluation/c7b.Rdata')
save(d1b,file = 'SimulationEvaluation/d1b.Rdata')
save(d2b,file = 'SimulationEvaluation/d2b.Rdata')
save(d3b,file = 'SimulationEvaluation/d3b.Rdata')
save(d4b,file = 'SimulationEvaluation/d4b.Rdata')
save(d5b,file = 'SimulationEvaluation/d5b.Rdata')
save(d6b,file = 'SimulationEvaluation/d6b.Rdata')
save(d7b,file = 'SimulationEvaluation/d7b.Rdata')


save(c1c,file = 'SimulationEvaluation/c1c.Rdata')
save(c2c,file = 'SimulationEvaluation/c2c.Rdata')
save(c3c,file = 'SimulationEvaluation/c3c.Rdata')
save(c4c,file = 'SimulationEvaluation/c4c.Rdata')
save(c5c,file = 'SimulationEvaluation/c5c.Rdata')
save(c6c,file = 'SimulationEvaluation/c6c.Rdata')
save(c7c,file = 'SimulationEvaluation/c7c.Rdata')
save(d1c,file = 'SimulationEvaluation/d1c.Rdata')
save(d2c,file = 'SimulationEvaluation/d2c.Rdata')
save(d3c,file = 'SimulationEvaluation/d3c.Rdata')
save(d4c,file = 'SimulationEvaluation/d4c.Rdata')
save(d5c,file = 'SimulationEvaluation/d5c.Rdata')
save(d6c,file = 'SimulationEvaluation/d6c.Rdata')
save(d7c,file = 'SimulationEvaluation/d7c.Rdata')


for(i in 1:7){
  for(j in c('c','d')){
    load(paste0('SimulationEvaluation/',j,i,'b','.Rdata'))
  }
}


c6b$recommend2
R = list()
for(file in c(paste0('d',1:7,'b'))){# ,paste0('d',1:7,'b')
  rc = get(file)[['recommend2']]
  r = NULL
  for(i in 1:ncol(rc)){
    r = c(r,paste(rc[,i,drop=F][!is.na(rc[,i,drop=F]),,drop = F] %>% row.names(), collapse = ';'))
  }
  names(r) <- colnames(rc)
  R[[file]] = r 
}

XX = list()
for(i in 1:length(R)){
  x = as.data.frame(R[[i]]) %>% rownames_to_column('diff.otu.mode') %>% separate(diff.otu.mode, c('diff.otu.mode','signal density','impact'),sep = '_') 
  colnames(x)[4] = 'recommendation'
  cat(names(R)[i],'\n')
  print(x)
  x[!duplicated(x[,c('recommendation','diff.otu.mode','signal density','impact')]),]
  XX[[names(R)[i]]] = x
}

DF = list()
for(i in 1:length(XX)){
  df = XX[[i]] %>% tidyr::spread(impact,recommendation) #%>% unite('id',1:2)
  DF[[names(XX)[i]]] = df
}

DFF = DF[[1]]
for(i in 2:7){
  DFF = merge(DFF,DF[[i]], by =c('diff.otu.mode','signal density'))
}

DFF = DFF[c(2,3,1,5,6,4),]
DFF1 = rbind(c('','',rep('Effect size',3),rep('Sample size',3),rep('Depth',3),rep('Taxa number',3),rep('Compostional',3),rep('Compositional+top5',3),rep('Depth confounding',3)),DFF)
# write.csv(DFF1, '~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/recommendation_top1_vaginal.csv', row.names = F)

rbind(rc, r)
mat = matrix(NA, nrow = )



star3.c1b = lapply(seq_len(ncol(c1b$recommend)), function(i) names(c1b$recommend[,i][grep('^\\*\\*\\*$',c1b$recommend[,i])]))
star2.c1b = lapply(seq_len(ncol(c1b$recommend)), function(i) names(c1b$recommend[,i][grep('^\\*\\*$',c1b$recommend[,i])]))
names(star2.c1b) = colnames(c1b$recommend)
names(star3.c1b) = colnames(c1b$recommend)

length(grep('\\*\\*\\*',c1b$recommend[1,]))

# for(i in 1:7){
#   for(j in c('c','d')){
#     for(k in c('c','b')){
#       saveRDS(get(paste0(j,i,k)), file = paste0('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/',j,i,k,'.Rds'))
#     }
#   }
# }

cal <- function(data = c1b, filename = 'summary',name = 'EffectSize_',measure = 'FDR'){
  d1bFDR = data[[filename]] %>% filter(measures %in% measure) %>% 
    dplyr::select(c('methods','measures','median','value')) %>% 
    group_by(methods) %>% 
    summarise(medianOFmedian=median(median), maxOFmedian = max(median)) %>% na.omit() %>% column_to_rownames('methods')
  colnames(d1bFDR) = paste0(name,colnames(d1bFDR))
  return(df =d1bFDR)
}


cal.max = function(data = c1b, name = 'EffectSize_',measure = c('FDR')){
  d1bFDR = data[['res2']] %>% filter(measures %in% measure) %>% 
    dplyr::select(c('methods','value')) %>% 
    group_by(methods) %>%
    summarise(max(value) *0.05) %>% column_to_rownames('methods')
  colnames(d1bFDR) = paste0(name,colnames(d1bFDR))
  return(df =d1bFDR)
}

## median FDR/TPR, time, NAs,
## for dinary data vaginal
if(covariate.type == 'continuous'){
  delete = c('eBay(Wilcoxon)','eBay(t-test)','RAIDA','DESeq2','Wrench+DESeq2','edgeR','Wrench+edgeR','Wrench+metagenomeSeq','metagenomeSeq','Aldex2(t-test)','Aldex2(Wilcoxon)','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon','RioNorm2','LINDA2')
}else{
  #delete = c('GMPR+DESeq2','GMPR+edgeR','eBay(Wilcoxon)','Wrench+metagenomeSeq','Aldex2(Wilcoxon)','DESeq2','edgeR','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon','LINDA2','LinDA','RioNorm2')#'Wrench+DESeq2','Wrench+edgeR''RioNorm2',
  delete = c('GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon','LINDA2','LinDA','RioNorm2')#'Wrench+DESeq2','Wrench+edgeR''RioNorm2',
  supp = c('GMPR+DESeq2','GMPR+edgeR','eBay(Wilcoxon)','Wrench+metagenomeSeq','Aldex2(Wilcoxon)','DESeq2','edgeR')
}


FDR.binary = list()
for (j in c('c','d')) {
  for(i in 1:7){
    df= cal(data = get(paste0(j,i,'b')), filename = 'res.df2') %>% dplyr::select(EffectSize_medianOFmedian)
    names = c(`1` = 'EffectSize_',`2`= 'SampleSize_',`3` ='Depth_', `4` ='TaxaNumber_', `5` = 'Compositional_', `6` = 'CompositionalTop5_', `7` = 'DepthConfounding_')
    names1 = c('c' = 'stool', 'd' = 'vaginal')
    colnames(df) = paste0(names[[as.character(i)]],'medianFDR.',names1[[j]])
    df = rownames_to_column(df,'methods')
    FDR.binary[[paste0(j,i,'c')]] = df
  }
}
FDR.binary <- reduce_list(data = FDR.binary) %>% 
  filter(!(methods %in% delete)) %>%
  column_to_rownames('methods') 
FDR.level <- apply(FDR.binary, 2, function(x) ifelse(x <= 1, 'Good', ifelse(x>2,'Poor','Intermediate')))%>% as.data.frame()
dim(FDR.level)

maxFDR.binary = list()
for (j in c('c','d')) {
  for(i in 1:7){
    cat(paste0(j,i,'b'))
    df= cal(data = get(paste0(j,i,'b')), measure = 'FDR', filename = 'res.df2')%>% dplyr::select(EffectSize_maxOFmedian)
    names = c(`1` = 'EffectSize_',`2`= 'SampleSize_',`3` ='Depth_', `4` ='TaxaNumber_', `5` = 'Compositional_', `6` = 'CompositionalTop5_', `7` = 'DepthConfounding_')
    names1 = c('c' = 'stool', 'd' = 'vaginal')
    colnames(df) = paste0(names[[as.character(i)]],'maxMedianFDR.',names1[[j]])
    df = rownames_to_column(df,'methods')
    maxFDR.binary[[paste0(j,i,'c')]] = df
  }
}

maxFDR.binary <- reduce_list(data = maxFDR.binary) %>% filter(!(methods %in% delete)) %>% column_to_rownames('methods') 
maxFDR.level <- apply(maxFDR.binary, 2, function(x) ifelse(x <= 2, 'Good', ifelse(x>4,'Poor','Intermediate')))%>% as.data.frame()
dim(maxFDR.level)
# maxFDR.binary = list()
# for(i in 1:7){
#   df= cal.max(data = get(paste0('d',i,'b')), measure = 'FDR')
#   names = c(`1` = 'EffectSize_',`2`= 'SampleSize_',`3` ='Depth_', `4` ='TaxaNumber_', `5` = 'Compositional_', `6` = 'CompositionalTop5_', `7` = 'DepthConfounding_')
#   colnames(df) = paste0(names[[as.character(i)]],'maxFDR.vaginal')
#   df = rownames_to_column(df,'methods')
#   maxFDR.binary[[paste0('d',i,'b')]] = df
# }
# maxFDR.binary <- reduce_list(maxFDR.binary) %>% column_to_rownames('methods')
# maxFDR.level.vaginal <- apply(maxFDR.binary, 2, function(x) ifelse(x < 0.2, 'Good', ifelse(x>0.4,'Poor','Intermediate'))) %>% as.data.frame()


## median TPR(0.6-0.8)
detach(package:plyr) # or will produce one 
median_TPR.sum = list()
for(j in c('c','d')){
  for(i in 1:7){
    # df = get(paste0(j,i,'b'))[['summary.TPR']] %>% 
    #   filter(measures =='TPR') %>% 
    #   mutate(TPR_score = median/max_median) %>% 
    #   dplyr::select(c('methods','TPR_score')) %>% group_by(methods) %>% 
    #   summarise(median_TPR = median(TPR_score, na.rm = T)) %>% filter(!(methods %in% delete))
    
    # df = get(paste0(j,i,'b'))[['tprscore']] # score sum
    df = get(paste0(j,i,'b'))[['tpr.rank']] 
    df = apply(df, 1, function(x) median(x)) %>% as.data.frame() %>% rownames_to_column('methods')
    
    colnames(df) = c('methods',paste0(j,i,'b'))
    names = c(`1` = 'EffectSize_',`2`= 'SampleSize_',`3` ='Depth_', `4` ='TaxaNumber_', `5` = 'Compositional_', `6` = 'CompositionalTop5_', `7` = 'DepthConfounding_')
    colnames(df)[2] = paste0(names[[as.character(i)]],'medianTPR.',names1[[j]])
    
    
    median_TPR.sum[[paste0(j,i,'b')]] = df
  }
  
}
median_TPR.sum.df <- reduce_list(median_TPR.sum) %>% column_to_rownames('methods')
# median_TPR.level = apply(median_TPR.sum.df, 2, function(x) ifelse(x > 180, 'Good', ifelse(x<=90,'Poor','Intermediate')))
median_TPR.level = apply(median_TPR.sum.df, 2, function(x) ifelse(x > 10, 'Good', ifelse(x<=5,'Poor','Intermediate')))



## time was calculated by mean value of vaginal binary dataset
times = list()
for(j in c('d')){
  for(i in 1:7){
    time = get(paste0(j,i,'b'))[['summary']] %>% filter(measures %in% c('time') & !(methods %in% delete)) %>% mutate(time = value) %>% dplyr::select(c('methods','time'))# 
    colnames(time)[2] = paste0(j,i,'b')
    times[[paste0(j,i,'b')]] = time
  }
}
time = reduce_list(times) %>% column_to_rownames('methods') %>% rowMeans() %>% as.data.frame() 
colnames(time) = 'time'
time.level <- apply(time, 2, function(x) ifelse(x < 30, 'Good', ifelse(x>60,'Poor','Intermediate'))) %>% as.data.frame()
colnames(time.level) = 'Speed'
time.level 


na.df.binary = list()
for(j in c('c','d')){
  for(i in 1:7){
    na = get(paste0(j,i,'b'))[['na.sum']] %>% filter(!(methods %in% delete)) %>% mutate(na = value) %>% dplyr::select(c('methods','na'))# 
    head(na)
    colnames(na)[2] = paste0(j,i,'b')
    na.df.binary[[paste0(j,i,'b')]] = na
  }
}

na = reduce_list(na.df.binary) %>% column_to_rownames('methods') %>% rowMeans()%>% as.data.frame() 
colnames(na) = 'Failure rate'
na.binary.level = apply(na, 2, function(x) ifelse(x < 0.05, 'Good', ifelse(x>0.15,'Poor','Intermediate')))
dim(na.binary.level)
# na.df.binary = cbind(d1b$na.sum$value,d2b$na.sum$value,d3b$na.sum$value,d4b$na.sum$value,d5b$na.sum$value,d6b$na.sum$value,d7b$na.sum$value,
#                      c1b$na.sum$value,c2b$na.sum$value,c3b$na.sum$value,c4b$na.sum$value,c5b$na.sum$value,c6b$na.sum$value,c7b$na.sum$value)
# na.df.continuous = rbind(d0c$na.sum,d1c$na.sum,d2c$na.sum,d3c$na.sum,d4c$na.sum,d5c$na.sum,d6c$na.sum,d7c$na.sum,
#                          c0c$na.sum,c1c$na.sum,c2c$na.sum,c3c$na.sum,c4c$na.sum,c5c$na.sum,c6c$na.sum,c7c$na.sum)
# na.binary = aggregate(value ~ methods, na.df.binary, function(x) (sum(x)/8)) %>% dplyr::rename(`Failure rate`=value)%>% column_to_rownames('methods')
# na.continuous = aggregate(value ~ methods, na.df.continuous, function(x) (sum(x)/8))%>% dplyr::rename(`Failure rate`=value)%>% column_to_rownames('methods')
# na.binary.level = apply(na.binary, 2, function(x) ifelse(x < 0.1, 'Good', ifelse(x>0.25,'Poor','Intermediate')))
# na.continuous.level = apply(na.continuous, 2, function(x) ifelse(x < 0.1, 'Good', ifelse(x>0.25,'Poor','Intermediate')))



## bias DAs 
load('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/StoolbinaryBiasDAs.Rdata')
Stool.bias = aggregate(ct ~ methods+covariate.types + depth.conf.factors, data = m, function(x) median(x)) %>% filter(!(methods %in% delete))
Stool.bias.binary.none = Stool.bias[Stool.bias$covariate.types=='binary' & Stool.bias$depth.conf.factors=='None',] 
Stool.bias.binary.none$score = 'Intermediate'
Stool.bias.binary.none[Stool.bias.binary.none$ct<=0.05,'score'] = 'Good'
Stool.bias.binary.none[Stool.bias.binary.none$ct>0.1,'score'] = 'Poor'
Stool.bias.binary.none =Stool.bias.binary.none[,c('methods','score')]
colnames(Stool.bias.binary.none)[2] = 'BiasDAs_NoneConfounding.stool'
dim(Stool.bias.binary.none)

Stool.bias.binary.conf = Stool.bias[Stool.bias$covariate.types=='binary' & Stool.bias$depth.conf.factors!='None',] %>% filter(!(methods %in% delete))
Stool.bias.binary.conf$score = 'Intermediate'
Stool.bias.binary.conf[Stool.bias.binary.conf$ct<=0.05,'score'] = 'Good'
Stool.bias.binary.conf[Stool.bias.binary.conf$ct>0.1,'score'] = 'Poor'
Stool.bias.binary.conf =Stool.bias.binary.conf[,c('methods','score')]
colnames(Stool.bias.binary.conf)[2] = 'BiasDAs_Confounding.stool'



load('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/VaginalbinaryBiasDAs.Rdata')
Vaginal.bias = aggregate(ct ~ methods+covariate.types + depth.conf.factors, data = m, function(x) median(x)) %>% filter(!(methods %in% delete))
Vaginal.bias.binary.none = Vaginal.bias[Vaginal.bias$covariate.types=='binary' & Vaginal.bias$depth.conf.factors=='None',] 
Vaginal.bias.binary.none$score = 'Intermediate'
Vaginal.bias.binary.none[Vaginal.bias.binary.none$ct<=0.05,'score'] = 'Good'
Vaginal.bias.binary.none[Vaginal.bias.binary.none$ct>0.1,'score'] = 'Poor'
Vaginal.bias.binary.none =Vaginal.bias.binary.none[,c('methods','score')]
colnames(Vaginal.bias.binary.none)[2] = 'BiasDAs_NoneConfounding.vaginal'


Vaginal.bias.binary.conf = Vaginal.bias[Vaginal.bias$covariate.types=='binary' & Vaginal.bias$depth.conf.factors!='None',]  %>% filter(!(methods %in% delete))
Vaginal.bias.binary.conf$score = 'Intermediate'
Vaginal.bias.binary.conf[Vaginal.bias.binary.conf$ct<=0.05,'score'] = 'Good'
Vaginal.bias.binary.conf[Vaginal.bias.binary.conf$ct>0.1,'score'] = 'Poor'
Vaginal.bias.binary.conf =Vaginal.bias.binary.conf[,c('methods','score')]
colnames(Vaginal.bias.binary.conf)[2] = 'BiasDAs_Confounding.vaginal'


load('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/stability/correlationR.Rdata')
head(corr.sum)
colnames(corr.sum) = c('methods','Stability')
corr.sum$methods = as.character(corr.sum$methods)
corr.sum$methods = gsub('Wilcox','TSS+Wilcoxon',corr.sum$methods)
corr.sum$methods = gsub('^Rarefy$','Rarefy+Wilcoxon',corr.sum$methods)
corr.sum$methods = gsub('ANCOMBC','ANCOM-BC',corr.sum$methods)
corr.sum$methods = gsub('glmquassi','GLM(quasipoisson)',corr.sum$methods)
corr.sum$methods = gsub('DESeq2.gmpr','GMPR+DESeq2',corr.sum$methods)
corr.sum$methods = gsub('DESeq2.Wrench','Wrench+DESeq2',corr.sum$methods)
corr.sum$methods = gsub('edgeR.gmpr','GMPR+edgeR',corr.sum$methods)
corr.sum$methods = gsub('edgeR.Wrench','Wrench+edgeR',corr.sum$methods)
corr.sum$methods = gsub('MSeq2.Wrench','Wrench+metagenomeSeq',corr.sum$methods)
corr.sum$methods = gsub('MSeq2','metagenomeSeq',corr.sum$methods)
corr.sum$methods = gsub('eBayW','eBay(Wilcoxon)',corr.sum$methods)
corr.sum$methods = gsub('eBayt','eBay(t-test)',corr.sum$methods)
corr.sum$methods = gsub('BBinomial','Beta-binomial',corr.sum$methods)
corr.sum$methods = gsub('Aldex2we','Aldex2(t-test)',corr.sum$methods)
corr.sum$methods = gsub('^Aldex2$','Aldex2(Wilcoxon)',corr.sum$methods)
corr.sum$methods = gsub('^ttest$',"TSS+t-test",corr.sum$methods)
corr.sum$methods = gsub('^Rarefyttest$','Rarefy+t-test',corr.sum$methods)
corr.sum$methods = gsub('CLRBC',"LinDA",corr.sum$methods)

corr.sum  = corr.sum %>% column_to_rownames('methods')
# corr.sum = corr.sum[rownames(heat.binary),,drop = F] 
str(corr.sum)
corr.sum  = apply(corr.sum, 2, function(x) ifelse(x >0.9, 'Good', ifelse(x<0.7,'Poor','Intermediate'))) %>% as.data.frame() 
corr.sum = corr.sum[rownames(FDR.level),,drop = F]
dim(corr.sum)
heat.binary = full_join(FDR.level %>% rownames_to_column('methods'), as.data.frame(median_TPR.level) %>% rownames_to_column('methods')) %>%
  full_join(maxFDR.level %>% rownames_to_column('methods')) %>%
  full_join(time.level%>% rownames_to_column('methods')) %>% full_join(as.data.frame(na.binary.level) %>% rownames_to_column('methods')) %>%
  full_join(Vaginal.bias.binary.conf) %>% full_join(Vaginal.bias.binary.none) %>% 
  full_join(Stool.bias.binary.conf) %>% full_join(Stool.bias.binary.none) %>% full_join(corr.sum %>% rownames_to_column('methods')) %>% filter(!(methods %in% delete))%>%column_to_rownames('methods')
## for supplementary
heat.binary = full_join(FDR.level %>% rownames_to_column('methods'), as.data.frame(median_TPR.level) %>% rownames_to_column('methods')) %>%
  full_join(maxFDR.level %>% rownames_to_column('methods')) %>%
  full_join(time.level%>% rownames_to_column('methods')) %>% full_join(as.data.frame(na.binary.level) %>% rownames_to_column('methods')) %>%
  full_join(Vaginal.bias.binary.conf) %>% full_join(Vaginal.bias.binary.none) %>% 
  full_join(Stool.bias.binary.conf) %>% full_join(Stool.bias.binary.none) %>% full_join(corr.sum %>% rownames_to_column('methods')) %>% filter((methods %in% supp))%>%column_to_rownames('methods')


heat.binary$`Complex design` = NA
heat.binary[rownames(heat.binary) %in% c('GMPR+DESeq2','GMPR+edgeR','LDM','ANCOM-BC',
                                         'GLM(quasipoisson)','Beta-binomial','DESeq2','edgeR'),'Complex design'] = 'Good'
heat.binary[rownames(heat.binary) %in% c('DACOMP','Wrench+DESeq2','Wrench+edgeR'),'Complex design'] = 'Intermediate'
heat.binary[rownames(heat.binary) %in% c('metagenomeSeq','Aldex2(Wilcoxon)','Aldex2(t-test)','mbzinb','eBay(t-test)', 'eBay(Wilcoxon)','RAIDA', 'Rarefy+Wilcoxon','TSS+Wilcoxon','TSS+t-test','Wrench+metagenomeSeq','Rarefy+t-test'),'Complex design'] = 'Poor'


heat = apply(heat.binary,2, function(x) as.factor(x))
rownames(heat) = rownames(heat.binary)


dd = heat
dd[dd=="Poor"] = 1
dd[dd=="Intermediate"] = 2
dd[dd=="Good"] = 3
str(dd)
dd = apply(dd, 2, function(x)(as.numeric(as.character(x)))) %>% as.data.frame()
rownames(dd) = rownames(heat)
rownames(dd);colnames(dd)
d = dist(dd, method = "euclidean")
d1 = hclust(d, method="ward.D2")
ord.cl = rownames(dd)[d1$order]
# pheatmap::pheatmap(t(dd),legend = TRUE, cluster_rows = F, 
#                    cluster_cols = T,clustering_method = 'ward.D2', color =c('red','yellow','green'),
#                    annotation_legend = F,show_colnames = T,show_rownames =T,fontsize_row =6,fontsize_col =6)

d = melt(heat)
colnames(d)[1:2] = c('methods','level')
table(d$methods, d$value)
d$level = gsub('vaginal','low',d$level)
d$level = gsub('stool','high',d$level)

head(d)
d[,1] <- as.factor(d[,1])
d[,2] <- as.factor(d[,2])
medianfdr =c(as.character(unique(d$level)[grep('medianFDR.low',unique(d$level))]),
             as.character(unique(d$level)[grep('MedianFDR.low',unique(d$level))]),
             as.character(unique(d$level)[grep('medianTPR.low',unique(d$level))]),
             as.character(unique(d$level)[grep('_Confounding.low',unique(d$level))]),      
             as.character(unique(d$level)[grep('NoneConfounding.low',unique(d$level))]),
             as.character(unique(d$level)[grep('medianFDR.high',unique(d$level))]),
             as.character(unique(d$level)[grep('MedianFDR.high',unique(d$level))]),
             as.character(unique(d$level)[grep('medianTPR.high',unique(d$level))]),
             as.character(unique(d$level)[grep('_Confounding.high',unique(d$level))]),      
             as.character(unique(d$level)[grep('NoneConfounding.high',unique(d$level))]),
             as.character(unique(d$level)[grep('Failure rate',unique(d$level))]),
                           as.character(unique(d$level)[grep('Speed',unique(d$level))]),'Complex design','Stability')
d[,2] <- factor(d[,2], levels = medianfdr)
d[,1] <- factor(d[,1], levels = ord.cl)
d$color = 'black'
d[grep('medianTPR.high$',d[,2]),'color'] = brewer.pal(9,'BrBG')[2]
d[grep('medianFDR.high$',d[,2]),'color'] =brewer.pal(9,'BrBG')[4]
d[grep('medianTPR.low$',d[,2]),'color'] =brewer.pal(9,'PuBu')[8]
d[grep('medianFDR.low$',d[,2]),'color'] =brewer.pal(9,'PuBu')[6]
d$color <- as.factor(d$color)
cols1 = c("#BF812D" =brewer.pal(9,'BrBG')[2], "#F6E8C3" = brewer.pal(9,'BrBG')[4], "#045A8D" =brewer.pal(9,'PuBu')[8],"#FFFF33"=brewer.pal(9,'Set1')[6],'black'='black')
cols = c('Good' =brewer.pal(9,'Set1')[3], 'Intermediate' = brewer.pal(8,'Set2')[6], 'Poor' =brewer.pal(9,'Set1')[1])
cols1 =c(rep(brewer.pal(9,'BrBG')[1],7),
         rep(brewer.pal(9,'BrBG')[2],7),
         rep(brewer.pal(9,'BrBG')[3],7),
         rep(brewer.pal(8,'Set1')[4],2),
         rep(brewer.pal(9,'PuBu')[8],7), 
         rep(brewer.pal(9,'PuBu')[6],7),
         rep(brewer.pal(9,'PuBu')[4],7),
         rep(brewer.pal(8,'Set1')[4],2),
         rep('black',4))
ggplot(d, aes(x = level, y = methods)) +
  geom_tile(colour = "white", aes(fill = value)) +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(color="black", size =12, angle = 90, hjust = 1, vjust = 0.5),
             axis.text.y = element_text(size = 12, colour = cols1),
             axis.title = element_text(color="black", size = 12)) + labs(x = '',y = '', fill = '', title = 'Binary') +coord_flip()
ggsave('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/SimulationEvaluation/heatmap_binary.pdf', width = 9, height = 12, dpi = 100)
ggsave('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/SimulationEvaluation/heatmap_binary_supplementary.pdf', width = 6, height = 12, dpi = 100)
getwd()










## ========================Continuous =================================
##============ for continuous 
detach(package:plyr) # or will produce one 
delete = c("Wrench+metagenomeSeq","Wrench+DESeq2","edgeR","Wrench+edgeR","DESeq2","Aldex2(t-test)","eBay(t-test)","CLRBC")
delete = c('eBay(Wilcoxon)','eBay(t-test)','RAIDA','metagenomeSeq',"Rarefy+Wilcoxon","TSS+Wilcoxon",'Wrench+DESeq2','Wrench+edgeR','Wrench+metagenomeSeq','Aldex2(t-test)','Aldex2(Wilcoxon)','CLRBC')

# FDR-continuous
FDR.continuous = list()
for (j in c('c','d')) {
  for(i in 1:7){
    df= cal(data = get(paste0(j,i,'c')), measure = 'FDR')%>% dplyr::select(EffectSize_medianOFmedian)
    names = c(`1` = 'EffectSize_',`2`= 'SampleSize_',`3` ='Depth_', `4` ='TaxaNumber_', `5` = 'Compositional_', `6` = 'CompositionalTop5_', `7` = 'DepthConfounding_')
    names1 = c('c' = 'stool', 'd' = 'vaginal')
    colnames(df) = paste0(names[[as.character(i)]],'medianFDR.',names1[[j]])
    df = df %>% rownames_to_column('methods')
    FDR.continuous[[paste0(j,i,'c')]] = df
  }
}
FDR.continuous <- reduce_list(data = FDR.continuous) %>% filter(!(methods %in% delete)) %>% column_to_rownames('methods') 
FDR.level <- apply(FDR.continuous, 2, function(x) ifelse(x < 0.05, 'Good', ifelse(x>0.2,'Poor','Intermediate')))%>% as.data.frame()

## median TPR(0.6-0.8)
detach(package:plyr) # or will produce one 
median_TPR.sum = list()
for(j in c('c','d')){
  for(i in 1:7){
    df = get(paste0(j,i,'c'))[['summary.TPR']] %>% 
      filter(measures =='TPR') %>% 
      mutate(TPR_score = median/max_median) %>% 
      dplyr::select(c('methods','TPR_score')) %>% group_by(methods) %>% 
      summarise(median_TPR = median(TPR_score, na.rm = T)) %>% filter(!(methods %in% delete))
    names = c(`1` = 'EffectSize_',`2`= 'SampleSize_',`3` ='Depth_', `4` ='TaxaNumber_', `5` = 'Compositional_', `6` = 'CompositionalTop5_', `7` = 'DepthConfounding_')
    names1 = c('c' = 'stool', 'd' = 'vaginal')
    colnames(df)[2] = paste0(names[[as.character(i)]],'medianTPR.',names1[[j]])
    median_TPR.sum[[paste0(j,i,'c')]] = df
  }
}

median_TPR.sum.df <- reduce_list(median_TPR.sum) %>% column_to_rownames('methods')
median_TPR.level = apply(median_TPR.sum.df, 2, function(x) ifelse(x > 0.7, 'Good', ifelse(x<=0.4,'Poor','Intermediate')))


## time was calculated by mean value of vaginal binary dataset
time = d1c[['res.df2']] %>% filter(d1c$res.df2$measures %in% c('time')) %>% mutate(time = value) %>% dplyr::select(c('methods','time')) %>% column_to_rownames('methods')
time$time = time$time/max(time$time)
time.level <- apply(time, 2, function(x) ifelse(x < 0.1, 'Good', ifelse(x>0.7,'Poor','Intermediate'))) %>% as.data.frame()
colnames(time.level) = 'Speed'


na.df.continuous = list()
for(j in c('c','d')){
  for(i in 1:7){
    na = get(paste0(j,i,'b'))[['na.sum']] %>% filter(!(methods %in% delete)) %>% mutate(na = value) %>% dplyr::select(c('methods','na'))# 
    head(na)
    colnames(na)[2] = paste0(j,i,'b')
    na.df.continuous[[paste0(j,i,'b')]] = na
  }
}

na = reduce_list(na.df.continuous) %>% column_to_rownames('methods') %>% rowMeans()%>% as.data.frame() 
colnames(na) = 'Failure rate'
na.continuous.level = apply(na, 2, function(x) ifelse(x < 0.05, 'Good', ifelse(x>0.15,'Poor','Intermediate')))
dim(na.continuous.level)




## bias DAs 
load('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/VaginalcontinuousBiasDAs.Rdata')
Stool.bias = aggregate(ct ~ methods+covariate.types + depth.conf.factors, data = m, function(x) median(x)) 
Stool.bias.continuous.none = Stool.bias[Stool.bias$covariate.types=='continuous' & Stool.bias$depth.conf.factors=='None',] 
Stool.bias.continuous.none$score = 'Intermediate'
Stool.bias.continuous.none[Stool.bias.continuous.none$ct<=0.05,'score'] = 'Good'
Stool.bias.continuous.none[Stool.bias.continuous.none$ct>0.1,'score'] = 'Poor'
Stool.bias.continuous.none =Stool.bias.continuous.none[,c('methods','score')]
colnames(Stool.bias.continuous.none)[2] = 'BiasDAs_NoneConfounding.stool'


Stool.bias.continuous.conf = Stool.bias[Stool.bias$covariate.types=='continuous' & Stool.bias$depth.conf.factors!='None',] 
Stool.bias.continuous.conf$score = 'Intermediate'
Stool.bias.continuous.conf[Stool.bias.continuous.conf$ct<=0.05,'score'] = 'Good'
Stool.bias.continuous.conf[Stool.bias.continuous.conf$ct>0.1,'score'] = 'Poor'
Stool.bias.continuous.conf =Stool.bias.continuous.conf[,c('methods','score')]
colnames(Stool.bias.continuous.conf)[2] = 'BiasDAs_Confounding.stool'



load('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/VaginalcontinuousBiasDAs.Rdata')

Vaginal.bias = aggregate(ct ~ methods+covariate.types + depth.conf.factors, data = m, function(x) median(x)) 
Vaginal.bias.continuous.none = Vaginal.bias[Vaginal.bias$covariate.types=='continuous' & Vaginal.bias$depth.conf.factors=='None',] 
Vaginal.bias.continuous.none$score = 'Intermediate'
Vaginal.bias.continuous.none[Vaginal.bias.continuous.none$ct<=0.05,'score'] = 'Good'
Vaginal.bias.continuous.none[Vaginal.bias.continuous.none$ct>0.1,'score'] = 'Poor'
Vaginal.bias.continuous.none =Vaginal.bias.continuous.none[,c('methods','score')]
colnames(Vaginal.bias.continuous.none)[2] = 'BiasDAs_NoneConfounding.vaginal'


Vaginal.bias.continuous.conf = Vaginal.bias[Vaginal.bias$covariate.types=='continuous' & Vaginal.bias$depth.conf.factors!='None',] 
Vaginal.bias.continuous.conf$score = 'Intermediate'
Vaginal.bias.continuous.conf[Vaginal.bias.continuous.conf$ct<=0.05,'score'] = 'Good'
Vaginal.bias.continuous.conf[Vaginal.bias.continuous.conf$ct>0.1,'score'] = 'Poor'
Vaginal.bias.continuous.conf =Vaginal.bias.continuous.conf[,c('methods','score')]
colnames(Vaginal.bias.continuous.conf)[2] = 'BiasDAs_Confounding.vaginal'




heat.continuous = full_join(FDR.level.vaginal %>% rownames_to_column('methods'), FDR.level.stool %>% rownames_to_column('methods')) %>% 
  full_join(as.data.frame(median_TPR.level) %>% rownames_to_column('methods')) %>%
  full_join(time.level%>% rownames_to_column('methods')) %>%
  full_join(as.data.frame(na.continuous.level) %>% rownames_to_column('methods')) %>%
  full_join(Vaginal.bias.continuous.conf) %>% full_join(Vaginal.bias.continuous.none) %>% 
  full_join(Stool.bias.continuous.conf) %>% full_join(Stool.bias.continuous.none)  %>% filter(!(methods %in% delete))%>%column_to_rownames('methods')
heat.continuous$`Failure rate`[is.na(heat.continuous$`Failure rate`)] = 'Good'
heat = apply(heat.continuous,2, function(x) as.factor(x))
rownames(heat) = rownames(heat.continuous)

d = melt(heat)
d$X2 <- as.factor(d$X2)
d$X1 <- as.factor(d$X1)
medianfdr =c(as.character(unique(d$X2)[grep('medianFDR.vaginal',unique(d$X2))]),
             as.character(unique(d$X2)[grep('medianTPR.vaginal',unique(d$X2))]),
             as.character(unique(d$X2)[grep('medianFDR.stool',unique(d$X2))]),
             as.character(unique(d$X2)[grep('medianTPR.stool',unique(d$X2))]),
             as.character(unique(d$X2)[grep('_Confounding',unique(d$X2))]),      
             as.character(unique(d$X2)[grep('NoneConfounding',unique(d$X2))]),
             as.character(unique(d$X2)[grep('Failure rate',unique(d$X2))]),
             as.character(unique(d$X2)[grep('Speed',unique(d$X2))]))
d$X2 <- factor(d$X2, levels = medianfdr)
d$X1 <- factor(d$X1, levels = c('LDM','ANCOM-BC','Omnibus','GLM(quasipoisson)',
                                'Rarefy+Spearman','Beta-binomial','DACOMP','TSS+Spearman',    
                                'GMPR+DESeq2','DESeq2','GMPR+edgeR','edgeR'))
d$color = 'black'
d[grep('medianTPR.stool$',d$X2),'color'] = brewer.pal(9,'BrBG')[2]
d[grep('medianFDR.stool$',d$X2),'color'] =brewer.pal(9,'BrBG')[4]
d[grep('medianTPR.vaginal$',d$X2),'color'] =brewer.pal(9,'PuBu')[8]
d[grep('medianFDR.vaginal$',d$X2),'color'] =brewer.pal(9,'PuBu')[6]
d$color <- as.factor(d$color)
cols1 = c("#BF812D" =brewer.pal(9,'BrBG')[2], "#F6E8C3" = brewer.pal(9,'BrBG')[4], "#045A8D" =brewer.pal(9,'PuBu')[8],"#FFFF33"=brewer.pal(9,'Set1')[6],'black'='black')
cols = c('Good' =brewer.pal(9,'Set1')[3], 'Intermediate' = brewer.pal(8,'Set2')[6], 'Poor' =brewer.pal(9,'Set1')[1])
cols1 =c(rep(brewer.pal(9,'BrBG')[1],7),
         rep(brewer.pal(9,'BrBG')[2],7),
         rep(brewer.pal(9,'PuBu')[8],7), 
         rep(brewer.pal(9,'PuBu')[6],7),
         rep(brewer.pal(8,'Set1')[4],4),
         rep('black',3))

ggplot(d, aes(x = X2, y = X1)) +
  geom_tile(colour = "white", aes(fill = value)) +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(color="black", size =12, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12, colour = cols1),
        axis.title = element_text(color="black", size = 12)) + labs(x = '',y = '', fill = '') +coord_flip() +
  labs(title = 'Continuous')
ggsave('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/heatmap.continuous.pdf', width = 9, height = 8, dpi = 100)

