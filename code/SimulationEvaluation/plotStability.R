pkg =  c('ggpubr','microbiome',"eBay","modeest","ANCOMBC","aod","phyloseq","reshape","corncob","MASS","readr","DESeq2", "ALDEx2", "metagenomeSeq", "edgeR", "GUniFrac", "grDevices", "dirmult", "exactRankTests","nlme", "dplyr", "magrittr", "tidyr", "protoclust", "ggplot2", "compositions","rmutil","tibble","reticulate","dacomp","LDM","Wrench")
sapply(pkg, require, character = TRUE)
setwd('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/stability/')

methods <- c('mbzinb','ttest','ttest.gmpr','ttest.Wrench','Wilcox.Wrench','Wilcox.gmpr','CLRBC','eBayt','eBayW','Rarefy', 'Aldex2', 'RAIDA', 'DACOMP', 'LDM', 'glmquassi','BBinomial','ANCOMBC','DESeq2', 'DESeq2.Wrench', 'DESeq2.gmpr','Wilcox','edgeR', 'edgeR.Wrench', 'edgeR.gmpr','MSeq2','MSeq2.Wrench','Rarefyttest','Aldex2we')
P = EST = DF = list()
for(method in methods){
  files = NULL
  for(i in 1:100){
    file = paste0("allC5f_summarymatrix",i,"-",method,".Rdata")
    files = c(files, file)
    }
  p = est = NULL
  df2 = NULL
  
  for(file in files){
    load(paste0('0/',file))
    res_seqs0 = res_seqs
    names0 = names(res_seqs0)
    names0 = names0[grep('binary',names0)] %>% na.omit()
    
    load(paste0('40/',file))
    res_seqs40 = res_seqs
    names40 = names(res_seqs40)
    names40 = names40[grep('binary',names40)] %>% na.omit()
    names = intersect(names0, names40)
    for(name in names){
      df0 = res_seqs0[[name]]
      df40 = res_seqs40[[name]]
      cat(nrow(df0),':', nrow(df40),'\n')
      df = inner_join(df0 %>% dplyr::select(otu.id, fdr) %>% dplyr::rename(nofilter = fdr),
                      df40 %>% dplyr::select(otu.id, fdr) %>% dplyr::rename(filter = fdr))
      df1 = df %>% mutate(data = name)
      df2 = rbind(df2, df1)
      cor = cor.test(df$nofilter,df$filter, method = 'spearman')
      p = c(p,cor$p.value)
      est = c(est,cor$estimate)
    }
  }
  DF[[method]] = df2
  P[[method]] = p
  EST[[method]] = est
}

# save(P,file= 'correlationP.Rdata')
# save(EST,file= 'correlationR.Rdata')
# save(DF,file= 'correlationDataframe.Rdata')

corr= NULL
for(i in 1:length(EST)){
  cor = as.data.frame(EST[[i]])
  colnames(cor) = 'R'
  cor$method = names(EST)[i]
  corr = rbind(corr, cor)
}
head(corr)

corr.sum = aggregate(R ~ ., data =corr, function(x) mean(x))

# save(corr.sum,file= 'correlationR.Rdata')

# make scatterplot for each method
# load('correlationDataframe.Rdata')
head(DF$mbzinb)
methods <- c('mbzinb','ttest','eBayt','eBayW','Rarefy', 'Aldex2', 'RAIDA', 'DACOMP', 'LDM', 'glmquassi','BBinomial','ANCOMBC','DESeq2', 'DESeq2.Wrench', 'DESeq2.gmpr','Wilcox','edgeR', 'edgeR.Wrench', 'edgeR.gmpr','MSeq2','MSeq2.Wrench','Rarefyttest','Aldex2we')
Methods <- c('mbzinb','TSS+t-test','eBay(t-test)','eBay(Wilcoxon)','Rarefy+Wilcoxon', 'Aldex2(Wilcoxon)', 'RAIDA', 'DACOMP', 'LDM', 'GLM(quasipoisson)','Beta-binomial','ANCOM-BC',
             'DESeq2', 'Wrench+DESeq2', 'GMPR+DESeq2','TSS+Wilcoxon','edgeR', 'Wrench+edgeR', 'GMPR+edgeR','metagenomeSeq','Wrench+metagenomeSeq','Rarefy+t-test','Aldex2(t-test)')
for(i in 1:length(methods)){
  method = methods[i]
  data = DF[[method]]
  data$data = gsub('loglinearSub_L1nOTU_L5binarynone|nSam_L2|^D3|none|unbalanced','',data$data)
  # each has 18 facets
  unq = unique(data$data)
  data$diff.otu.modes = 'rare'
  data[grep('abundant',data$data),'diff.otu.modes']= 'abundant'
  
  data$signaldensity = 'low'
  data[grep('medium',data$data),'signaldensity']= 'medium'
  data[grep('high',data$data),'signaldensity']= 'high'
  data$data = gsub(method,'',data$data)
  data$data = gsub('low|medium|high|abundant|rare','',data$data)
  data = data %>% dplyr::select(otu.id, nofilter, filter, diff.otu.modes, signaldensity, data)
  #graphics.off()
  p1 = ggplot(data%>% filter(diff.otu.modes =='abundant'), aes(x = nofilter, y = filter)) +
    geom_point(size = 0.1) +theme_bw() +
    facet_wrap(data ~ signaldensity, ncol = 5)+
    labs(x = 'No filter', y = 'Filter by 40% prevelance') +
    theme(axis.text.x = element_text(color="black", size =26),
          axis.text.y = element_text(color="black", size = 26),
          axis.title = element_text(color="black", size = 30),
          strip.text = element_text(size = 30),
          strip.background = element_rect(fill="white",color = "black", size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black",size= 1),
          legend.position = 'none',
          legend.title=element_text(size=30),
          legend.text = element_text(size=30),
          plot.title = element_text(size=22))+ ggtitle(paste0('abundant','_',Methods[i]))
  p2 = ggplot(data%>% filter(diff.otu.modes =='rare'), aes(x = nofilter, y = filter)) +
    geom_point(size = 0.1) +theme_bw() +
    facet_wrap(data ~ signaldensity, ncol = 5)+
    labs(x = 'No filter', y = 'Filter by 40% prevelance') +
    theme(axis.text.x = element_text(color="black", size =26),
          axis.text.y = element_text(color="black", size = 26),
          axis.title = element_text(color="black", size = 30),
          strip.text = element_text(size = 30),
          strip.background = element_rect(fill="white",color = "black", size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black",size= 1),
          legend.position = 'none',
          legend.title=element_text(size=30),
          legend.text = element_text(size=30),
          plot.title = element_text(size=22)) +ggtitle(paste0('rare','_',Methods[i]))
  pp = ggarrange(p1, p2, nrow = 2)
  ggsave(file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/SimulationEvaluation/stability/',method,'.pdf'), width = 35, height = 40, dpi = 100)
  
}



