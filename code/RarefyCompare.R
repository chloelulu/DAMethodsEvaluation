pkg = c('dplyr','tidyr','reshape','RColorBrewer','ggplot2','ggpubr','stringi','scales','tibble')
suppressPackageStartupMessages(sapply(pkg, require, character = T))

source('~/Documents/Mayo_project/Code/DailyCode.R')
load('/Users/m216453/Documents/Mayo_Research/SemiSimulation/0925/sim/new/sim/allC70_res.Rdata')
rarefy = melt(res)
table(rarefy$X14)
load('/Users/m216453/Documents/Mayo_Research/SemiSimulation/0925/sim/new/sim/allC7_res.Rdata')
rarefyF = melt(res)
colnames(rarefyF) <- colnames(rarefy) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
rarefy$methods = sub('^Omnibus$','mbzinb',rarefy$methods)
table(rarefy$methods)


rarefy = rarefy %>% dplyr::filter( measures %in% c('FDR','TPR') & depth.conf.factors %in% c('DL1','DL2','DL3')) %>%
  dplyr::select(c('diff.otu.modes','covariate.types','depth.conf.factors', 'diff.otu.pcts','measures', 'methods', 'value','iters')) %>% droplevels()
rarefyF = rarefyF %>% filter(measures %in% c('FDR','TPR') & depth.conf.factors %in% c('DL1','DL2','DL3'))%>%
  dplyr::select(c('diff.otu.modes','covariate.types','depth.conf.factors', 'diff.otu.pcts','measures', 'methods', 'value','iters'))%>% droplevels()
levels(rarefy$depth.conf.factors) <- levels(rarefyF$depth.conf.factors) <- c('Depth confounding +','Depth confounding ++','Depth confounding +++')

rarefy_FDR  = rarefy[rarefy$measures =='FDR',]
rarefy_FDR$value = rarefy_FDR$value
rarefy_noFDR  = rarefy[rarefy$measures !='FDR',]
rarefy = rbind(rarefy_FDR, rarefy_noFDR)
rarefy <- rarefy[!is.na(rarefy$value),] %>% mutate(grp = 'Yes')


rarefyF_FDR  = rarefyF[rarefyF$measures =='FDR',]
rarefyF_FDR$value = rarefyF_FDR$value
rarefyF_noFDR  = rarefyF[rarefyF$measures !='FDR',]
rarefyF = rbind(rarefyF_FDR, rarefyF_noFDR)
rarefyF <- rarefyF[!is.na(rarefyF$value),] %>% mutate(grp = 'No')


res2 = rbind(rarefy, rarefyF) %>% filter(covariate.types =='binary') %>% droplevels()
table(res2$methods, res2$grp)
res2[res2$grp =='No' & res2$methods =='Rarefyttest','grp'] = 'Yes'
res2[res2$grp =='Yes' & res2$methods =='Rarefyttest','methods'] = 'ttest'
res2 = res2 %>% filter(!(methods %in% c('Rarefy','CLRBC')))



levels(res2$diff.otu.modes) = c("Abundant", "Rare")
levels(res2$diff.otu.pcts) = c("Low density", "Medium density","High density")

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
res2$methods = as.character(res2$methods)
res2$methods[res2$methods =='Wilcox' & res2$covariate.types =='binary'] = 'TSS+Wilcoxon'
res2$methods[res2$methods =='Wilcox' & res2$covariate.types =='continuous'] = 'TSS+Spearman'
res2$methods[res2$methods =='Rarefy' & res2$covariate.types =='binary'] = 'Rarefy+Wilcoxon'
res2$methods[res2$methods =='Rarefy' & res2$covariate.types =='continuous'] = 'Rarefy+Spearman'
res2$methods = gsub('ANCOMBC','ANCOM-BC',res2$methods)
res2$methods = gsub('glmquassi','GLM(quasipoisson)',res2$methods)
res2$methods = gsub('DESeq2.gmpr','GMPR+DESeq2',res2$methods)
res2$methods = gsub('DESeq2.Wrench','Wrench+DESeq2',res2$methods)
res2$methods = gsub('edgeR.gmpr','GMPR+edgeR',res2$methods)
res2$methods = gsub('edgeR.Wrench','Wrench+edgeR',res2$methods)
res2$methods = gsub('MSeq2.Wrench','Wrench+metagenomeSeq',res2$methods)
res2$methods = gsub('MSeq2','metagenomeSeq',res2$methods)
res2$methods = gsub('eBayW','eBay(Wilcoxon)',res2$methods)
res2$methods = gsub('eBayt','eBay(t-test)',res2$methods)
res2$methods = gsub('BBinomial','Beta-binomial',res2$methods)
res2$methods = gsub('Aldex2we','Aldex2(t-test)',res2$methods)
res2$methods = gsub('^Aldex2$','Aldex2(Wilcoxon)',res2$methods)
res2$methods = gsub('^Rarefyttest$','Rarefy+t-test',res2$methods)
res2$methods = gsub('^ttest$','TSS+t-test',res2$methods)


res.df2 <- data_summary(data= res2, formula = paste0('value ~ covariate.types +diff.otu.modes+depth.conf.factors+ diff.otu.pcts+ methods + measures'))

measure = 'TPR'
measure1 = 'FDR'
measure2 = 'na'
ylab = 'observed FDR/expected FDR(0.05)'
grid.formula = 'diff.otu.pcts ~ diff.otu.modes'
delete = c('DESeq2','Wrench+DESeq2','edgeR','Wrench+edgeR','Wrench+metagenomeSeq','eBay(t-test)','Aldex2(t-test)','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon','RioNorm2','LINDA2','LinDA')

res.df2 = res.df2 %>% filter(measures %in% c('TPR','FDR') & !(methods %in% delete))
thw =theme(axis.text.x = element_text(color="black", size =26),
           axis.text.y = element_text(color="black", size = 26),
           axis.title = element_text(color="black", size = 30),
           strip.text = element_text(size = 24),
           strip.background = element_rect(fill="white",color = "black", size = 1),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_rect(colour = "black",size= 1),
           legend.title=element_text(size=30),
           legend.text = element_text(size=30),
           plot.title = element_text(size=22))
head(res.df2)
size = 6
table(res2$methods)
newlevels = c(as.character(unique(res2$methods)[grep('edgeR',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('\\+DESeq2$',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('metagenomeSeq',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('eBay',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('Aldex',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('LDM',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('RAIDA',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('mbzinb',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('^TSS',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('DACOMP',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('Rarefy',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('ANCOM-BC',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('Beta-binomial',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('^GLM',unique(res2$methods))]),'DESeq2')
res2$methods = factor(res2$methods, levels = newlevels)
res2.fdr = res2[res2$measures =='FDR',]

m = NULL
for(i in unique(res2.fdr$methods)){
    df = res2.fdr %>% filter(methods ==i)
    md = aggregate(value~grp, data = df, function(x) median(x))
    mn = aggregate(value~grp, data = df, function(x) mean(x))
    
    if(md[md$grp =='No','value'] < md[md$grp =='Yes','value']){
      cat(i,' does not increase in median! \n')
    }
    
    if(mn[mn$grp =='No','value'] < mn[mn$grp =='Yes','value']){
      cat(i,' does not increase in mean! \n')
    }
    try({
      
    t = t.test(df$value ~ df$grp)
    if(t$p.value > 0.05){
      cat('No significant difference found in',i,'\n')
    }else{
      m = c(m, i)
    }
  })
}

df = res2.fdr %>% filter(!(methods %in% delete))
df$Improve = 'No improve'
df[df$methods %in% m, 'Improve'] = 'Improve' 
p1 = ggplot(df, aes(x =grp, y = value,  fill = methods)) +
  theme_bw() +
  geom_violin()+
  geom_boxplot(width=0.05,outlier.size = 0.01, fill = 'white') +
  scale_fill_manual(values = cols) +
  facet_wrap(Improve~methods, nrow = 2, scales = 'free') +
  scale_y_continuous(expand = c(0.1, 0, 0, 0)) +
  xlab('') +
  ylab('FDR')+
  labs(color = "", fill = '') +
  guides(fill = FALSE)+thw+
  geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5)


ggsave(paste0("/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/",'Stool_C1_GeneralBoxplot_rarefyCompareFDR.pdf'), width =30, height = 18, dpi = 100)

res2.tpr = res2[res2$measures =='TPR',]
ggplot(res2.tpr, aes(x =grp, y = value,  fill = methods)) +
  theme_bw() +
  geom_violin()+
  geom_boxplot(width=0.1,outlier.size = 0.01, fill = 'white') +
  scale_fill_manual(values = cols) +
  facet_wrap(.~methods, nrow = 4, scales = 'free') +
  scale_y_continuous(expand = c(0.1, 0, 0, 0)) +
  xlab('') +
  ylab('Power')+
  labs(color = "", fill = '') +
  guides(fill = FALSE)+thw+
  geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5)

for(i in unique(res2.tpr$methods)){
  df = res2.tpr %>% filter(methods ==i)
  md = aggregate(value~grp, data = df, function(x) median(x))
  mn = aggregate(value~grp, data = df, function(x) mean(x))
  
  if(md[md$grp =='No','value'] < md[md$grp =='Yes','value']){
    cat(i,' does not increase in median! \n')
  }

  if(mn[mn$grp =='No','value'] < mn[mn$grp =='Yes','value']){
    cat(i,' does not increase in mean! \n')
  }
  t = t.test(df$value ~ df$grp)
  if(t$p.value > 0.05){
    cat('No significant difference found in',i,'\n')
  }
}

ggsave(paste0("/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/",'Stool_C1_GeneralBoxplot_rarefyCompareTPR.pdf'), width =30, height = 18, dpi = 100)

