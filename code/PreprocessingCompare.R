source('~/Documents/Mayo_project/Code/DailyCode.R')
load('/Users/m216453/Documents/Mayo_Research/SemiSimulation/0925/sim/new/sim/nofilter/sim/allC1_res.Rdata')
filterF = melt(res)
load('/Users/m216453/Documents/Mayo_Research/SemiSimulation/0925/sim/new/sim/allC1_res.Rdata')
filter = melt(res)
colnames(filterF) <- colnames(filter) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')

filter = filter %>% filter( measures %in% c('FDR','TPR') & covariate.eff.means %in% c('L2','L3','L4')) %>%
  dplyr::select(c('diff.otu.modes','covariate.types','covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')) %>% droplevels()
filterF = filterF %>% filter(measures %in% c('FDR','TPR') & covariate.eff.means %in% c('L2','L3','L4'))%>%
  dplyr::select(c('diff.otu.modes','covariate.types','covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters'))%>% droplevels()
levels(filter$covariate.eff.means) <- levels(filterF$covariate.eff.means) <- c('Weak effect','Moderate effect','Strong effect')

filter_FDR  = filter[filter$measures =='FDR',]
filter_FDR$value = filter_FDR$value/0.05
filter_noFDR  = filter[filter$measures !='FDR',]
filter = rbind(filter_FDR, filter_noFDR)
filter <- filter[!is.na(filter$value),] %>% mutate(grp = 'Yes')


filterF_FDR  = filterF[filterF$measures =='FDR',]
filterF_FDR$value = filterF_FDR$value/0.05
filterF_noFDR  = filterF[filterF$measures !='FDR',]
filterF = rbind(filterF_FDR, filterF_noFDR)
filterF <- filterF[!is.na(filterF$value),]%>% mutate(grp = 'No')


res2 = rbind(filter, filterF) %>% filter(covariate.types =='binary' & methods != 'CLRBC') %>% droplevels()
levels(res2$diff.otu.modes) = c("Abundant", "Rare")
levels(res2$diff.otu.pcts) = c("Low density", "Medium density","High density")

cols <- c('ZicoSeq'=brewer.pal(9,'Paired')[5], 'Omnibus' = brewer.pal(9,'PuBuGn')[7], 
          'GLM(quasipoisson)' = brewer.pal(9,'PuRd')[3],'CLRBC'=brewer.pal(9,'Set1')[8],
          'RAIDA'=brewer.pal(11,'BrBG')[7],
          'eBay(Wilcoxon)' = brewer.pal(8,'RdPu')[4],'eBay(t-test)' = brewer.pal(8,'RdPu')[6], 'ANCOM-BC' = brewer.pal(8,'Set2')[7],
          'LDM'=brewer.pal(9,'Set1')[1], 'DACOMP'=brewer.pal(8,'Dark2')[6],'Aldex2(Wilcoxon)'=brewer.pal(8,'YlGn')[3],'Aldex2(t-test)'=brewer.pal(8,'YlGn')[2],
          'Beta-binomial'=brewer.pal(11,'BrBG')[2],
          'DESeq2'=brewer.pal(9,'Greys')[7], 'Wrench+DESeq2'= brewer.pal(9,'Greys')[3], 'GMPR+DESeq2'=brewer.pal(9,'Greys')[5],
          'TSS+Wilcoxon'=brewer.pal(9,'Purples')[7],'Rarefy+Wilcoxon'= brewer.pal(9,'Set1')[5], 
          'TSS+Spearman'=brewer.pal(9,'Purples')[7],'Rarefy+Spearman'=brewer.pal(9,'Set1')[5],
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

res.df2 <- data_summary(data= res2, formula = paste0('value ~ covariate.types +diff.otu.modes+covariate.eff.means+ diff.otu.pcts+ methods + measures+ label + legend'))

measure = 'TPR'
measure1 = 'FDR'
measure2 = 'na'
ylab = 'observed FDR/expected FDR(0.05)'
grid.formula = 'diff.otu.pcts ~ diff.otu.modes'
res.df2 = res.df2 %>% filter(measures %in% c('TPR','FDR'))
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
newlevels = c(as.character(unique(res2$methods)[grep('edgeR',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('DESeq2',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('metagenomeSeq',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('eBay',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('Aldex',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('LDM',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('RAIDA',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('Omnibus',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('DACOMP',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('Rarefy',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('Beta-binomial',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('^TSS',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('^GLM',unique(res2$methods))]),
              as.character(unique(res2$methods)[grep('ANCOM-BC',unique(res2$methods))]))
res2$methods = factor(res2$methods, levels = newlevels)
ggplot(res2[res2$measures =='FDR',], aes(x =grp, y = value,  fill = methods)) +
  theme_bw() +
  geom_violin()+
  geom_boxplot(width=0.1,fill = 'white', outlier.size = 0.01) +
  scale_fill_manual(values = cols) +
  facet_wrap(.~methods, nrow = 4, scales = 'free') +
  scale_y_continuous(expand = c(0.1, 0, 0, 0)) +
  xlab('') +
  ylab(ylab)+
  labs(color = "", fill = '') +
  guides(fill = FALSE)+thw+
  geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5)
ggsave(paste0("/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/",'Stool_C1_GeneralBoxplot_filterCompare.pdf'), width =30, height = 18, dpi = 100)

