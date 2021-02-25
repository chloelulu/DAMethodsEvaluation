source('~/Documents/Mayo_project/Code/DailyCode.R')
load('/Users/m216453/Documents/Mayo_Research/SemiSimulation/0925/sim/new/sim/allC1_Omnibus_winsorF_res.Rdata')
winsorF = melt(res)
load('/Users/m216453/Documents/Mayo_Research/SemiSimulation/0925/sim/new/sim/allC1_res.Rdata')
winsor = melt(res)
colnames(winsorF) <- colnames(winsor) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')

winsor = winsor %>% filter(methods =='Omnibus' & measures %in% c('FDR','TPR') & covariate.eff.means %in% c('L2','L3','L4')) %>%
  dplyr::select(c('diff.otu.modes','covariate.types','covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')) %>% droplevels()
winsorF = winsorF %>% filter(methods =='Omnibus' & measures %in% c('FDR','TPR') & covariate.eff.means %in% c('L2','L3','L4'))%>%
  dplyr::select(c('diff.otu.modes','covariate.types','covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters'))%>% droplevels()
levels(winsor$covariate.eff.means) <- levels(winsorF$covariate.eff.means) <- c('Weak effect','Moderate effect','Strong effect')

winsor_FDR  = winsor[winsor$measures =='FDR',]
winsor_FDR$value = winsor_FDR$value/0.05
winsor_noFDR  = winsor[winsor$measures !='FDR',]
winsor = rbind(winsor_FDR, winsor_noFDR)
winsor <- winsor[!is.na(winsor$value),]


winsorF_FDR  = winsorF[winsorF$measures =='FDR',]
winsorF_FDR$value = winsorF_FDR$value/0.05
winsorF_noFDR  = winsorF[winsorF$measures !='FDR',]
winsorF = rbind(winsorF_FDR, winsorF_noFDR)
winsorF <- winsorF[!is.na(winsorF$value),]

winsorF$methods = gsub('Omnibus','No winsorization',winsorF$methods)
winsor$methods = gsub('Omnibus','Winsorization',winsor$methods)

res2 = rbind(winsor, winsorF)
levels(res2$diff.otu.modes) = c("Abundant", "Rare")
levels(res2$diff.otu.pcts) = c("Low density", "Medium density","High density")




res.df2 <- data_summary(data= res2, formula = paste0('value ~ covariate.types +diff.otu.modes+covariate.eff.means+ diff.otu.pcts+ methods + measures'))

measure = 'TPR'
measure1 = 'FDR'
measure2 = 'na'
ylab = 'observed FDR/expected FDR(0.05)'
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
head(res.df2)
size = 6
cols1 = c('No winsorization'=brewer.pal(9,'Paired')[4],'Winsorization'=brewer.pal(9,'PuBuGn')[7])

for(covariate.type in c('continuous','binary')){
  obj1 <- ggplot(res.df2 %>% filter(measures == 'TPR' & covariate.types ==covariate.type) %>% droplevels(), aes(x = covariate.eff.means, y = value,  fill = methods)) +
    theme_bw() +
    scale_fill_manual(values = cols1) +
    geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7) +
    geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
    scale_y_continuous(limits = c(0,1),expand = c(0.1, 0, 0, 0)) +
    facet_grid(as.formula(grid.formula), scales = 'free_y')+
    labs(y = 'Power', x = '', color = "", fill = '') +
    thw+
    labs(title = covariate.type)
  obj2 <- ggplot(res.df2 %>% filter(measures == 'FDR' & covariate.types ==covariate.type), aes(x = covariate.eff.means, y = value,  fill = methods)) +
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
    labs(title = covariate.type)
  box1 <- ggplot(res2 %>% filter(measures == 'TPR' & covariate.types ==covariate.type) %>% droplevels(), aes(x = covariate.eff.means, y = value,  fill = methods)) +
    theme_bw() +
    scale_fill_manual(values = cols1) +
    stat_boxplot(geom = "errorbar", aes(width = 0.2), lwd = 0.2) +
    geom_boxplot(aes(width=0.8),outlier.size = 0.1, outlier.colour = 'grey30') +
    scale_y_continuous(limits = c(0,1),expand = c(0.1, 0, 0, 0)) +
    facet_grid(as.formula(grid.formula), scales = 'free_y')+
    labs(y = 'Power', x = '', color = "", fill = '') +
    thw+
    labs(title = covariate.type)
  box2 <- ggplot(res2 %>% filter(measures == 'FDR' & covariate.types ==covariate.type), aes(x = covariate.eff.means, y = value,  fill = methods)) +
    theme_bw() +
    stat_boxplot(geom = "errorbar", aes(width = 0.2), lwd = 0.2) +
    geom_boxplot(aes(width=0.8),outlier.size = 0.1, outlier.colour = 'grey30') +
    scale_fill_manual(values = cols1) +
    scale_y_continuous(expand = c(0.1, 0, 0, 0)) +
    facet_grid(as.formula(grid.formula), scales = 'free_y') +
    xlab('') +
    ylab(ylab)+
    labs(color = "", fill = '') +
    guides(fill = FALSE)+thw+
    geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5) +
    labs(title = covariate.type)
  p = ggarrange(box1, box2, common.legend = T)
  ggsave(paste0("/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/",'Stool_C1_',covariate.type,'_','boxplot_WinsorCompare.pdf'), width = 34, height = 20, dpi = 100)
  
  p = ggarrange(obj1, obj2, common.legend = T)
  ggsave(paste0("/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/",'Stool_C1_',covariate.type,'_','barplot_WinsorCompare.pdf'), width = 34, height = 20, dpi = 100)
}



ggplot(res2[res2$measures =='FDR',], aes(x =methods, y = value,  fill = methods)) +
  theme_bw() +
  geom_violin()+
  geom_boxplot(width=0.1,fill = 'white') +
  scale_fill_manual(values = cols1) +
  scale_y_continuous(expand = c(0.1, 0, 0, 0)) +
  xlab('') +
  ylab(ylab)+
  labs(color = "", fill = '') +
  guides(fill = FALSE)+thw+
  geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5)
ggsave(paste0("/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/",'Stool_C1_GeneralBoxplot_WinsorCompare.pdf'), width =8, height = 8, dpi = 100)

