pkg = c('dplyr','tidyr','reshape','RColorBrewer','ggplot2','ggpubr','stringi','scales','reshape2','magick','kableExtra','stringr')
suppressPackageStartupMessages(sapply(pkg, require, character = T))
setwd('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/')
source('/Users/m216453/Documents/Mayo_project/Code/DailyCode.R')

for(name in c('allD0')){
  if(name =='allC0'){
    name1 ='Stool'
  }else{
    name1='Vaginal'
  }
  # load(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/',name,'_res.Rdata'))
  # # load(paste0('~/Documents/Mayo_Research/SemiSimulation/0925/sim/new/sim/',name,'_res.Rdata'))
  # res0 <- melt(res)
  # colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs',
  #                     'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams',
  #                     'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
  getwd()
  # save(res0,file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/',name,'_melt.Rdata'))
  load(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/',name,'_melt.Rdata'))
  for(covariate.type in c('binary','continuous')){
    res1 = res0 %>% filter((nOTUs %in% c('nOTU_L1','nOTU_L3','nOTU_L5'))) %>% filter(nSams %in% c('nSam_L1','nSam_L2','nSam_L4') & measures =='FP') %>% droplevels() %>% 
      dplyr::select(nOTUs,nSams, covariate.types, depth.conf.factors, methods,value) %>% na.omit() 
    
    res1$depth.conf.factors = gsub('DL3','Depth confounding',res1$depth.conf.factors)
    res1$depth.conf.factors = gsub('none','None',res1$depth.conf.factors)
    
    levels(res1$nSams) <- c('sample=50','sample=100','sample=200')
    levels(res1$nOTUs) <- c('OTU=50','OTU=200','OTU=500')
    levels(res1$depth.conf.factors) <- c('None','Depth confounding')
    
    na = res1[is.na(res1$value),]
    res1$methods = as.character(res1$methods)
    res1$methods[res1$methods =='Wilcox.gmpr' & res1$covariate.types =='binary'] = 'GMPR+Wilcoxon'
    res1$methods[res1$methods =='Wilcox.Wrench' & res1$covariate.types =='binary'] = 'Wrench+Wilcoxon'
    res1$methods[res1$methods =='Wilcox' & res1$covariate.types =='binary'] = 'TSS+Wilcoxon'
    res1$methods[res1$methods =='Wilcox.gmpr' & res1$covariate.types =='continuous'] = 'GMPR+Spearman'
    res1$methods[res1$methods =='Wilcox' & res1$covariate.types =='continuous'] = 'TSS+Spearman'
    res1$methods[res1$methods =='Rarefyttest' & res1$covariate.types =='binary'] = 'Rarefy+t-test'
    res1$methods[res1$methods =='ttest' & res1$covariate.types =='binary'] = 'TSS+t-test'
    res1$methods[res1$methods =='ttest.gmpr' & res1$covariate.types =='binary'] = 'GMPR+t-test'
    res1$methods[res1$methods =='ttest.Wrench' & res1$covariate.types =='binary'] = 'Wrench+t-test'
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
    
    res1 = res1 %>% filter(!(methods  %in% c('LinDA','GMPR+Spearman','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon')))
    sub0 = c("Aldex2(Wilcoxon)","ANCOM-BC","Beta-binomial","DACOMP","eBay(t-test)",
             "GLM(quasipoisson)","GMPR+DESeq2","GMPR+edgeR","LDM","mbzinb",
             'DESeq2','Wrench+DESeq2','edgeR', 'Wrench+edgeR','Wrench+metagenomeSeq',
             'eBay(Wilcoxon)','Aldex2(t-test)',"RAIDA","Rarefy+t-test","Rarefy+Wilcoxon",
             "TSS+t-test","TSS+Wilcoxon","metagenomeSeq") # for binary
    sub0_sup = c('RioNorm2') # for supplementary binary
    sub1 = c('DESeq2','GMPR+DESeq2','edgeR', 'GMPR+edgeR','eBay(t-test)',
             'ANCOM-BC','DACOMP','LDM','GLM(quasipoisson)','Beta-binomial',
             'Rarefy+Spearman','TSS+Spearman') # for continuous
    
    
    if(covariate.type =='continuous'){
      sub = sub1
      res1 = res1 %>% filter(methods %in% sub1)
    }else{
      sub =sub0
      res1 = res1 %>% filter(methods %in% c(sub0,sub0_sup))
    }
    
    ## For kable of false discover FPs
    x1 = res1 %>%filter(covariate.types==covariate.type) %>% dplyr::select(depth.conf.factors,nOTUs, nSams, methods, value) #%>% filter(nOTUs =='OTU=50') %>% dplyr::select(-nOTUs)
    x1 = aggregate(value~., x1, function(x) mean(x))
    x1 = x1 %>% filter((methods %in% sub))
    x1$nOTUs = factor(x1$nOTUs, levels = c('OTU=50','OTU=200','OTU=500'))
    x1$nSams = factor(x1$nSams, levels = c('sample=50','sample=100','sample=200'))
    x1$depth.conf.factors = factor(x1$depth.conf.factors, levels = c('None','Depth confounding'))
    
    ## plot barplot for FPs
    m = x1
    letters1 = c(letters,rev(toupper(letters)))
    for(i in 1:length(unique(m$methods))){
      m[(m$methods ==unique(m$methods)[i]),'label'] = letters1[i]
    }
    m$legend = paste0(m$label,':',m$methods)
    cols <- c('ZicoSeq'=brewer.pal(9,'Paired')[5], 'Omnibus' = brewer.pal(9,'PuBuGn')[7], 'mbzinb' = brewer.pal(9,'PuBuGn')[8], 
              'GLM(quasipoisson)' = brewer.pal(9,'PuRd')[3],'CLRBC'=brewer.pal(9,'Set1')[8],
              'RAIDA'=brewer.pal(11,'BrBG')[7],'LinDA'=brewer.pal(9,'Set1')[8],
              'eBay(Wilcoxon)' = brewer.pal(8,'RdPu')[4],'eBay(t-test)' = brewer.pal(8,'RdPu')[6], 'ANCOM-BC' = brewer.pal(8,'Set2')[7],
              'LDM'=brewer.pal(9,'Set1')[1], 'DACOMP'=brewer.pal(8,'Dark2')[6],'Aldex2(Wilcoxon)'=brewer.pal(8,'YlGn')[3],'Aldex2(t-test)'=brewer.pal(8,'YlGn')[2],
              'Beta-binomial'=brewer.pal(11,'BrBG')[2],
              'DESeq2'=brewer.pal(9,'Greys')[7], 'Wrench+DESeq2'= brewer.pal(9,'Greys')[3], 'GMPR+DESeq2'=brewer.pal(9,'Greys')[5],
              'TSS+Wilcoxon'=brewer.pal(9,'Purples')[7],'Rarefy+Wilcoxon'= brewer.pal(9,'Set1')[5], 'GMPR+Wilcoxon'= brewer.pal(9,'Purples')[5], 'Wrench+Wilcoxon'= brewer.pal(9,'Purples')[3], 
              'TSS+Spearman'=brewer.pal(9,'Purples')[7],'Rarefy+Spearman'=brewer.pal(9,'Set1')[5], 'GMPR+Spearman'= brewer.pal(9,'Purples')[5], 
              'TSS+t-test'=brewer.pal(9,'GnBu')[7],'Rarefy+t-test'= brewer.pal(9,'GnBu')[5], 'GMPR+t-test'= brewer.pal(9,'GnBu')[5], 'Wrench+t-test'= brewer.pal(9,'GnBu')[3], 
              'edgeR'=brewer.pal(9,'Greens')[7], 'Wrench+edgeR'=brewer.pal(9,'Greens')[3], 'GMPR+edgeR'=brewer.pal(9,'Greens')[5], 
              'metagenomeSeq'=brewer.pal(9,'Blues')[7], 'Wrench+metagenomeSeq'=brewer.pal(9,'Blues')[5])
    
    cols1 = cols
    
    s = unique(m$legend)
    for(i in 1:length(cols1)){
      if(names(cols1)[i] %in% gsub('.*:','',s)){
        names(cols1)[i] = s[which((names(cols1)[i] == gsub('.*:','',s)))]
      }else{
        cat('This legend does not exist! \n')
      }
    }
    # save(m,file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/',name1,covariate.type,'BiasDAs.Rdata'))
  head(m)
    size = 6
    # if(covariate.type =='binary'){
    #   df = m %>% filter(covariate.types==covariate.type& (methods %in% c(sub0,sub0_sup))) %>% droplevels()
    # }else{
    #   df = m %>% filter(covariate.types==covariate.type& (methods %in% c(sub1)))
    # }
    m[m$value > 5,'shape'] = 'uncontrol'
    m[m$value <= 1,'shape'] = 'fdr_0.05'
    m[m$value < 2 & m$value>1,'shape'] = 'fdr_0.1'
    m[m$value <= 5 & m$value>2,'shape'] = 'fdr_0.2'
    shapes = c('fdr_0.05' =8, 'uncontrol' = 19, 'fdr_0.1' = 17,'fdr_0.2' = 18)
    
    thw =theme(axis.text.x = element_text(color="black", size =26),
               axis.text.y = element_text(color="black", size = 18),
               axis.title = element_text(color="black", size = 28),
               strip.text = element_text(size = 28),
               strip.background = element_rect(fill="white",color = "black", size = 1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(colour = "black",size= 1),
               legend.position = 'top',
               legend.title=element_text(size=26),
               legend.text = element_text(size=26),
               plot.title = element_text(size=36),
               plot.caption = element_text(size = 20))
    p <- ggplot(m,aes(x = nOTUs, y = value, fill = legend)) + # fill = legend shhould be here, or the geom_text label will have same position of label
      geom_bar(position = position_dodge(width = 0.8), stat="identity", width = 0.1) +
      geom_point(aes(color = legend, shape = shape),  size = 5, position = position_dodge(width = 0.8)) +
      geom_text(aes(label=label), position=position_dodge(width=0.8),size=5, vjust = -0.7) +
      facet_grid(nSams~depth.conf.factors, scales = 'free_y') +
      scale_fill_manual(values = c(rep('grey50',length(m$methods %>% unique())))) +
      scale_color_manual(values = cols1) +
      scale_shape_manual(values = shapes) +
      scale_y_continuous(trans = sqrt_trans(),
                         breaks = trans_breaks("sqrt", function(x) x^2),
                         labels = trans_format("sqrt", math_format(.x^2)))+
      
      # scale_y_continuous(breaks=c(1, 2, 5, max(df$ct))) +
      theme_bw() +thw +
      labs(y = 'FPs', x = 'OTU number', fill = "", shape = "", color = "",
           title = paste0("Default: no differential taxa (", covariate.type,"-",name1,')')) +
      guides(shape = FALSE, fill = FALSE)
    
    ggsave(paste0('result/SimulationEvaluation/',name1,covariate.type,'_','FPs_barplot_Null.pdf'),  width = 26, height = 20, dpi = 100)
    
    ## for kabel FPs
    x2 = x1 %>% unite('grp',c('depth.conf.factors','nOTUs','nSams'))%>% filter(methods !='RioNorm2')
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
    x4 = apply(x4, 2, function(x) ifelse(x<1, '***',ifelse(x<2, '**', ifelse(x<5,'*','x'))))
    x4 = as.data.frame(x4)
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
      for(i in 1:nrow(data)){
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
      add_header_above(c("Vaginal FPs" = 1, "None" =9," " = 1, "Depth confounding" = 9," " = 1," " = 1), bold = T) %>%
      row_spec(0,bold=TRUE)  %>%
      column_spec(2, background= count.ct(data = x4,2),bold = T) %>%
      column_spec(3, background= count.ct(data = x4,3),bold = T) %>%
      column_spec(4, background= count.ct(data = x4,4),bold = T) %>%
      column_spec(5, background= count.ct(data = x4,5),bold = T) %>% 
      column_spec(6, background= count.ct(data = x4,6),bold = T) %>% 
      column_spec(7, background= count.ct(data = x4,7),bold = T) %>% 
      column_spec(8, background= count.ct(data = x4,8),bold = T) %>% 
      column_spec(9, background= count.ct(data = x4,9),bold = T) %>% 
      column_spec(10, background= count.ct(data = x4,10),bold = T) %>% 
      column_spec(12, background= count.ct(data = x4,12),bold = T)%>% 
      column_spec(13, background= count.ct(data = x4,13),bold = T)%>% 
      column_spec(14, background= count.ct(data = x4,14),bold = T)%>% 
      column_spec(15, background= count.ct(data = x4,15),bold = T)%>% 
      column_spec(16, background= count.ct(data = x4,16),bold = T)%>% 
      column_spec(17, background= count.ct(data = x4,17),bold = T)%>% 
      column_spec(18, background= count.ct(data = x4,18),bold = T)%>%
      column_spec(19, background= count.ct(data = x4,19),bold = T) %>% 
      column_spec(20, background= count.ct(data = x4,20),bold = T) %>%
      # column_spec(21, background= count.ct(data = x4,21),bold = T) %>%
      # column_spec(22, background= count.ct(data = x4,22),bold = T) %>%
      column_spec(11, color = count.score(data = x4, j = 11,color1 = 'white', color2 = 'black'),
                  background = count.score(data = x4, j = 11,color1 = brewer.pal(8,'Dark2')[1], color2 = 'white')) %>%
      column_spec(21, color = count.score(data = x4, j = 21,color1 = 'white', color2 = 'black'),
                  background = count.score(data = x4, j = 21,color1 = brewer.pal(8,'Dark2')[1], color2 = 'white')) %>%
      column_spec(22, color = count.score(data = x4, j = 22,color1 = 'white', color2 = 'black', value = 54),
                  background = count.score(data = x4, j = 22,color1 = brewer.pal(8,'Dark2')[3], color2 = 'white', value = 54),bold = T) %>%
      save_kable(file = paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/SimulationEvaluation/',covariate.type,'_',name1,'_FPs_Null.pdf'))
    
    ## For FDIs
    
    head(res1)
    df = res1 %>%filter(covariate.types==covariate.type) %>% dplyr::select(depth.conf.factors,nOTUs, nSams, methods, value) 
    df[df$value>0,'ct'] = 1
    df[df$value==0,'ct'] = 0
    df = df %>% dplyr::select(-value) 
    formula = paste0('ct~ depth.conf.factors + nOTUs + nSams + methods')
    grid.formula = '. ~ nSams'
    m <- aggregate(as.formula(formula), df, function(x) sum(x[!is.na(x)])) %>% mutate(ct = ct/1000)
    kb = m
    head(kb)

    kb[kb$ct > 0.1,'shape'] = 'uncontrol'
    kb[kb$ct <= 0.05,'shape'] = 'fdr_0.05'
    kb[kb$ct <= 0.1 & kb$ct>0.05,'shape'] = 'fdr_0.1'
    kb[kb$ct <= 0.2 & kb$ct>0.1,'shape'] = 'fdr_0.2'
    shapes = c('fdr_0.05' =8, 'uncontrol' = 19, 'fdr_0.1' = 17,'fdr_0.2' = 18)
    
    
    head(kb)
    m = kb
    letters1 = c(letters,rev(toupper(letters)))
    for(i in 1:length(unique(m$methods))){
      m[(m$methods ==unique(m$methods)[i]),'label'] = letters1[i]
    }
    m$legend = paste0(m$label,':',m$methods)
    cols <- c('ZicoSeq'=brewer.pal(9,'Paired')[5], 'Omnibus' = brewer.pal(9,'PuBuGn')[7], 'mbzinb' = brewer.pal(9,'PuBuGn')[8], 
              'GLM(quasipoisson)' = brewer.pal(9,'PuRd')[3],'CLRBC'=brewer.pal(9,'Set1')[8],
              'RAIDA'=brewer.pal(11,'BrBG')[7],'LinDA'=brewer.pal(9,'Set1')[8],
              'eBay(Wilcoxon)' = brewer.pal(8,'RdPu')[4],'eBay(t-test)' = brewer.pal(8,'RdPu')[6], 'ANCOM-BC' = brewer.pal(8,'Set2')[7],
              'LDM'=brewer.pal(9,'Set1')[1], 'DACOMP'=brewer.pal(8,'Dark2')[6],'Aldex2(Wilcoxon)'=brewer.pal(8,'YlGn')[3],'Aldex2(t-test)'=brewer.pal(8,'YlGn')[2],
              'Beta-binomial'=brewer.pal(11,'BrBG')[2],
              'DESeq2'=brewer.pal(9,'Greys')[7], 'Wrench+DESeq2'= brewer.pal(9,'Greys')[3], 'GMPR+DESeq2'=brewer.pal(9,'Greys')[5],
              'TSS+Wilcoxon'=brewer.pal(9,'Purples')[7],'Rarefy+Wilcoxon'= brewer.pal(9,'Set1')[5], 'GMPR+Wilcoxon'= brewer.pal(9,'Purples')[5], 'Wrench+Wilcoxon'= brewer.pal(9,'Purples')[3], 
              'TSS+Spearman'=brewer.pal(9,'Purples')[7],'Rarefy+Spearman'=brewer.pal(9,'Set1')[5], 'GMPR+Spearman'= brewer.pal(9,'Purples')[5], 
              'TSS+t-test'=brewer.pal(9,'GnBu')[7],'Rarefy+t-test'= brewer.pal(9,'GnBu')[5], 'GMPR+t-test'= brewer.pal(9,'GnBu')[5], 'Wrench+t-test'= brewer.pal(9,'GnBu')[3], 
              'edgeR'=brewer.pal(9,'Greens')[7], 'Wrench+edgeR'=brewer.pal(9,'Greens')[3], 'GMPR+edgeR'=brewer.pal(9,'Greens')[5], 
              'metagenomeSeq'=brewer.pal(9,'Blues')[7], 'Wrench+metagenomeSeq'=brewer.pal(9,'Blues')[5])
    
    cols1 = cols
    
    s = unique(m$legend)
    for(i in 1:length(cols1)){
      if(names(cols1)[i] %in% gsub('.*:','',s)){
        names(cols1)[i] = s[which((names(cols1)[i] == gsub('.*:','',s)))]
      }else{
        cat('This legend does not exist! \n')
      }
    }
    # save(m,file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/',name1,covariate.type,'BiasDAs.Rdata'))
    head(m)
    
    size = 6

    p0 = ggplot(m,aes(x = nOTUs, y = ct, fill = legend)) + # fill = legend shhould be here, or the geom_text label will have same position of label
      geom_bar(position = position_dodge(width = 0.8), stat="identity", width = 0.1) +
      geom_point(aes(color = legend, shape = shape),  size = 5, position = position_dodge(width = 0.8)) +
      geom_text(aes(label=label), position=position_dodge(width=0.8),size=5, vjust = -0.7) + 
      facet_grid(nSams~depth.conf.factors, scales = 'free_y') +
      scale_fill_manual(values = c(rep('grey50',length(df$methods %>% unique())))) +
      scale_color_manual(values = cols1) +
      scale_shape_manual(values = shapes) +
      scale_y_continuous(limits = c(0,1),breaks=c(0.05, 0.1, 0.2, 1)) +
      theme_bw() +thw +
      labs(y = 'FDIs', x = 'OTU number', fill = "", shape = "", color = "",
           title = paste0("Default: no differential taxa (", covariate.type,"-",name1,')')) +
      geom_hline(yintercept = 0.05, colour = '#238823', linetype = 'dashed', size =  0.5) + # green
      geom_hline(yintercept = 0.1, colour = '#ffbf00', linetype = 'dashed', size =  0.5) + # yellow
      geom_hline(yintercept = 0.2, colour = '#e03531', linetype = 'dashed', size =  0.5) + # red
      guides(shape = FALSE, fill = FALSE)
    ggsave(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/SimulationEvaluation/',name1,'_',covariate.type,'_FDIs_barplot_Null.pdf'), width = 25, height =20)
    
    
    
    
    ## kable FDIs
    kb = m
    kb[kb$ct > 0.1,'shape'] = 'x'
    kb[kb$ct <= 0.05,'shape'] = '***'
    kb[kb$ct <= 0.1 & kb$ct>0.05,'shape'] = '**'
    kb[kb$ct <= 0.2 & kb$ct>0.1,'shape'] = '*'
    
    x1 = kb %>% dplyr::select(depth.conf.factors,nOTUs, nSams, methods, shape) %>% filter(methods !='RioNorm2')
    x1$nOTUs = factor(x1$nOTUs, levels = c('OTU=50','OTU=200','OTU=500'))
    x1$nSams = factor(x1$nSams, levels = c('sample=50','sample=100','sample=200'))
    x1$depth.conf.factors = factor(x1$depth.conf.factors, levels = c('None','Depth confounding'))
    x2 = x1 %>% unite('grp',c('depth.conf.factors','nOTUs','nSams'))
    x3 = x2 %>% spread(c('grp'), shape)
    head(x3)
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
      for(i in 1:nrow(data)){
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
      add_header_above(c("Vaginal FDIs" = 1, "None" =9," " = 1, "Depth confounding" = 9," " = 1," " = 1), bold = T) %>%
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
                  background = count.score(data = x4, j = 22,color1 = brewer.pal(8,'Dark2')[3], color2 = 'white', value = 54),bold = T) %>%
      save_kable(file = paste0(getwd(),'/result/SimulationEvaluation/',covariate.type,'_',name1,'_FDI_Null.pdf'))
    # as_image(width = 7, height = 3, file = paste0(getwd(),'/plot/',covariate.type,'_',name1,'.pdf')) 

  }
  
}





for(name in c('allC0')){
  if(name =='allC0'){
    name1 ='Stool'
  }else{
    name1='Vaginal'
  }
  # load(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/',name,'_res.Rdata'))
  # # load(paste0('~/Documents/Mayo_Research/SemiSimulation/0925/sim/new/sim/',name,'_res.Rdata'))
  # res0 <- melt(res)
  # colnames(res0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs',
  #                     'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams',
  #                     'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
  getwd()
  # save(res0,file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/',name,'_melt.Rdata'))
  load(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/',name,'_melt.Rdata'))
  for(covariate.type in c('binary','continuous')){
    res1 = res0 %>% filter((nOTUs %in% c('nOTU_L1','nOTU_L3','nOTU_L5'))) %>% filter(nSams %in% c('nSam_L1','nSam_L2','nSam_L4') & measures =='FP') %>% droplevels() %>% 
      dplyr::select(nOTUs,nSams, covariate.types, depth.conf.factors, methods,value) %>% na.omit() 
    
    res1$depth.conf.factors = gsub('DL3','Depth confounding',res1$depth.conf.factors)
    res1$depth.conf.factors = gsub('none','None',res1$depth.conf.factors)
    
    levels(res1$nSams) <- c('sample=50','sample=100','sample=200')
    levels(res1$nOTUs) <- c('OTU=50','OTU=200','OTU=500')
    levels(res1$depth.conf.factors) <- c('None','Depth confounding')
    
    na = res1[is.na(res1$value),]
    res1$methods = as.character(res1$methods)
    res1$methods[res1$methods =='Wilcox.gmpr' & res1$covariate.types =='binary'] = 'GMPR+Wilcoxon'
    res1$methods[res1$methods =='Wilcox.Wrench' & res1$covariate.types =='binary'] = 'Wrench+Wilcoxon'
    res1$methods[res1$methods =='Wilcox' & res1$covariate.types =='binary'] = 'TSS+Wilcoxon'
    res1$methods[res1$methods =='Wilcox.gmpr' & res1$covariate.types =='continuous'] = 'GMPR+Spearman'
    res1$methods[res1$methods =='Wilcox' & res1$covariate.types =='continuous'] = 'TSS+Spearman'
    res1$methods[res1$methods =='Rarefyttest' & res1$covariate.types =='binary'] = 'Rarefy+t-test'
    res1$methods[res1$methods =='ttest' & res1$covariate.types =='binary'] = 'TSS+t-test'
    res1$methods[res1$methods =='ttest.gmpr' & res1$covariate.types =='binary'] = 'GMPR+t-test'
    res1$methods[res1$methods =='ttest.Wrench' & res1$covariate.types =='binary'] = 'Wrench+t-test'
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
    
    res1 = res1 %>% filter(!(methods  %in% c('LinDA','GMPR+Spearman','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon')))
    sub0 = c("Aldex2(Wilcoxon)","ANCOM-BC","Beta-binomial","DACOMP","eBay(t-test)",
             "GLM(quasipoisson)","GMPR+DESeq2","GMPR+edgeR","LDM","mbzinb",
             'DESeq2','Wrench+DESeq2','edgeR', 'Wrench+edgeR','Wrench+metagenomeSeq',
             'eBay(Wilcoxon)','Aldex2(t-test)',"RAIDA","Rarefy+t-test","Rarefy+Wilcoxon",
             "TSS+t-test","TSS+Wilcoxon","metagenomeSeq") # for binary
    sub0_sup = c('RioNorm2') # for supplementary binary
    sub1 = c('DESeq2','GMPR+DESeq2','edgeR', 'GMPR+edgeR','eBay(t-test)',
             'ANCOM-BC','DACOMP','LDM','GLM(quasipoisson)','Beta-binomial',
             'Rarefy+Spearman','TSS+Spearman') # for continuous
    
    
    if(covariate.type =='continuous'){
      sub = sub1
      res1 = res1 %>% filter(methods %in% sub1)
    }else{
      sub =sub0
      res1 = res1 %>% filter(methods %in% c(sub0,sub0_sup))
    }
    
    ## For kable of false discover FPs
    x1 = res1 %>%filter(covariate.types==covariate.type) %>% dplyr::select(depth.conf.factors,nOTUs, nSams, methods, value) #%>% filter(nOTUs =='OTU=50') %>% dplyr::select(-nOTUs)
    x1 = aggregate(value~., x1, function(x) mean(x))
    x1 = x1 %>% filter((methods %in% sub))
    x1$nOTUs = factor(x1$nOTUs, levels = c('OTU=50','OTU=200','OTU=500'))
    x1$nSams = factor(x1$nSams, levels = c('sample=50','sample=100','sample=200'))
    x1$depth.conf.factors = factor(x1$depth.conf.factors, levels = c('None','Depth confounding'))
    
    ## plot barplot for FPs
    m = x1
    letters1 = c(letters,rev(toupper(letters)))
    for(i in 1:length(unique(m$methods))){
      m[(m$methods ==unique(m$methods)[i]),'label'] = letters1[i]
    }
    m$legend = paste0(m$label,':',m$methods)
    cols <- c('ZicoSeq'=brewer.pal(9,'Paired')[5], 'Omnibus' = brewer.pal(9,'PuBuGn')[7], 'mbzinb' = brewer.pal(9,'PuBuGn')[8], 
              'GLM(quasipoisson)' = brewer.pal(9,'PuRd')[3],'CLRBC'=brewer.pal(9,'Set1')[8],
              'RAIDA'=brewer.pal(11,'BrBG')[7],'LinDA'=brewer.pal(9,'Set1')[8],
              'eBay(Wilcoxon)' = brewer.pal(8,'RdPu')[4],'eBay(t-test)' = brewer.pal(8,'RdPu')[6], 'ANCOM-BC' = brewer.pal(8,'Set2')[7],
              'LDM'=brewer.pal(9,'Set1')[1], 'DACOMP'=brewer.pal(8,'Dark2')[6],'Aldex2(Wilcoxon)'=brewer.pal(8,'YlGn')[3],'Aldex2(t-test)'=brewer.pal(8,'YlGn')[2],
              'Beta-binomial'=brewer.pal(11,'BrBG')[2],
              'DESeq2'=brewer.pal(9,'Greys')[7], 'Wrench+DESeq2'= brewer.pal(9,'Greys')[3], 'GMPR+DESeq2'=brewer.pal(9,'Greys')[5],
              'TSS+Wilcoxon'=brewer.pal(9,'Purples')[7],'Rarefy+Wilcoxon'= brewer.pal(9,'Set1')[5], 'GMPR+Wilcoxon'= brewer.pal(9,'Purples')[5], 'Wrench+Wilcoxon'= brewer.pal(9,'Purples')[3], 
              'TSS+Spearman'=brewer.pal(9,'Purples')[7],'Rarefy+Spearman'=brewer.pal(9,'Set1')[5], 'GMPR+Spearman'= brewer.pal(9,'Purples')[5], 
              'TSS+t-test'=brewer.pal(9,'GnBu')[7],'Rarefy+t-test'= brewer.pal(9,'GnBu')[5], 'GMPR+t-test'= brewer.pal(9,'GnBu')[5], 'Wrench+t-test'= brewer.pal(9,'GnBu')[3], 
              'edgeR'=brewer.pal(9,'Greens')[7], 'Wrench+edgeR'=brewer.pal(9,'Greens')[3], 'GMPR+edgeR'=brewer.pal(9,'Greens')[5], 
              'metagenomeSeq'=brewer.pal(9,'Blues')[7], 'Wrench+metagenomeSeq'=brewer.pal(9,'Blues')[5])
    
    cols1 = cols
    
    s = unique(m$legend)
    for(i in 1:length(cols1)){
      if(names(cols1)[i] %in% gsub('.*:','',s)){
        names(cols1)[i] = s[which((names(cols1)[i] == gsub('.*:','',s)))]
      }else{
        cat('This legend does not exist! \n')
      }
    }
    # save(m,file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/',name1,covariate.type,'BiasDAs.Rdata'))
    head(m)
    size = 6
    # if(covariate.type =='binary'){
    #   df = m %>% filter(covariate.types==covariate.type& (methods %in% c(sub0,sub0_sup))) %>% droplevels()
    # }else{
    #   df = m %>% filter(covariate.types==covariate.type& (methods %in% c(sub1)))
    # }
    m[m$value > 5,'shape'] = 'uncontrol'
    m[m$value <= 1,'shape'] = 'fdr_0.05'
    m[m$value < 2 & m$value>1,'shape'] = 'fdr_0.1'
    m[m$value <= 5 & m$value>2,'shape'] = 'fdr_0.2'
    shapes = c('fdr_0.05' =8, 'uncontrol' = 19, 'fdr_0.1' = 17,'fdr_0.2' = 18)
    
    thw =theme(axis.text.x = element_text(color="black", size =26),
               axis.text.y = element_text(color="black", size = 18),
               axis.title = element_text(color="black", size = 28),
               strip.text = element_text(size = 28),
               strip.background = element_rect(fill="white",color = "black", size = 1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(colour = "black",size= 1),
               legend.position = 'top',
               legend.title=element_text(size=26),
               legend.text = element_text(size=26),
               plot.title = element_text(size=36),
               plot.caption = element_text(size = 20))
    p <- ggplot(m,aes(x = nOTUs, y = value, fill = legend)) + # fill = legend shhould be here, or the geom_text label will have same position of label
      geom_bar(position = position_dodge(width = 0.8), stat="identity", width = 0.1) +
      geom_point(aes(color = legend, shape = shape),  size = 5, position = position_dodge(width = 0.8)) +
      geom_text(aes(label=label), position=position_dodge(width=0.8),size=5, vjust = -0.7) +
      facet_grid(nSams~depth.conf.factors, scales = 'free_y') +
      scale_fill_manual(values = c(rep('grey50',length(m$methods %>% unique())))) +
      scale_color_manual(values = cols1) +
      scale_shape_manual(values = shapes) +
      scale_y_continuous(trans = sqrt_trans(),
                         breaks = trans_breaks("sqrt", function(x) x^2),
                         labels = trans_format("sqrt", math_format(.x^2)))+
      
      # scale_y_continuous(breaks=c(1, 2, 5, max(df$ct))) +
      theme_bw() +thw +
      labs(y = 'FPs', x = 'OTU number', fill = "", shape = "", color = "",
           title = paste0("Default: no differential taxa (", covariate.type,"-",name1,')')) +
      # geom_hline(yintercept = 1, colour = '#238823', linetype = 'dashed', size =  0.5) + # green
      # geom_hline(yintercept = 2, colour = '#ffbf00', linetype = 'dashed', size =  0.5) + # yellow
      # geom_hline(yintercept = 5, colour = '#e03531', linetype = 'dashed', size =  0.5) + # red
      guides(shape = FALSE, fill = FALSE)
    
    ggsave(paste0('result/SimulationEvaluation/',name1,covariate.type,'_','FPs_barplot_Null.pdf'),  width = 26, height = 20, dpi = 100)
    
    ## for kabel FPs
    x2 = x1 %>% unite('grp',c('depth.conf.factors','nOTUs','nSams'))%>% filter(methods !='RioNorm2')
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
    x4 = apply(x4, 2, function(x) ifelse(x<1, '***',ifelse(x<2, '**', ifelse(x<5,'*','x'))))
    x4 = as.data.frame(x4)
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
      for(i in 1:nrow(data)){
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
      add_header_above(c("Stool FPs" = 1, "None" =9," " = 1, "Depth confounding" = 9," " = 1," " = 1), bold = T) %>%
      row_spec(0,bold=TRUE)  %>%
      column_spec(2, background= count.ct(data = x4,2),bold = T) %>%
      column_spec(3, background= count.ct(data = x4,3),bold = T) %>%
      column_spec(4, background= count.ct(data = x4,4),bold = T) %>%
      column_spec(5, background= count.ct(data = x4,5),bold = T) %>% 
      column_spec(6, background= count.ct(data = x4,6),bold = T) %>% 
      column_spec(7, background= count.ct(data = x4,7),bold = T) %>% 
      column_spec(8, background= count.ct(data = x4,8),bold = T) %>% 
      column_spec(9, background= count.ct(data = x4,9),bold = T) %>% 
      column_spec(10, background= count.ct(data = x4,10),bold = T) %>% 
      column_spec(12, background= count.ct(data = x4,12),bold = T)%>% 
      column_spec(13, background= count.ct(data = x4,13),bold = T)%>% 
      column_spec(14, background= count.ct(data = x4,14),bold = T)%>% 
      column_spec(15, background= count.ct(data = x4,15),bold = T)%>% 
      column_spec(16, background= count.ct(data = x4,16),bold = T)%>% 
      column_spec(17, background= count.ct(data = x4,17),bold = T)%>% 
      column_spec(18, background= count.ct(data = x4,18),bold = T)%>%
      column_spec(19, background= count.ct(data = x4,19),bold = T) %>% 
      column_spec(20, background= count.ct(data = x4,20),bold = T) %>%
      # column_spec(21, background= count.ct(data = x4,21),bold = T) %>%
      # column_spec(22, background= count.ct(data = x4,22),bold = T) %>%
      column_spec(11, color = count.score(data = x4, j = 11,color1 = 'white', color2 = 'black'),
                  background = count.score(data = x4, j = 11,color1 = brewer.pal(8,'Dark2')[1], color2 = 'white')) %>%
      column_spec(21, color = count.score(data = x4, j = 21,color1 = 'white', color2 = 'black'),
                  background = count.score(data = x4, j = 21,color1 = brewer.pal(8,'Dark2')[1], color2 = 'white')) %>%
      column_spec(22, color = count.score(data = x4, j = 22,color1 = 'white', color2 = 'black', value = 54),
                  background = count.score(data = x4, j = 22,color1 = brewer.pal(8,'Dark2')[3], color2 = 'white', value = 54),bold = T) %>%
      save_kable(file = paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/SimulationEvaluation/',covariate.type,'_',name1,'_FPs_Null.pdf'))
    
    ## For FDIs
    
    head(res1)
    df = res1 %>%filter(covariate.types==covariate.type) %>% dplyr::select(depth.conf.factors,nOTUs, nSams, methods, value) 
    df[df$value>0,'ct'] = 1
    df[df$value==0,'ct'] = 0
    df = df %>% dplyr::select(-value) 
    formula = paste0('ct~ depth.conf.factors + nOTUs + nSams + methods')
    grid.formula = '. ~ nSams'
    m <- aggregate(as.formula(formula), df, function(x) sum(x[!is.na(x)])) %>% mutate(ct = ct/1000)
    kb = m
    head(kb)
    
    kb[kb$ct > 0.1,'shape'] = 'uncontrol'
    kb[kb$ct <= 0.05,'shape'] = 'fdr_0.05'
    kb[kb$ct <= 0.1 & kb$ct>0.05,'shape'] = 'fdr_0.1'
    kb[kb$ct <= 0.2 & kb$ct>0.1,'shape'] = 'fdr_0.2'
    shapes = c('fdr_0.05' =8, 'uncontrol' = 19, 'fdr_0.1' = 17,'fdr_0.2' = 18)
    
    
    head(kb)
    m = kb
    letters1 = c(letters,rev(toupper(letters)))
    for(i in 1:length(unique(m$methods))){
      m[(m$methods ==unique(m$methods)[i]),'label'] = letters1[i]
    }
    m$legend = paste0(m$label,':',m$methods)
    cols <- c('ZicoSeq'=brewer.pal(9,'Paired')[5], 'Omnibus' = brewer.pal(9,'PuBuGn')[7], 'mbzinb' = brewer.pal(9,'PuBuGn')[8], 
              'GLM(quasipoisson)' = brewer.pal(9,'PuRd')[3],'CLRBC'=brewer.pal(9,'Set1')[8],
              'RAIDA'=brewer.pal(11,'BrBG')[7],'LinDA'=brewer.pal(9,'Set1')[8],
              'eBay(Wilcoxon)' = brewer.pal(8,'RdPu')[4],'eBay(t-test)' = brewer.pal(8,'RdPu')[6], 'ANCOM-BC' = brewer.pal(8,'Set2')[7],
              'LDM'=brewer.pal(9,'Set1')[1], 'DACOMP'=brewer.pal(8,'Dark2')[6],'Aldex2(Wilcoxon)'=brewer.pal(8,'YlGn')[3],'Aldex2(t-test)'=brewer.pal(8,'YlGn')[2],
              'Beta-binomial'=brewer.pal(11,'BrBG')[2],
              'DESeq2'=brewer.pal(9,'Greys')[7], 'Wrench+DESeq2'= brewer.pal(9,'Greys')[3], 'GMPR+DESeq2'=brewer.pal(9,'Greys')[5],
              'TSS+Wilcoxon'=brewer.pal(9,'Purples')[7],'Rarefy+Wilcoxon'= brewer.pal(9,'Set1')[5], 'GMPR+Wilcoxon'= brewer.pal(9,'Purples')[5], 'Wrench+Wilcoxon'= brewer.pal(9,'Purples')[3], 
              'TSS+Spearman'=brewer.pal(9,'Purples')[7],'Rarefy+Spearman'=brewer.pal(9,'Set1')[5], 'GMPR+Spearman'= brewer.pal(9,'Purples')[5], 
              'TSS+t-test'=brewer.pal(9,'GnBu')[7],'Rarefy+t-test'= brewer.pal(9,'GnBu')[5], 'GMPR+t-test'= brewer.pal(9,'GnBu')[5], 'Wrench+t-test'= brewer.pal(9,'GnBu')[3], 
              'edgeR'=brewer.pal(9,'Greens')[7], 'Wrench+edgeR'=brewer.pal(9,'Greens')[3], 'GMPR+edgeR'=brewer.pal(9,'Greens')[5], 
              'metagenomeSeq'=brewer.pal(9,'Blues')[7], 'Wrench+metagenomeSeq'=brewer.pal(9,'Blues')[5])
    
    cols1 = cols
    
    s = unique(m$legend)
    for(i in 1:length(cols1)){
      if(names(cols1)[i] %in% gsub('.*:','',s)){
        names(cols1)[i] = s[which((names(cols1)[i] == gsub('.*:','',s)))]
      }else{
        cat('This legend does not exist! \n')
      }
    }
    # save(m,file = paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/',name1,covariate.type,'BiasDAs.Rdata'))
    head(m)
    
    size = 6
    
    p0 = ggplot(m,aes(x = nOTUs, y = ct, fill = legend)) + # fill = legend shhould be here, or the geom_text label will have same position of label
      geom_bar(position = position_dodge(width = 0.8), stat="identity", width = 0.1) +
      geom_point(aes(color = legend, shape = shape),  size = 5, position = position_dodge(width = 0.8)) +
      geom_text(aes(label=label), position=position_dodge(width=0.8),size=5, vjust = -0.7) + 
      facet_grid(nSams~depth.conf.factors, scales = 'free_y') +
      scale_fill_manual(values = c(rep('grey50',length(df$methods %>% unique())))) +
      scale_color_manual(values = cols1) +
      scale_shape_manual(values = shapes) +
      scale_y_continuous(limits = c(0,1),breaks=c(0.05, 0.1, 0.2, 1)) +
      theme_bw() +thw +
      labs(y = 'FDIs', x = 'OTU number', fill = "", shape = "", color = "",
           title = paste0("Default: no differential taxa (", covariate.type,"-",name1,')')) +
      geom_hline(yintercept = 0.05, colour = '#238823', linetype = 'dashed', size =  0.5) + # green
      geom_hline(yintercept = 0.1, colour = '#ffbf00', linetype = 'dashed', size =  0.5) + # yellow
      geom_hline(yintercept = 0.2, colour = '#e03531', linetype = 'dashed', size =  0.5) + # red
      guides(shape = FALSE, fill = FALSE)
    ggsave(paste0('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/SimulationEvaluation/',name1,'_',covariate.type,'_FDIs_barplot_Null.pdf'), width = 25, height =20)
    
    
    
    
    ## kable FDIs
    kb = m
    kb[kb$ct > 0.1,'shape'] = 'x'
    kb[kb$ct <= 0.05,'shape'] = '***'
    kb[kb$ct <= 0.1 & kb$ct>0.05,'shape'] = '**'
    kb[kb$ct <= 0.2 & kb$ct>0.1,'shape'] = '*'
    
    x1 = kb %>% dplyr::select(depth.conf.factors,nOTUs, nSams, methods, shape) %>% filter(methods !='RioNorm2')
    x1$nOTUs = factor(x1$nOTUs, levels = c('OTU=50','OTU=200','OTU=500'))
    x1$nSams = factor(x1$nSams, levels = c('sample=50','sample=100','sample=200'))
    x1$depth.conf.factors = factor(x1$depth.conf.factors, levels = c('None','Depth confounding'))
    x2 = x1 %>% unite('grp',c('depth.conf.factors','nOTUs','nSams'))
    x3 = x2 %>% spread(c('grp'), shape)
    head(x3)
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
      for(i in 1:nrow(data)){
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
      add_header_above(c("Stool FDIs" = 1, "None" =9," " = 1, "Depth confounding" = 9," " = 1," " = 1), bold = T) %>%
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
                  background = count.score(data = x4, j = 22,color1 = brewer.pal(8,'Dark2')[3], color2 = 'white', value = 54),bold = T) %>%
      save_kable(file = paste0(getwd(),'/result/SimulationEvaluation/',covariate.type,'_',name1,'_FDI_Null.pdf'))
    # as_image(width = 7, height = 3, file = paste0(getwd(),'/plot/',covariate.type,'_',name1,'.pdf')) 
    
  }
  
}
















getwd()
load('sim/allC1_nonWinsor_res.Rdata')
re0 = melt(res)
colnames(re0) <- c('depth.mus', 'diff.otu.modes', 'model', 'nSub', 'nOTUs', 
                    'covariate.types', 'depth.conf.factors', 'diff.otu.directs', 'confounder.types', 'nSams', 
                    'covariate.eff.means', 'diff.otu.pcts','measures', 'methods', 'value','iters')
re0$value <- re0$value/0.05
re2 <- data_summary(re0, formula = 'value ~ diff.otu.modes + covariate.types + covariate.eff.means +diff.otu.pcts +measures+ methods')
name0 = 'covariate.eff.means'
grid.formula = 'diff.otu.modes ~ diff.otu.pcts'
ggplot(re2 %>% filter(measures =='FDR' & covariate.types=='binary' & covariate.eff.means %in% c('L2','L3','L5')), aes(x = covariate.eff.means, y = value,  fill = methods)) +
  theme_bw() +
  geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7) +
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.7, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
  scale_fill_manual(values = cols) +
  # scale_y_continuous(limits = c(0,1),expand = c(0.1, 0, 0, 0)) + 
  facet_grid(diff.otu.pcts~diff.otu.modes, scales = 'free_y')+
  labs(y = 'obsevered FDR/expected FDR(0.05)', x = stri_trans_totitle(name), color = "Methods") +
  thw +
  geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5)



