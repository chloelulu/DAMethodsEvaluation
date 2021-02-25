setwd("~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/realdata/result_preprocessing/")
files = list.files(getwd(), pattern = 'Rdata$')
library(tidyverse);library(reshape2);library(RColorBrewer);library(scales);library(patchwork)

## calculate # of findings in each method
FDR = list()
for(file in files){
  load(file)
  cat(file,'\n')
  try({
    fdr <- res.sum  %>% column_to_rownames('taxa') %>% dplyr::select(contains('fdr.')) #%>% replace(is.na(.), 1)
    fdr[,names(which(sapply(fdr, class)!='numeric'))] <- as.numeric(as.character(fdr[,names(which(sapply(fdr, class)!='numeric'))]))
    fdr[is.na(fdr)] <- 1
    fdr <- colSums(fdr <= 0.05, na.rm = T) %>% as.data.frame()
    colnames(fdr) = file
    fdr <- fdr %>% rownames_to_column('methods')
  }
  )
  FDR[[file]] <- fdr
}
  

res.sum0 <- Reduce(
  function(x, y, ...) {
    if (is.null(x)){
      y
    } else if (is.null(y)){
      x
    }
    else{
      merge(x, y, all = TRUE, ...)
    }},
  FDR
)
res.sum0 = melt(res.sum0)



setwd("/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/realdata/result_preprocessing/v1")
files = list.files(getwd(), pattern = 'Rdata$')
files  = files[!(files %in% c("FengQ_2015_res.Rdata","HanniganGD_2017_res.Rdata","ThomasAM_2018a_res.Rdata","ThomasAM_2018b_res.Rdata","VogtmannE_2016_res.Rdata","YuJ_2015_res.Rdata","ZellerG_2014_res.Rdata",
                            "KarlssonFH_2013_res.Rdata","QinJ_2012_res.Rdata"))]
FDR = list()
for(file in files){
  load(file)
  cat(file,'\n')
  try({
    fdr <- res.sum  %>% column_to_rownames('taxa') %>% dplyr::select(contains('fdr.')) #%>% replace(is.na(.), 1)
    fdr[,names(which(sapply(fdr, class)!='numeric'))] <- as.numeric(as.character(fdr[,names(which(sapply(fdr, class)!='numeric'))]))
    fdr[is.na(fdr)] <- 1
    fdr <- colSums(fdr <= 0.05, na.rm = T) %>% as.data.frame()
    colnames(fdr) = gsub('.Rdata','',file)
    fdr <- fdr %>% rownames_to_column('methods')
    
    
  }
  )
  FDR[[file]] <- fdr
}


res.sum1 <- reduce_list(FDR)
head(res.sum1)
res.sum1 = melt(res.sum1)


FDR1 = rbind(res.sum0, res.sum1)
FDR1$methods = gsub('fdr.','',FDR1$methods)
FDR1 = FDR1 %>% filter(!(methods %in% c('CLRBC','DESeq2','DESeq2.Wrench', 'DESeq2.gmpr','edgeR', 'edgeR.gmpr','edgeR.Wrench')))

FDR1$methods = as.character(FDR1$methods)
FDR1$methods[FDR1$methods =='Wilcox'] = 'TSS+Wilcoxon'
FDR1$methods[FDR1$methods =='Rarefy'] = 'Rarefy+Wilcoxon'
FDR1$methods = gsub('ANCOMBC','ANCOM-BC',FDR1$methods)
FDR1$methods = gsub('glmquassi','GLM(quasipoisson)',FDR1$methods)
FDR1$methods = gsub('DESeq2.gmpr','GMPR+DESeq2',FDR1$methods)
FDR1$methods = gsub('DESeq2.Wrench','Wrench+DESeq2',FDR1$methods)
FDR1$methods = gsub('edgeR.gmpr','GMPR+edgeR',FDR1$methods)
FDR1$methods = gsub('edgeR.Wrench','Wrench+edgeR',FDR1$methods)
FDR1$methods = gsub('MSeq2.Wrench','Wrench+metagenomeSeq',FDR1$methods)
FDR1$methods = gsub('MSeq2','metagenomeSeq',FDR1$methods)
FDR1$methods = gsub('eBayW','eBay(Wilcoxon)',FDR1$methods)
FDR1$methods = gsub('eBayt','eBay(t-test)',FDR1$methods)
FDR1$methods = gsub('BBinomial','Beta-binomial',FDR1$methods)
FDR1$methods = gsub('Aldex2we','Aldex2(t-test)',FDR1$methods)
FDR1$methods = gsub('^Aldex2$','Aldex2(Wilcoxon)',FDR1$methods)
FDR1$methods = gsub('^Rarefyttest$','Rarefy+t-test',FDR1$methods)
FDR1$methods = gsub('^ttest$','TSS+t-test',FDR1$methods)
FDR1$methods = gsub('^Omnibus$','mbzinb',FDR1$methods)

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

head(FDR1)
mn = aggregate(value~methods, FDR1 , function(x) median(x[!is.na(x)]))
mn = mn[order(mn$value),]
FDR1 <- within(FDR1, methods <- factor(methods, levels=rev(mn$methods)))
head(FDR1)
ggplot(FDR1 %>% na.omit(),aes(x = methods, y = value)) +
  stat_boxplot(geom = "errorbar", width = 0.5,lwd=0.2) +
  # geom_boxjitter(aes(fill= methods),outlier.color = NA, jitter.shape = 21, 
  #                jitter.color = NA, jitter.alpha = 0.5,
  #                jitter.width = 0.2, errorbar.draw =T) +
  geom_boxplot(aes(fill = methods),outlier.size =1, outlier.colour = 'grey', width = 0.9) +
  # geom_jitter(aes(color = methods),position=position_jitter(0.3),size = 0.7, alpha  = 0.7) +
  theme_bw() +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  scale_y_continuous(trans = sqrt_trans(),
                     breaks = trans_breaks("sqrt", function(x) x^2),
                     labels = trans_format("sqrt", math_format(.x^2))) +
  theme(axis.text.x = element_text(color="black", size =22, angle = 90, vjust = 0.4, hjust = 0.95),
        axis.text.y = element_text(color="black", size = 22),
        axis.title = element_text(color="black", size = 22),
        strip.text = element_text(size = 22),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        legend.title=element_text(size=22),
        legend.text = element_text(size=22)) +
  ylab('Findings') + xlab('') + guides(fill = F, color = F)
# ggsave(file = '/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/Findings_preprocessed.pdf', width = 10, height = 8, dpi = 100)




## found by 80% of the methods
setwd("/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/realdata/result_preprocessing/")
files = list.files(getwd(), pattern = 'Rdata$')
files  = files[!(files %in% c("FengQ_2015_res.Rdata","HanniganGD_2017_res.Rdata","ThomasAM_2018a_res.Rdata","ThomasAM_2018b_res.Rdata","VogtmannE_2016_res.Rdata","YuJ_2015_res.Rdata","ZellerG_2014_res.Rdata",
                              "KarlssonFH_2013_res.Rdata","QinJ_2012_res.Rdata"))]

pct = 0.2

FDRs <- TPRs <- list()
for(file in files){
  try({
    load(file)
    fdr.df <- res.sum  %>% column_to_rownames('taxa') %>% dplyr::select(contains('fdr.')) #%>% replace(is.na(.), 1)
    fdr.df[,names(which(sapply(fdr.df, class)!='numeric'))] <- as.numeric(as.character(fdr.df[,names(which(sapply(fdr.df, class)!='numeric'))]))
    
    load(paste0('v1/',file))
    fdr.df1 <- res.sum  %>% column_to_rownames('taxa') %>% dplyr::select(contains('fdr.')) #%>% replace(is.na(.), 1)
    fdr.df1[,names(which(sapply(fdr.df1, class)!='numeric'))] <- as.numeric(as.character(fdr.df1[,names(which(sapply(fdr.df1, class)!='numeric'))]))
  
    fdr.df2 <- full_join(fdr.df %>% rownames_to_column('id'), fdr.df1 %>% rownames_to_column('id')) %>% column_to_rownames('id')
    
    fdr.df2[is.na(fdr.df2)] <- 1
    colnames(fdr.df2) = gsub('fdr.','',colnames(fdr.df2))
    fdr.df2 <- fdr.df2 %>% dplyr::select(-colnames(fdr.df2)[grep('DESeq2|CLRBC|edgeR',colnames(fdr.df2))])
  })
  
  fdr.df2$truth <- FALSE
  for(i in 1:nrow(fdr.df2)){
    if((sum(fdr.df2[i,] <= 0.05)/ncol(fdr.df2)) > pct){
      fdr.df2[i,'truth'] = TRUE
    }
  }
  fdrs <- tprs <- NULL
  for(j in 1:(ncol(fdr.df2)-1)){
    tp <- sum(fdr.df2[,j] <= 0.05 & fdr.df2$truth==TRUE, na.rm = T)
    tn <- sum(fdr.df2[,j] > 0.05 & fdr.df2$truth==FALSE, na.rm = T)#sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(fdr.df2[,j] > 0.05 &  fdr.df2$truth ==TRUE, na.rm = T)
    fp <- sum(fdr.df2[,j] <= 0.05 &  fdr.df2$truth ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    fdrs <- c(fdrs, fdr)
    tprs <- c(tprs, tpr)
  }
  names(fdrs) = names(tprs) = colnames(fdr.df2)[!(colnames(fdr.df2) %in% 'truth')]
  FDRs[[file]] <- fdrs
  TPRs[[file]] <- tprs
}

# save(FDRs,file = paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/realdata/FDRs_realdata_',pct,'.Rdata'))
# save(TPRs,file = paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/realdata/TPRs_realdata_',pct,'.Rdata'))

pct = 0.2
load(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/realdata/FDRs_realdata_',pct,'.Rdata'))
load(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/realdata/TPRs_realdata_',pct,'.Rdata'))
func <- function(FDRs){
  FDR1 <- as.data.frame(do.call(rbind, FDRs)) %>% melt(.)
  colnames(FDR1)[1]= 'methods'
  FDR1$methods = as.character(FDR1$methods)
  FDR1$methods[FDR1$methods =='Wilcox'] = 'TSS+Wilcoxon'
  FDR1$methods[FDR1$methods =='Rarefy'] = 'Rarefy+Wilcoxon'
  FDR1$methods = gsub('ANCOMBC','ANCOM-BC',FDR1$methods)
  FDR1$methods = gsub('glmquassi','GLM(quasipoisson)',FDR1$methods)
  FDR1$methods = gsub('DESeq2.gmpr','GMPR+DESeq2',FDR1$methods)
  FDR1$methods = gsub('DESeq2.Wrench','Wrench+DESeq2',FDR1$methods)
  FDR1$methods = gsub('edgeR.gmpr','GMPR+edgeR',FDR1$methods)
  FDR1$methods = gsub('edgeR.Wrench','Wrench+edgeR',FDR1$methods)
  FDR1$methods = gsub('MSeq2.Wrench','Wrench+metagenomeSeq',FDR1$methods)
  FDR1$methods = gsub('MSeq2','metagenomeSeq',FDR1$methods)
  FDR1$methods = gsub('eBayW','eBay(Wilcoxon)',FDR1$methods)
  FDR1$methods = gsub('eBayt','eBay(t-test)',FDR1$methods)
  FDR1$methods = gsub('BBinomial','Beta-binomial',FDR1$methods)
  FDR1$methods = gsub('Aldex2we','Aldex2(t-test)',FDR1$methods)
  FDR1$methods = gsub('^Aldex2$','Aldex2(Wilcoxon)',FDR1$methods)
  FDR1$methods = gsub('^Rarefyttest$','Rarefy+t-test',FDR1$methods)
  FDR1$methods = gsub('^ttest$','TSS+t-test',FDR1$methods)
  FDR1$methods = gsub('^Omnibus$','mbzinb',FDR1$methods)
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
  
  mn = aggregate(value~methods, FDR1 , function(x) median(x[!is.na(x)]))
  mn = mn[order(mn$value),]
  FDR1 <- within(FDR1, methods <- factor(methods, levels=mn$methods))
  delete = c('DESeq2','Wrench+DESeq2','edgeR','Wrench+edgeR','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon','RioNorm2','LINDA2','LinDA')
  
  FDR1 <- FDR1 %>% filter(!(methods %in% delete))
}
data1 <- func(FDRs = FDRs)
data2 <- func(FDRs = TPRs)

p1 = ggplot(data1,aes(x = methods, y = value)) +
  stat_boxplot(geom = "errorbar", width = 0.5,lwd=0.2) +
  geom_boxplot(aes(fill = methods),outlier.size = 0.7, outlier.colour = 'grey', width = 0.9) +
  theme_bw() +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(color="black", size =22, angle = 90, vjust = 0.2, hjust = 1),
        axis.text.y = element_text(color="black", size = 22),
        axis.title = element_text(color="black", size = 22),
        strip.text = element_text(size = 22),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        legend.title=element_text(size=22),
        legend.text = element_text(size=22)) +
  ylab('FDR') + xlab('') + guides(fill = F) +
  geom_hline(yintercept = 0.05, colour = brewer.pal(8,'Set1')[3], linetype = 'dashed', size = 0.5) +
  geom_hline(yintercept = 0.1, colour = brewer.pal(8,'Set2')[6], linetype = 'dashed', size = 0.5)  +
  geom_hline(yintercept = 0.2, colour = brewer.pal(8,'Set1')[1], linetype = 'dashed', size = 0.5) 


head(data2)
(aggregate(value ~ methods, data2, function(x) median(x)) %>% filter(value <=0.7 &value >0.5))[,1]
(aggregate(value ~ methods, data2, function(x) mean(x)) %>% filter(value >0.6))[,1]

p2 = ggplot(data2,aes(x = methods, y = value)) +
  stat_boxplot(geom = "errorbar", width = 0.5,lwd=0.2) +
  geom_boxplot(aes(fill = methods),outlier.size =0.7, outlier.colour = 'grey', width = 0.9) +
  theme_bw() +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(color="black", size =22, angle = 90, vjust = 0.2, hjust = 1),
        axis.text.y = element_text(color="black", size = 22),
        axis.title = element_text(color="black", size = 22),
        strip.text = element_text(size = 22),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        legend.title=element_text(size=22),
        legend.text = element_text(size=22)) +
  ylab('Power') + xlab('') + guides(fill = F) 
p1 +p2
# ggsave(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/plot/FDRs_realdata_',pct,'.pdf'), width = 12, height = 8, dpi = 100)
ggsave(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/FDRs_realdata_',pct,'.pdf'), width = 12, height = 8, dpi = 100)
