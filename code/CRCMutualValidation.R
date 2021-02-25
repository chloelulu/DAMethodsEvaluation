setwd('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/')
source("code/wrappers.R")
source('code/zeroinfl.plus.daa.R')
source('code/zeroinfl.plus.github.R')
source('code/raida.R') # The default installed RAIDA package has bugs, use this author suggested version
pkg <- c('ANCOMBC','mbzinb','ZicoSeq','phyloseq','microbiome','aod','reshape','MASS','GMPR','corncob','readr','DESeq2', 'ALDEx2', 'metagenomeSeq', 'edgeR', 'GUniFrac', 'grDevices', 'dirmult', 'exactRankTests','nlme', 'dplyr', 
         'magrittr', 'tidyr', 'protoclust', 'ggplot2', 'compositions','rmutil','tibble','reticulate','dacomp','LDM','Wrench','RioNorm2','ggiraphExtra')
lapply(pkg, require, character.only = TRUE)
rdata <- list.files('data/CRCdata')
calculate <- function(dat, method){
  wrapper <- match.fun(methods_funs[[method]])
  cat(paste(method, '\n'))
  out <- wrapper(dat, FDR=F)
  sig <- out$sig
  if(is.na(sig)){
    sig <- 0
  }
  names(sig) <- method
  res <- as.data.frame(out$res) %>% dplyr::select('otu.id','pval','fdr')
  colnames(res) <- c('taxa',paste0('pval.',method),paste0('fdr.',method))
  # write.csv(res,paste0(name,'_',paste0(method,'_noFil.csv')), row.names = F)
  return(list(res = res, sig = sig))
}


files = list.files('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/data/CRCdata/',pattern = 'Rdata$')
setwd('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/data/CRCdata/')
mn = NULL
for(file in files){
  load(file)
  cat(mean(sample_sums(phy)),'\n')
  mn = c(mn,mean(sample_sums(phy)))
  phy1 = prune_samples(sample_sums(phy)>(quantile(sample_sums(phy), probs = seq(0,1,0.25))[2] %>% as.numeric() %>% round()),phy)
  cat(file,':',quantile(sample_sums(phy), probs = seq(0,1,0.25))[2] %>% as.numeric() %>% round(),';',table(phy1@sam_data$grp),';',mean(sample_sums(phy)),'\n')
  }
min(mn)


# dir = '/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/CRC/'

files <- c("FengQ_2015.Rdata","HanniganGD_2017.Rdata","ThomasAM_2018a.Rdata","ThomasAM_2018b.Rdata","VogtmannE_2016.Rdata","YuJ_2015.Rdata","ZellerG_2014.Rdata") # CRC datasets
files <- c("KarlssonFH_2013.Rdata","QinJ_2012.Rdata") # T2D datasets
# res.df = list()
# for(file in files){
#   load(paste0(dir,file))
#   res = NULL
#   for(i in 1:length(res.sum)){
#     res = cbind(res, res.sum[[i]])
#   }
#   res = as.data.frame(res)
#   res = res[grep('^s_',res$V1),]
#   colnames(res) = names(res.sum)
#   res.df[[file]] = res
# }
# save(res.df, file='result/CRC_DA_res_0203.Rdata')


title = gsub('.Rdata','',files)
res.sum <- sig.sum <- list()
prev =0;minp = 0
for(i in 1:length(title)){
  name = title[i]
  load(files[i])
  # table(phy@sam_data$study_condition)
  # load(paste0('result/',files[i]))
  # load(paste0('curated/','TettAJ_2016.Rdata'))
  # phy = prune_samples(sample_sums(phy)>30000,phy)
  dep <- sample_sums(phy)
  diff.seq.p <- summary(aov(dep ~ phy@sam_data$grp))[[1]][1, 'Pr(>F)']
  # sort(dep)
  if (!is.na(diff.seq.p) & diff.seq.p <= 0.05) {
    cat("Signficant sequencing depth confounding!\n")
    cat("For parametric test with sequence depth adjustment, please be cautious about the results!\n")
    cat("There may be potential residual sequence depth confounding! Let's do rarefraction!\n")
    if(mean(sample_sums(phy))>1000000){
      # for the curated datasets, metagenomic data 
      if(min(sample_sums(phy)) < 30000){
        phy = rarefy_even_depth(phy, sample.size = 30000)
      }else{
        phy = rarefy_even_depth(phy)
      }
    }else{
      phy = rarefy_even_depth(phy, sample.size = round(as.numeric(quantile(sample_sums(phy), probs = seq(0,1,0.25))[2])))
    }
    # save(phy, file = paste0('result/',name,'_Rarefied.RData'))
  }else{
    cat('No depth confounding was found!\n')
  }
  otu.tab.sim <- otu_table(phy)@.Data
  meta.dat <- as.data.frame(as.matrix(sample_data(phy)))[,'grp',drop = F]
  covariates <- meta.dat$grp# %>% sample() # sample for shuffle purpose
  names(covariates) <- rownames(meta.dat)
  cat('1.Loading data finshed! \n')
  ##-- Preprocessing
  ##--- Normalization
  ## GMPR
  size.factor <- GMPR(otu.tab.sim)
  gmpr.size.factor <- size.factor[!is.na(size.factor)]
  gmpr.comm <- otu.tab.sim[,names(gmpr.size.factor)]
  gmpr.counts <- t(t(gmpr.comm) / gmpr.size.factor)
  cat('2.GMPR finshed! \n')
  ## gmpr otu table prevelance filtration
  prop.gmpr <- t(t(gmpr.counts) / colSums(gmpr.counts))
  prop.gmpr <- prop.gmpr[rowSums(prop.gmpr!=0) > prev * ncol(prop.gmpr), , drop=FALSE]
  gmpr.counts <- gmpr.counts[rownames(prop.gmpr), , drop=FALSE]
  
  ## gmpr otu table minimal abundance filtration
  prop.gmpr <- prop.gmpr[rowMaxs(prop.gmpr) > minp, , drop=FALSE]
  gmpr.counts <- gmpr.counts[rownames(prop.gmpr), , drop=FALSE] # used for analysis
  
  gmpr.covariates <- covariates[!is.na(size.factor)]
  gmpr.meta.dat <- meta.dat[rownames(meta.dat) %in% colnames(gmpr.counts),, drop=FALSE]
  
  ## Wrench:Does not support continuous covarites !!! needs to be cleaned up and recheck!
  W <- try(wrench(otu.tab.sim, condition=covariates))
  if(inherits(W, "try-error")){
    compositionalFactors = normalizationFactors = Wrench.nf = Wrench.covariates = Wrench.truth =Wrench.meta.dat = NULL
    cat('3.Wrench ERROR! \n')
  }else{
    compositionalFactors <- W$ccf[!is.na(W$ccf)]
    normalizationFactors <- W$nf[!is.na(W$nf)]
    Wrench.nf <- t(t(otu.tab.sim[,names(normalizationFactors)]) / normalizationFactors)
    ##  otu table prevelance filtration
    prop.Wrench <- t(t(Wrench.nf) / colSums(Wrench.nf))
    prop.Wrench <- prop.Wrench[rowSums(prop.Wrench!=0) > prev * ncol(prop.Wrench), , drop=FALSE]
    Wrench.nf <- Wrench.nf[rownames(prop.Wrench), , drop=FALSE]
    
    ## gmpr otu table minimal abundance filtration
    prop.Wrench <- prop.Wrench[rowMaxs(prop.Wrench) > minp, , drop=FALSE]
    Wrench.nf <- Wrench.nf[rownames(prop.Wrench), , drop=FALSE] # used for analysis
    
    compositionalFactors <- compositionalFactors[colnames(Wrench.nf)]
    Wrench.covariates <- covariates[which(!is.na(compositionalFactors)) | which(!is.na(normalizationFactors))]
    Wrench.meta.dat <- meta.dat[rownames(meta.dat) %in% colnames(Wrench.nf),, drop=FALSE]
    cat('3.Wrench finshed! \n')
  }
  
  ## raw otu table prevelance filtration
  prop <- t(t(otu.tab.sim) / colSums(otu.tab.sim))
  prop <- prop[rowSums(prop!=0) > prev * ncol(prop), , drop=FALSE]
  otu.tab.sim <- otu.tab.sim[rownames(prop), , drop=FALSE]
  ## raw otu tableminimal abundance filtration
  prop <- prop[rowMaxs(prop) > minp, , drop=FALSE]
  otu.tab.sim <- otu.tab.sim[rownames(prop), , drop=FALSE] # used for analysis
  cat('4.Preprocessing finshed! \n')
  
  dat <- list(counts = otu.tab.sim, prop = prop, covariates = covariates, meta.dat = meta.dat,
              Wrench.nf = Wrench.nf, compositionalFactors= compositionalFactors, normalizationFactors = normalizationFactors, Wrench.covariates = Wrench.covariates,Wrench.meta.dat = Wrench.meta.dat,
              gmpr.comm=gmpr.comm, gmpr.counts = gmpr.counts, gmpr.size.factor = gmpr.size.factor, gmpr.covariates=gmpr.covariates, gmpr.meta.dat = gmpr.meta.dat)

  methods_funs <- list('ZicoSeq'= "ZicoSeq.wrapper", 'Rarefy'= "Rarefy.wrapper", 'Aldex2'= "Aldex2.wrapper", 'mbzinb'="mbzinb.wrapper",
                       'RAIDA'="RAIDA.wrapper", 'DACOMP' = "dacomp.wrapper", 'LDM'="ldm.wrapper", 'RioNorm2'="RioNorm.wrapper",'glmquassi'='glmquassi.wrapper',
                       'BBinomial' = 'BBinomial.wrapper','ANCOMBC'="ANCOMBC.wrapper", 'ttest' ='ttest.wrapper','Rarefyttest' = 'Rarefyttest.wrapper',
                       'DESeq2'="DESeq.wrapper", 'DESeq2.Wrench'="DESeq.Wrench.wrapper", 'DESeq2.gmpr'="DESeq.gmpr.wrapper",
                       'Wilcox'='wilcox.wrapper', 'edgeR'="edgeR.wrapper", 'edgeR.gmpr'="edgeR.gmpr.wrapper", 'edgeR.Wrench'="edgeR.Wrench.wrapper",
                       'MSeq2.Wrench'="metagenomeSeq2.Wrench.wrapper", 'MSeq2'="metagenomeSeq2.wrapper",'CLRBC'="CLRBC.wrapper",'eBayt'="eBayt.wrapper",'eBayW'="eBayW.wrapper",'Aldex2we'="Aldex2we.wrapper")
  methods <- c(
    # 'eBayt','eBayW','Aldex2we','Rarefy','Aldex2','RAIDA','DACOMP','LDM','glmquassi','BBinomial','ANCOMBC','Wilcox','MSeq2.Wrench', 'MSeq2',
               'mbzinb', 'CLRBC','ttest','Rarefyttest')
  
  res <- sig <- NULL
  for(method in methods){
    tryCatch({
    cal <- calculate(dat,method = method)
    res[[method]] <- cal$res
    sig[[method]] <- cal$sig
    },error =function(e){cat(method,' ERROR : ',conditionMessage(e), "\n")})
  }
      
  multi_full <- Reduce(
    function(x, y, ...) {
      if (is.null(x)){
        y
      } else if (is.null(y)){
        x
      }
      else{
        merge(x, y, all = TRUE, ...)
      }},
    res
  )
  sig.sum[[name]] <- sig
  res.sum[[name]] <- multi_full
}


# save(sig.sum, file='/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/CRC_DA_sig1.Rdata')
# save(res.sum, file='/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/CRC_DA_res1.Rdata')

# load('result/CRC_DA_res.Rdata')


# Below are the biological markers
taxa.T2D = c("s__Streptococcus_mutans","s__Streptococcus_vestibularis","s__Lactobacillus_gasseri","s__Lactococcus_raffinolactis","s__Clostridium_symbiosum",
             "s__Clostridium_citroniae", 's__Clostridium_asparagiforme',"s__Clostridium_bolteae","s__Clostridium_bartlettii","s__Clostridium_clostridioforme",
             "s__Clostridium_methylpentosum","s__Lactobacillus_vaginalis","s__Lactobacillus_oris","s__Actinomyces_turicensis","s__Actinomyces_viscosus") # T2D
taxa.CRC = c("s__Fusobacterium_nucleatum","s__Fusobacterium_necrophorum","s__Fusobacterium_mortiferum","s__Fusobacterium_periodonticum","s__Peptostreptococcus_stomatis",
             "s__Peptostreptococcus_anaerobius", "s__Prevotella_bivia", "s__Prevotella_intermedia","s__Prevotella_stercorea","s__Gemella_haemolysans", 
             "s__Gemella_morbillorum","s__Streptococcus_gallolyticus", "s__Streptococcus_constellatus", "s__Streptococcus_oligofermentans", "s__Streptococcus_peroris") # CRC
setwd("/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/")
CRC_T2D = NULL
for(name in c('CRC','T2D')){
  load(paste0('result/',name,'_DA_res.Rdata'))
  res = NULL
  for(i in 1:length(res.sum)){
    cat(names(res.sum)[i],'\n')
    multi_full1 = res.sum[[i]] %>% dplyr::select(c('taxa',colnames(.)[grep('fdr',colnames(.))])) %>% column_to_rownames('taxa')
    if(name =='CRC'){taxa = taxa.CRC}else{taxa=taxa.T2D}
    multi_full3 = multi_full1[taxa,]
    multi_full3 = multi_full3[grep('^s__',rownames(multi_full3)),]
    multi_full3[is.na(multi_full3)] = 1
    res = cbind(res,colSums(multi_full3 <= 0.1))
  }
  colnames(res) = names(res.sum)
  CRC_T2D = cbind(CRC_T2D, res)
 
}

rownames(CRC_T2D) = gsub('fdr.','', rownames(CRC_T2D))
CRC_T2D = as.data.frame(CRC_T2D) %>% rownames_to_column('methods')

CRC_T2D$methods[CRC_T2D$methods =='Wilcox'] = 'TSS+Wilcoxon'
CRC_T2D$methods[CRC_T2D$methods =='Rarefy'] = 'Rarefy+Wilcoxon'
CRC_T2D$methods = gsub('ANCOMBC','ANCOM-BC',CRC_T2D$methods)
CRC_T2D$methods = gsub('glmquassi','GLM(quasipoisson)',CRC_T2D$methods)
CRC_T2D$methods = gsub('DESeq2.gmpr','GMPR+DESeq2',CRC_T2D$methods)
CRC_T2D$methods = gsub('DESeq2.Wrench','Wrench+DESeq2',CRC_T2D$methods)
CRC_T2D$methods = gsub('edgeR.gmpr','GMPR+edgeR',CRC_T2D$methods)
CRC_T2D$methods = gsub('edgeR.Wrench','Wrench+edgeR',CRC_T2D$methods)
CRC_T2D$methods = gsub('MSeq2.Wrench','Wrench+metagenomeSeq',CRC_T2D$methods)
CRC_T2D$methods = gsub('MSeq2','metagenomeSeq',CRC_T2D$methods)
CRC_T2D$methods = gsub('eBayW','eBay(Wilcoxon)',CRC_T2D$methods)
CRC_T2D$methods = gsub('eBayt','eBay(t-test)',CRC_T2D$methods)
CRC_T2D$methods = gsub('BBinomial','Beta-binomial',CRC_T2D$methods)
CRC_T2D$methods = gsub('Aldex2we','Aldex2(t-test)',CRC_T2D$methods)
CRC_T2D$methods = gsub('^Aldex2$','Aldex2(Wilcoxon)',CRC_T2D$methods)



## load mbzinb and other methods
CRC_T2D2 = NULL
for(name in c('CRC','T2D')){
  load(paste0('result/',name,'_DA_res1.Rdata'))
  res = NULL
  for(i in 1:length(res.sum)){
    cat(names(res.sum)[i],'\n')
    multi_full1 = res.sum[[i]] %>% dplyr::select(c('taxa',colnames(.)[grep('fdr',colnames(.))])) %>% column_to_rownames('taxa')
    if(name =='CRC'){taxa = taxa.CRC}else{taxa=taxa.T2D}
    multi_full3 = multi_full1[taxa,]
    multi_full3 = multi_full3[grep('^s__',rownames(multi_full3)),]
    multi_full3[is.na(multi_full3)] = 1
    res = cbind(res,colSums(multi_full3 <= 0.1))
  }
  colnames(res) = names(res.sum)
  CRC_T2D2 = cbind(CRC_T2D2, res)
}
rownames(CRC_T2D2) = gsub('fdr.','', rownames(CRC_T2D2))
CRC_T2D2 = as.data.frame(CRC_T2D2) %>% rownames_to_column('methods')
CRC_T2D2$methods = gsub('^ttest$','TSS+t-test',CRC_T2D2$methods)
CRC_T2D2$methods = gsub('^Rarefyttest$','Rarefy+t-test',CRC_T2D2$methods)


# combine new with old
CRC_T2D = rbind(CRC_T2D, CRC_T2D2) %>% filter(!(methods %in% c('CLRBC','Omnibus')))

CRC_T2D = CRC_T2D %>% column_to_rownames('methods')


CRC = CRC_T2D[,c(1:7),drop = F] # max = 6
T2D = CRC_T2D[,c(8:9),drop = F] # max = 5






x = rownames(CRC_T2D)[apply(CRC_T2D,2,which.max)] %>% table() %>% sort()
barplot(x)
rownames(CRC)[apply(CRC,2,which.max)] %>% table()
rownames(T2D)[apply(T2D,2,which.max)] %>% table()

library(ggradar)

d = melt(CRC_T2D %>% rownames_to_column('methods'))
d$variable = as.factor(d$variable)
d <- within(d, variable <- factor(variable, levels=c(names(sort(colSums(CRC))),names(sort(colSums(T2D))))))
unique(d$methods)
names(sort(-colsSums(CRC_T2D)))
d <- within(d, methods <- factor(methods, levels=names(sort(-rowSums(CRC_T2D)))))
# d$variable <- gsub('_','.',d$variable)
display.brewer.all()
col = c("HanniganGD_2017" = brewer.pal(8,'YlGn')[7],
        "ThomasAM_2018a"= brewer.pal(8,'YlGn')[8],
        "VogtmannE_2016"= brewer.pal(8,'YlGn')[6],
        "ThomasAM_2018b"= brewer.pal(8,'YlGn')[5],
        "ZellerG_2014"= brewer.pal(8,'YlGn')[4],
        "FengQ_2015" = brewer.pal(8,'YlGn')[3],
        "YuJ_2015" = brewer.pal(8,'YlGn')[2],
        "KarlssonFH_2013"= brewer.pal(8,'OrRd')[2], # T2D
        "QinJ_2012"= brewer.pal(8,'OrRd')[1]) # T2D
ggBar(d,aes(x=methods,fill=variable,y=value),stat="identity",polar=TRUE,width=1, addlabel =F,labelsize = 3,
      color="black",size=0.1,reverse=F,interactive=F)  +theme_bw() +
  scale_fill_manual(values = col) +
  theme(axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(color = 'black', size = 20, vjust =0.2, hjust = 0.95, angle = 300),
        axis.ticks=element_blank(),plot.margin = unit(c(2,2,2,2), "cm"),
        legend.title=element_text(size=20),
        legend.text = element_text(size=20),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) + 
  labs(fill = 'Datasets')
length(unique(d$methods))
# ggsave('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/CRC_T2D_ggBar.pdf', width = 15, height = 15, dpi = 100)
ggsave('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/CRC_T2D_ggBar.pdf', width = 15, height = 15, dpi = 100)



## Making upset plot 
pkg <- c('ComplexHeatmap','ZicoSeq','phyloseq','aod','reshape','MASS','GMPR','corncob','ggforce','readr','DESeq2', 'ALDEx2', 'metagenomeSeq', 'edgeR', 'GUniFrac', 'grDevices', 'dirmult', 'exactRankTests','nlme', 'dplyr', 'magrittr', 'tidyr', 'protoclust', 'ggplot2', 'compositions','rmutil','tibble','reticulate','dacomp','LDM','Wrench','RioNorm2')# ,'RAIDA'
suppressPackageStartupMessages(lapply(pkg, require, character.only = TRUE))
setwd("~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/CRC")
load('CRC_DA_res.Rdata')
load('../../result/CRC_DA_res.Rdata')
## apply Upsetplot to the CRC dataset
# each dataset among 21 methods
names = names(res.sum)
Up_plot <- function(res.sum, name, fdr = 0.05){
  Up <- res.sum[[name]] %>% as.data.frame()
  Up$fdr.Aldex2 <- as.numeric(as.character(Up$fdr.Aldex2))
  Up$pval.Aldex2 <- as.numeric(as.character(Up$pval.Aldex2))
  fdr = names(res.sum[[name]])[grep('^fdr',names(res.sum[[name]]))][!(names(res.sum[[name]])[grep('^fdr',names(res.sum[[name]]))] %in% c("fdr.ZicoSeq"))]
  Up = Up %>% dplyr::select(c('taxa',fdr))
  colnames(Up) = gsub('fdr.','',colnames(Up))
  Up[is.na(Up)] = 1
  
  taxa.T2D = c("s__Streptococcus_mutans","s__Streptococcus_vestibularis","s__Lactobacillus_gasseri","s__Lactococcus_raffinolactis","s__Clostridium_symbiosum",
               "s__Clostridium_citroniae", 's__Clostridium_asparagiforme',"s__Clostridium_bolteae","s__Clostridium_bartlettii","s__Clostridium_clostridioforme",
               "s__Clostridium_methylpentosum","s__Lactobacillus_vaginalis","s__Lactobacillus_oris","s__Actinomyces_turicensis","s__Actinomyces_viscosus") # T2D
  taxa.CRC = c("s__Fusobacterium_nucleatum","s__Fusobacterium_necrophorum","s__Fusobacterium_mortiferum","s__Fusobacterium_periodonticum","s__Peptostreptococcus_stomatis",
               "s__Peptostreptococcus_anaerobius", "s__Prevotella_bivia", "s__Prevotella_intermedia","s__Prevotella_stercorea","s__Gemella_haemolysans", 
               "s__Gemella_morbillorum","s__Streptococcus_gallolyticus", "s__Streptococcus_constellatus", "s__Streptococcus_oligofermentans", "s__Streptococcus_peroris") # CRC
  taxa.s = taxa.CRC
  
  Up1 = Up %>% filter(taxa %in% taxa.s)
  
  if(nrow(Up1) >1){
    str(Up)
    Up1 = Up1[rowSums(Up1[ ,-1]<=0.05)>0,] %>% droplevels() %>% na.omit()
    colnames(Up1) = gsub('metagenomeSeq2','MSeq2',colnames(Up1))
    ## For methods as left rowname in plot
    Up = list()
    for (j in 2:ncol(Up1)){
      idx2 = (Up1[,j] <= 0.05)
      Up[[colnames(Up1)[j]]] = Up1[,1][idx2]
    }
    
    Up.lt = make_comb_mat(Up)
  }else{
    Up.lt = NULL
    cat('all methods can not detect the biological taxa! \n')
  }
  return(list(Up.lt = Up.lt, name = name))
}

plt = Up_plot(res.sum, name = names[1])
# pdf(file = paste0(plt$name,'UpSetPlot_MarkerOnly.pdf'))
UpSet(plt$Up.lt)
# dev.off() 

plt = Up_plot(res.sum, name = names[2])
# pdf(file = paste0(plt$name,'UpSetPlot_MarkerOnly.pdf'))
UpSet(plt$Up.lt)
# dev.off() 

plt = Up_plot(res.sum, name = names[3])
# pdf(file = paste0(plt$name,'UpSetPlot_MarkerOnly.pdf'))
UpSet(plt$Up.lt)
# dev.off() 


plt = Up_plot(res.sum, name = names[4])
#pdf(file = paste0(plt$name,'UpSetPlot_MarkerOnly.pdf'))
UpSet(plt$Up.lt)
# dev.off() 


plt = Up_plot(res.sum, name = names[5])
# pdf(file = paste0(plt$name,'UpSetPlot_MarkerOnly.pdf'))
UpSet(plt$Up.lt)
# dev.off() 

plt = Up_plot(res.sum, name = names[6])
# pdf(file = paste0(plt$name,'UpSetPlot_MarkerOnly.pdf'))
UpSet(plt$Up.lt)
# dev.off() 

plt = Up_plot(res.sum, name = names[7])
# pdf(file = paste0(plt$name,'UpSetPlot_MarkerOnly.pdf'))
UpSet(plt$Up.lt)
# dev.off() 

