setwd('~/Documents/Mayo_Research/CLR_BC/')
setwd('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/')
files = list.files(path = 'data/',pattern = 'Rdata$')
pkg <- c('eBay','ANCOMBC',"aod",'mbzinb',"phyloseq","reshape","corncob","MASS","readr","DESeq2", "ALDEx2", "metagenomeSeq", "edgeR", "GUniFrac", "grDevices", "dirmult", "exactRankTests","nlme", "dplyr", "magrittr", "tidyr", "protoclust", "ggplot2", "compositions","rmutil","tibble","reticulate","dacomp","LDM","Wrench","RioNorm2")
lapply(pkg, require, character.only = TRUE)
source(file.path('code/', 'zeroinfl.plus.daa.R'))
source(file.path('code/', 'zeroinfl.plus.github.R'))
source('code/wrappers.R')
source('code/raida.R')

## local computer can not install ANCOMBC

calculate <- function(dat, method){
  wrapper <- match.fun(methods_funs[[method]])
  cat(paste(method, '\n'))
  out <- wrapper(dat,FDR=F, cutoff = 0.05)
  sig <- out$sig
  if(is.na(sig)){
    sig <- 0
  }
  names(sig) <- method
  res <- as.data.frame(out$res) %>% dplyr::select('otu.id','pval','fdr')
  colnames(res) <- c('taxa',paste0('pval.',method),paste0('fdr.',method))
  return(list(res = res, sig = sig))
}

process <- function(data, prev = 0.1, minp = 0.002){
  prop <- t(t(data) / colSums(data))
  prop <- prop[rowSums(prop!=0) > prev* ncol(prop), , drop=FALSE]
  data <- data[rownames(prop), , drop=FALSE]
  
  prop <- prop[rowMaxs(prop) > minp, , drop=FALSE]
  data <- data[rownames(prop), , drop=FALSE]
  return(list(counts = data, prop = prop))
}


## https://github.com/knightlab-analyses/reference-frames/tree/master/data contains two types of data: 
## 1) oral-collapsed-table.biom is Genus level collapse dat; 2) oral_trimmed_deblur.biom is the trimmed otutable
## collapse is appllied in th absolute count analysis; oral_trimmed_deblur.biom is the ASV based(trimmed) applied by all methods in the manuscript
## This is the Genus collapse file: differentials.csv
oral_trimmed_metadata <- read_delim("~/Documents/Mayo_Research/CLR_BC/data/oral_trimmed_metadata.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
oral_trimmed_metadata$flowcount = (oral_trimmed_metadata$`flow cell 5min 1` + oral_trimmed_metadata$`flow cell 5min 2`)/2

table = phyloseq::import_biom('~/Documents/Mayo_Research/CLR_BC/data/oral-collapsed-table_json.biom')
table = table@.Data[,oral_trimmed_metadata$`#SampleID`,drop=F]

## absolute count
# qmt_table = sweep(table, MARGIN=2, oral_trimmed_metadata$flowcount, `*`)
# qmt_table[1:5,1:5]
# qmt_prop = t(t(qmt_table) / colSums(qmt_table))
# d =as.data.frame(t(qmt_prop)) %>% rownames_to_column('id')%>% inner_join(meta.dat %>% rownames_to_column('id'))
# d1 = d[,c('grp',"k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Haemophilus")]
# t.test(d1[,2]~d1$grp)
# d1 = d[,c('grp',"k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Actinomyces")]
# t.test(d1[,2]~d1$grp)

meta.dat = ((oral_trimmed_metadata) %>% column_to_rownames('#SampleID'))[,'brushing_event',drop = F]
colnames(meta.dat) = 'grp'

covariates <- as.factor(meta.dat$grp)
names(covariates) <- rownames(meta.dat)
minp = 0.002;prev = 0.1
##-- Preprocessing
##--- Normalization
## GMPR
otu.tab.sim = table
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
  cat('4.Wrench ERROR! \n')
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
  # Wrench.truth <- truth[truth$otu.id %in% rownames(Wrench.nf),]
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

methods_funs <-  list('ZicoSeq'= "ZicoSeq.wrapper", 'Rarefy'= "Rarefy.wrapper", 'Aldex2'= "Aldex2.wrapper", 'Omnibus'="Omnibus.wrapper", 'ANCOM2'="ANCOM2.wrapper",
                      'RAIDA'="RAIDA.wrapper", 'DACOMP' = "dacomp.wrapper", 'LDM'="ldm.wrapper", 'RioNorm2'="RioNorm.wrapper",'glmquassi'='glmquassi.wrapper',
                      'BBinomial' = 'BBinomial.wrapper','ANCOMBC'="ANCOMBC.wrapper", 'MSeq'="metagenomeSeq.wrapper",
                      'DESeq2'="DESeq.wrapper", 'DESeq2.Wrench'="DESeq.Wrench.wrapper", 'DESeq2.gmpr'="DESeq.gmpr.wrapper",
                      'Wilcox'='wilcox.wrapper' , 'Wilcox.Wrench'='wilcox.Wrench.wrapper' , 'Wilcox.gmpr'='wilcox.gmpr.wrapper' ,
                      'edgeR'="edgeR.wrapper", 'edgeR.gmpr'="edgeR.gmpr.wrapper", 'edgeR.Wrench'="edgeR.Wrench.wrapper",
                      'MSeq2.Wrench'="metagenomeSeq2.Wrench.wrapper", 'MSeq2'="metagenomeSeq2.wrapper",'MSeq2.gmpr'="metagenomeSeq2.gmpr.wrapper",
                      'CLRBC'="CLRBC.wrapper",'eBayt'="eBayt.wrapper",'eBayW'="eBayW.wrapper",'Aldex2we'="Aldex2we.wrapper",
                      'ttest' ='ttest.wrapper','ttest.gmpr' = 'ttest.gmpr.wrapper','ttest.Wrench' = 'ttest.Wrench.wrapper','mbzinb' = "mbzinb.wrapper","Rarefyttest"="Rarefyttest.wrapper")

methods <- c('mbzinb','MSeq2','eBayW','eBayt',
             'Rarefy', 'Wilcox', "Rarefyttest", 'ttest',
             'Aldex2','Aldex2we','RAIDA', 'DACOMP', 
             'LDM', 'glmquassi','BBinomial','ANCOMBC',
             'MSeq2.Wrench')
res <- sig <- NULL
for(method in methods){
  tryCatch({
    cal <- calculate(dat,method = method)
    res[[method]] <- cal$res
    sig[[method]] <- cal$sig
  },error = function(e){cat(method,' ERROR : ',conditionMessage(e), "\n")})
}

# save(res,file = 'saliva_res.Rata')
# save(sig,file = 'saliva_sig.Rata')

# Unable to run on local for some methods, result are retrieved from cluster
# setwd('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/')
load('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/realdataanalysis/saliva_res.Rdata')
load('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/realdataanalysis/saliva_sig.Rdata')
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

multi_full1 <- multi_full[,c('taxa',colnames(multi_full)[grep('fdr',colnames(multi_full))])]
head(multi_full1)
multi_full1[is.na(multi_full1)] = 1
multi_full1 = multi_full1 %>% column_to_rownames('taxa')
sum.df = multi_full1[rowSums(multi_full1 < 0.05),]
sum.df = multi_full1[c("k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Haemophilus",
                       "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Actinomyces"),]
colnames(sum.df)[grep('fdr',colnames(sum.df))]
sum.df = sum.df[,c(colnames(sum.df)[grep('pval|fdr',colnames(sum.df))]),drop = F]
colnames(sum.df) = gsub('fdr.','',colnames(sum.df) )

Haemophilus = sum.df[1,] %>% melt()
Actinomyces = sum.df[2,] %>% melt()

df = melt(sum.df %>% rownames_to_column('taxa'))
df$taxa = gsub('.*g__','',df$taxa)
df$variable = gsub('pval\\.|fdr\\.','',df$variable)
head(df)

df$variable = as.character(df$variable)
df$variable[df$variable =='Wilcox.gmpr'] = 'GMPR+Wilcoxon'
df$variable[df$variable =='Wilcox.Wrench'] = 'Wrench+Wilcoxon'
df$variable[df$variable =='Wilcox'] = 'TSS+Wilcoxon'
df$variable[df$variable =='Wilcox.gmpr'] = 'GMPR+Spearman'
df$variable[df$variable =='Wilcox'] = 'TSS+Spearman'
df$variable[df$variable =='ttest.gmpr'] = 'GMPR+t-test'
df$variable[df$variable =='ttest.Wrench'] = 'Wrench+t-test'
df$variable[df$variable =='Rarefyttest'] = 'Rarefy+t-test'
df$variable[df$variable =='ttest'] = 'TSS+t-test'
df$variable[df$variable =='Rarefy'] = 'Rarefy+Wilcoxon'
df$variable[df$variable =='Rarefy'] = 'Rarefy+Spearman'
df$variable = gsub('ANCOMBC','ANCOM-BC',df$variable)
df$variable = gsub('glmquassi','GLM(quasipoisson)',df$variable)
df$variable = gsub('DESeq2.gmpr','GMPR+DESeq2',df$variable)
df$variable = gsub('DESeq2.Wrench','Wrench+DESeq2',df$variable)
df$variable = gsub('edgeR.gmpr','GMPR+edgeR',df$variable)
df$variable = gsub('edgeR.Wrench','Wrench+edgeR',df$variable)
df$variable = gsub('MSeq2.Wrench','Wrench+metagenomeSeq',df$variable)
df$variable = gsub('MSeq2','metagenomeSeq',df$variable)
df$variable = gsub('eBayW','eBay(Wilcoxon)',df$variable)
df$variable = gsub('eBayt','eBay(t-test)',df$variable)
df$variable = gsub('BBinomial','Beta-binomial',df$variable)
df$variable = gsub('Aldex2we','Aldex2(t-test)',df$variable)
df$variable = gsub('^Aldex2$','Aldex2(Wilcoxon)',df$variable)
df$variable = gsub('^CLRBC$','LinDA',df$variable)
df0 = df[df$taxa =='Haemophilus',]
df0 = df0[order(df0$value),]
df0
df1 <- within(df, variable <- factor(variable, levels=df0$variable))

ggplot(df1, aes(x = variable, y = value, fill = taxa)) + 
  geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7) +
  scale_fill_brewer(palette = 'Set1') + 
  geom_hline(yintercept = 0.05, colour = brewer.pal(8, 'Set1')[3], linetype = 'dashed', size = 1) +
  geom_hline(yintercept = 0.1, colour = brewer.pal(8, 'Set2')[7], linetype = 'dashed', size = 1) + 
  theme_bw()+
  theme(axis.text.x = element_text(color="black", size =26, angle = 90, vjust = 0.1, hjust = 1),
        axis.text.y = element_text(color="black", size = 26),
        legend.text = element_text(size = 26),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(color="black", size = 26)) +
  labs(y = 'FDR-adjusted p value', fill = '')
ggsave(filename = '/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/salivaFDR.pdf', width = 10, height = 8, dpi = 100)
ggplot(df1 %>% filter(taxa =='Haemophilus'), aes(x = variable, y = value, fill = taxa)) + 
  geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7) +
  scale_fill_brewer(palette = 'Set1') + 
  geom_hline(yintercept = 0.05, colour = brewer.pal(8, 'Set1')[3], linetype = 'dashed', size = 1) +
  geom_hline(yintercept = 0.1, colour = brewer.pal(8, 'Set2')[7], linetype = 'dashed', size = 1) + 
  theme_bw()+
  theme(axis.text.x = element_text(color="black", size =26, angle = 90, vjust = 0.1, hjust = 1),
        axis.text.y = element_text(color="black", size = 26),
        legend.text = element_text(size = 26),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(color="black", size = 26)) +
  labs(y = 'FDR-adjusted p value', fill = '')
ggsave(filename = '/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/salivaFDR_Haemophilus.pdf', width = 10, height = 7, dpi = 100)
  

