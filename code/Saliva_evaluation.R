setwd('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/')
files = list.files(path = 'data/',pattern = 'Rdata$')
pkg <- c('ANCOMBC','eBay','microbiome',"aod","phyloseq","reshape","corncob","MASS","readr","DESeq2", "ALDEx2", "metagenomeSeq", "edgeR", "GUniFrac", "grDevices", "dirmult", "exactRankTests","nlme", "dplyr", "magrittr", "tidyr", "protoclust", "ggplot2", "compositions","rmutil","tibble","reticulate","dacomp","LDM","Wrench","RioNorm2")
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
qmt_table = sweep(table, MARGIN=2, oral_trimmed_metadata$flowcount, `*`)
qmt_table[1:5,1:5]
qmt_prop = t(t(qmt_table) / colSums(qmt_table))
d =as.data.frame(t(qmt_prop)) %>% rownames_to_column('id')%>% inner_join(meta.dat %>% rownames_to_column('id'))
d1 = d[,c('grp',"k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Haemophilus")]
t.test(d1[,2]~d1$grp)
d1 = d[,c('grp',"k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Actinomyces")]
t.test(d1[,2]~d1$grp)


meta.dat = ((oral_trimmed_metadata) %>% column_to_rownames('#SampleID'))[,'brushing_event',drop = F]
colnames(meta.dat) = 'grp'
meta.dat$grp <- as.factor(meta.dat$grp)
covariates <- as.factor(meta.dat$grp)
names(covariates) <- rownames(meta.dat)

otu.tab.sim <- table

## raw otu table prevelance filtration
prop <- t(t(otu.tab.sim) / colSums(otu.tab.sim)) 
prop <- prop[rowSums(prop!=0) > 0.1 * ncol(prop), , drop=FALSE]
otu.tab.sim <- otu.tab.sim[rownames(prop), , drop=FALSE]
## raw otu tableminimal abundance filtration
prop <- prop[rowMaxs(prop) > 0.002, , drop=FALSE]
otu.tab.sim <- otu.tab.sim[rownames(prop), , drop=FALSE] # used for analysis
cat('4.Preprocessing finshed! \n')
dim(otu.tab.sim)
# minimal reads filtering
otu.tab.sim <- otu.tab.sim[,colSums(otu.tab.sim)>1000]
prop <- prop[,colnames(otu.tab.sim)]
meta.dat <- meta.dat[colnames(otu.tab.sim),,drop =F]
covariates <- covariates[rownames(meta.dat)]
table(covariates)

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
prop.gmpr <- prop.gmpr[rowSums(prop.gmpr!=0) > 0.1 * ncol(prop.gmpr), , drop=FALSE]
gmpr.counts <- gmpr.counts[rownames(prop.gmpr), , drop=FALSE]

## gmpr otu table minimal abundance filtration
prop.gmpr <- prop.gmpr[rowMaxs(prop.gmpr) > 0.002, , drop=FALSE]
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
  prop.Wrench <- prop.Wrench[rowSums(prop.Wrench!=0) > 0.1 * ncol(prop.Wrench), , drop=FALSE]
  Wrench.nf <- Wrench.nf[rownames(prop.Wrench), , drop=FALSE]
  
  ## gmpr otu table minimal abundance filtration
  prop.Wrench <- prop.Wrench[rowMaxs(prop.Wrench) > 0.002, , drop=FALSE]
  Wrench.nf <- Wrench.nf[rownames(prop.Wrench), , drop=FALSE] # used for analysis
  
  compositionalFactors <- compositionalFactors[colnames(Wrench.nf)]
  Wrench.covariates <- covariates[which(!is.na(compositionalFactors)) | which(!is.na(normalizationFactors))]
  Wrench.meta.dat <- meta.dat[rownames(meta.dat) %in% colnames(Wrench.nf),, drop=FALSE]
  cat('3.Wrench finshed! \n')
}

dat <- list(counts = otu.tab.sim, prop = prop, covariates = covariates, meta.dat = meta.dat,
            Wrench.nf = Wrench.nf, compositionalFactors= compositionalFactors, normalizationFactors = normalizationFactors, Wrench.covariates = Wrench.covariates,Wrench.meta.dat = Wrench.meta.dat,
            gmpr.comm=gmpr.comm, gmpr.counts = gmpr.counts, gmpr.size.factor = gmpr.size.factor, gmpr.covariates=gmpr.covariates,gmpr.meta.dat = gmpr.meta.dat)

methods_funs <- list('ZicoSeq'= "ZicoSeq.wrapper", 'Rarefy'= "Rarefy.wrapper", 'Aldex2'= "Aldex2.wrapper", 'Omnibus'="Omnibus.wrapper", 'ANCOM2'="ANCOM2.wrapper",
                     'RAIDA'="RAIDA.wrapper", 'DACOMP' = "dacomp.wrapper", 'LDM'="ldm.wrapper",'glmquassi'='glmquassi.wrapper',
                     'BBinomial' = 'BBinomial.wrapper','ANCOMBC'="ANCOMBC.wrapper", 'MSeq'="metagenomeSeq.wrapper",
                     'DESeq2'="DESeq.wrapper", 'DESeq2.Wrench'="DESeq.Wrench.wrapper", 'DESeq2.gmpr'="DESeq.gmpr.wrapper",
                     'Wilcox'='wilcox.wrapper' ,
                     'edgeR'="edgeR.wrapper", 'edgeR.gmpr'="edgeR.gmpr.wrapper", 'edgeR.Wrench'="edgeR.Wrench.wrapper",
                     'MSeq2.Wrench'="metagenomeSeq2.Wrench.wrapper", 'MSeq2'="metagenomeSeq2.wrapper",
                     'CLRBC'="CLRBC.wrapper",'eBayt'="eBayt.wrapper",'eBayW'="eBayW.wrapper",'Aldex2we'="Aldex2we.wrapper")
methods <- c('Omnibus','eBayt','eBayW','Rarefy', 'Aldex2',
             'Aldex2we', 'RAIDA', 'DACOMP', 'LDM', 'glmquassi',
             'BBinomial','ANCOMBC','Wilcox','MSeq2','MSeq2.Wrench')
res <- sig <- NULL
for(method in methods){
  tryCatch({
    cal <- calculate(dat,method = method)
    res[[method]] <- cal$res
    sig[[method]] <- cal$sig
  },error = function(e){cat(method,' ERROR : ',conditionMessage(e), "\n")})
}


multi_full <- reduce_list(res)
# save(multi_full,file = '~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/Saliva.Rdata')
load('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/Saliva.Rdata')
cutoff = 0.05
head(multi_full)
multi_full1 = multi_full %>% dplyr::select(c('taxa',colnames(multi_full)[grep('fdr',colnames(multi_full))])) %>% column_to_rownames('taxa')
multi_full1[is.na(multi_full1)] = 1
sum.df = multi_full1[rowSums(multi_full1 < cutoff),,drop=F]
sum.df = multi_full1[c("k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Haemophilus",
                       "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Actinomyces"),]
rownames(sum.df) = c('Haemophilus','Actinomyces')
colnames(sum.df) = gsub('fdr.','',colnames(sum.df))
sum.df$MSeq2 = c(1,1)
sum.df$MSeq2.Wrench = c(1,1)

sum.df2 = multi_full1[c("k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Haemophilus"),]
rownames(sum.df2) = c('Haemophilus')
sum.df2$MSeq2 = 1
sum.df2$MSeq2.Wrench =1



sig.df = sum.df[,which(sum.df[1,] <= 0.05)]


glm.res = glmquassi.wrapper0(dat,FDR = F)$res
ANCOMBC.res = ANCOMBC.wrapper0(dat,FDR = F)$res
Omnibus.res = Omnibus.wrapper0(dat,FDR = F)$res

head(Omnibus.res);head(glm.res);head(ANCOMBC.res)
res.fc = full_join(Omnibus.res %>% dplyr::select('fcs','otu.id') %>% dplyr::rename(Omnibus = fcs), glm.res %>% dplyr::select('fc','otu.id')%>% dplyr::rename(glm = fc)) %>% 
  full_join(ANCOMBC.res%>% dplyr::select('fcs','otu.id')%>% dplyr::rename(ANCOMBC = fcs)) %>% 
  filter(otu.id =="k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Haemophilus")
is.na(res.fc)

Omnibus.res %>% filter(otu.id =="k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Haemophilus")
head(res.fc)


## Since most methods are not absolute abundance
oral_trimmed_metadata <- read_delim("~/Documents/Mayo_Research/CLR_BC/data/oral_trimmed_metadata.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
oral_trimmed_metadata$flowcount = (oral_trimmed_metadata$`flow cell 5min 1` + oral_trimmed_metadata$`flow cell 5min 2`)/2

table = phyloseq::import_biom('~/Documents/Mayo_Research/CLR_BC/data/oral-collapsed-table_json.biom')
table = table@.Data[,oral_trimmed_metadata$`#SampleID`,drop=F]







