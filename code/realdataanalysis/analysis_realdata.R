setwd('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/')
files = list.files(path = 'data/',pattern = 'Rdata$')
pkg <- c('ANCOMBC',"aod",'mbzinb',"phyloseq","reshape","corncob","MASS","readr","DESeq2", "ALDEx2", "metagenomeSeq", "edgeR", "GUniFrac", "grDevices", "dirmult", "exactRankTests","nlme", "dplyr", "magrittr", "tidyr", "protoclust", "ggplot2", "compositions","rmutil","tibble","reticulate","dacomp","LDM","Wrench","RioNorm2")
lapply(pkg, require, character.only = TRUE)
source(file.path('code/', 'zeroinfl.plus.daa.R'))
source(file.path('code/', 'zeroinfl.plus.github.R'))
source('code/wrappers.R')
source('code/raida.R')
calculate <- function(dat, method){
  wrapper <- match.fun(methods_funs[[method]])
  cat(paste(method, '\n'))
  out <- wrapper(dat, FDR =F,cutoff = 0.05)
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

Result = list()
files = files[files !='LouisS_2016.Rdata']

dep.conf = NULL
for(i in 1:length(files)){
  tryCatch({
    file = files[i]
    load(file)
    dep <- sample_sums(phy)
    diff.seq.p <-cor.test(dep, as.numeric(phy@sam_data$grp), alternative = "two.sided", method = 'spearman')$p.value
    })
  if(diff.seq.p <= 0.05){
    dep.conf = c(dep.conf, file)
  }
}


for(i in 1:length(files)){
  sink(paste0('result/realdata/',gsub('.Rdata','',file),'otuput.txt'))
  tryCatch({
    file = files[i]
    load(paste0('data/curated/',file))
    dep <- sample_sums(phy)
    diff.seq.p <- cor.test(dep, as.numeric(phy@sam_data$grp), alternative = "two.sided", method = 'spearman')$p.value
    # diff.seq.p <- summary(aov(dep ~ phy@sam_data$grp))[[1]][1, 'Pr(>F)']
    if (!is.na(diff.seq.p) & diff.seq.p <= 0.05) {
      cat("Signficant sequencing depth confounding!\n")
      cat("For parametric test with sequence depth adjustment, please be cautious about the results!\n")
      cat("There may be potential residual sequence depth confounding! Let's do rarefraction!\n")
      phy = rarefy_even_depth(phy, sample.size = 2000)
    }
    
    #-- original simulated otu.tab
    # phy = prune_samples(sample_sums(phy)>=1000, phy)
    otu.tab.sim <- otu_table(phy)@.Data
    covariates <- sample_data(phy)$grp
    
    ## shuffle data

    if(all.equal(rownames(sample_data(phy)),colnames(otu.tab.sim))){
      colnames(otu.tab.sim) = names(covariates) = paste0('sample',1:ncol(otu.tab.sim))
    }
    meta.dat <- as.data.frame(covariates)
    rownames(meta.dat) <- colnames(otu.tab.sim)
    colnames(meta.dat) <- 'grp'
    
    ## raw otu table prevelance filtration
    prop <- t(t(otu.tab.sim) / colSums(otu.tab.sim)) 
    prop <- prop[rowSums(prop!=0) > 0.1 * ncol(prop), , drop=FALSE]
    otu.tab.sim <- otu.tab.sim[rownames(prop), , drop=FALSE]
    ## raw otu tableminimal abundance filtration
    prop <- prop[rowMaxs(prop) > 0.002, , drop=FALSE]
    otu.tab.sim <- otu.tab.sim[rownames(prop), , drop=FALSE] # used for analysis
    cat('4.Preprocessing finshed! \n')
    
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

    add = function(res, name){
      colnames(res)[2:3] = paste0(name,'.',colnames(res)[2:3])
      return(res)
    }
    

    methods_funs <- list('ZicoSeq'= "ZicoSeq.wrapper", 'Rarefy'= "Rarefy.wrapper", 'Aldex2'= "Aldex2.wrapper", 'Omnibus'="Omnibus.wrapper",
                         'RAIDA'="RAIDA.wrapper", 'DACOMP' = "dacomp.wrapper", 'LDM'="ldm.wrapper", 'RioNorm2'="RioNorm.wrapper",'glmquassi'='glmquassi.wrapper',
                         'BBinomial' = 'BBinomial.wrapper','ANCOMBC'="ANCOMBC.wrapper", 
                         'DESeq2'="DESeq.wrapper", 'DESeq2.Wrench'="DESeq.Wrench.wrapper", 'DESeq2.gmpr'="DESeq.gmpr.wrapper",
                         'Wilcox'='wilcox.wrapper', 'edgeR'="edgeR.wrapper", 'edgeR.gmpr'="edgeR.gmpr.wrapper", 'edgeR.Wrench'="edgeR.Wrench.wrapper",
                         'MSeq2.Wrench'="metagenomeSeq2.Wrench.wrapper", 'MSeq2'="metagenomeSeq2.wrapper",
                         'CLRBC'="CLRBC.wrapper",'eBayt'="eBayt.wrapper",'eBayW'="eBayW.wrapper",'Aldex2we'="Aldex2we.wrapper")
    
    methods <- c('eBayt','eBayW','Aldex2we','Rarefy','Aldex2','mbzinb','RAIDA','DACOMP','LDM','glmquassi','BBinomial','ANCOMBC','DESeq2','DESeq2.Wrench','DESeq2.gmpr','Wilcox','edgeR', 'edgeR.gmpr', 'edgeR.Wrench','MSeq2.Wrench', 'MSeq2')
    res <- sig <- NULL
    for(method in methods){
      tryCatch({
        cal <- calculate(dat,method = method)
        res[[method]] <- cal$res
        sig[[method]] <- cal$sig
      },error =function(e){cat(method,' ERROR : ',conditionMessage(e), "\n")})
    }
    
    res.sum <- Reduce(
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
    save(res.sum, file=paste0('result/realdata/',gsub('.Rdata','',file), "_res.Rdata"))
    Result[[file]] = res.sum 
  }, error =function(e){cat(paste0(file),conditionMessage(e), "\n")})
  sink()
}
save(Result, 'result/realdata/Result.Rdata')

