edgeR.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  
  if(length(dat$covariates %>% unique())==2){
    cat('Covariate is 2-factor levels! \n')
    dat$covariates <- as.factor(dat$covariates)
  }else{
    dat$covariates <- dat$covariates
  }
  
  d <- DGEList(counts=dat$counts, group=dat$covariates)
  d <- edgeR::calcNormFactors(d) # TMM normalization
  design <- model.matrix(~dat$covariates)
  rownames(design) <- colnames(d)
  d <- estimateDisp(d, design)# also has robust = T mode
  
  if(length(table(dat$covariates))>2){ # see edgeR guide 
    # GLM likelihood ratio test
    fit <- glmFit(d,design) 
    res <- glmLRT(fit,coef=2)$table %>% mutate(fdr = p.adjust(PValue), pval = PValue) %>% dplyr::select(c('pval','fdr'))
    rownames(res) = otu.name
    
    ## GLM quassi-likelihood ratio test
    # d=DGEList(counts=dat$counts, group=dat$covariates)
    # X <- ns(dat$covariates, df=3)
    # design <- model.matrix(~ X)
    # y <- estimateDisp(d, design)
    # fit <- glmQLFit(y, design, robust=TRUE) # negative binomial GLM for each tag and produces an object of class DGEGLM with some new components
    # res <- glmQLFTest(fit, coef=2:4)$table %>% #quasi-likelihood (QL) F-test
    #   rownames_to_column('otu.id') %>% 
    #   mutate(fdr = p.adjust(PValue), pval = PValue) %>% 
    #   dplyr::select(c('otu.id','pval','fdr')) %>% column_to_rownames('otu.id')
  }else{ # see edgeR guide Casestudy 4.1
    ## classic procedure comparing 2/+ grps
    fit=exactTest(d)
    res=topTags(fit, n=nrow(dat$counts)) %>% as.data.frame()
    colnames(res)[3:4] = c('pval','fdr')
  }
  res <- res[,c('pval','fdr'), drop = F] %>% rownames_to_column('otu.id') 
  sig <- sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

edgeR.gmpr.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$gmpr.counts))) {
    rownames(dat$gmpr.counts) <- paste0('O', 1:nrow(dat$gmpr.counts))
    otu.name <- rownames(dat$gmpr.counts) 
  } else {
    otu.name <- rownames(dat$gmpr.counts) 
  }
  
  if(length(dat$gmpr.covariates %>% unique())==2){
    cat('Covariate is 2-factor levels! \n')
    dat$gmpr.covariates <- as.factor(dat$gmpr.covariates)
  }else{
    dat$gmpr.covariates <- dat$gmpr.covariates
  }
  d=DGEList(counts=dat$gmpr.counts, group=dat$gmpr.covariates)
  ## ignore normalization, none: the normalization factors are set to 1
  d=edgeR::calcNormFactors(d, method = 'none')
  design <- model.matrix(~dat$gmpr.covariates)
  d=estimateDisp(d, design)
  if(length(table(dat$gmpr.covariates))>2){
    # GLM approach
    fit <- glmFit(d,design)
    res <- glmLRT(fit,coef=2)$table %>% mutate(fdr = p.adjust(PValue), pval = PValue) %>% dplyr::select(c('pval','fdr'))
    rownames(res) = otu.name
  }else{
    ## the exact test in edgeR
    fit=exactTest(d)
    res=topTags(fit, n=nrow(dat$gmpr.counts)) %>% as.data.frame()
    colnames(res)[3:4] = c('pval','fdr')
  }
  res <- res[,c('pval','fdr'), drop = F] %>% rownames_to_column('otu.id') 
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

edgeR.Wrench.wrapper <- function(dat) { 
  t1 = Sys.time()
  # see link: https://bioconductor.org/packages/devel/bioc/vignettes/Wrench/inst/doc/vignette.html
  if (is.null(rownames(dat$Wrench.nf))) {
    rownames(dat$Wrench.nf) <- paste0('O', 1:nrow(dat$Wrench.nf))
    otu.name <- rownames(dat$Wrench.nf) 
  } else {
    otu.name <- rownames(dat$Wrench.nf) 
  }
  
  d=edgeR::DGEList(counts=dat$Wrench.nf,
                   group = as.matrix(dat$Wrench.covariates),
                   norm.factors=dat$compositionalFactors)
  ## the classic approach in edgeR
  d=estimateCommonDisp(d)
  d=estimateTagwiseDisp(d)
  fit=exactTest(d)
  res=topTags(fit, n=nrow(dat$Wrench.nf)) %>% as.data.frame()
  colnames(res)[3:4] = c('pval','fdr') 
  res <- res[,c('pval','fdr'), drop = F] %>% rownames_to_column('otu.id') 
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

metagenomeSeq2.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  mgs = newMRexperiment(counts=dat$counts)
  mgs = metagenomeSeq::cumNorm(mgs) # calculate the scaling factors, run cumNorm
  pd = pData(mgs)
  mod = model.matrix(~1 + dat$covariates, data = pd)
  if(length(table(dat$covariates))>2){
    #fitZig 
    settings = zigControl(maxit=100,verbose=TRUE)
    fit = fitZig(obj = mgs,mod=mod,control=settings) #by default, the normalizing factors for obj=MGS are included in the model
    res <- MRcoefs(fit, coef=colnames(mod)[2],by=colnames(mod)[2], number = dim(fData(mgs))[1], group=2)
  }else{
    fit = fitFeatureModel(mgs,mod) # MSeq recommned this over the MSeq1, which uses fitZig()
    # res = MRcoefs(fit)
    res = MRfulltable(fit, number = nrow(assayData(mgs)$counts)) 
  }
  colnames(res)[colnames(res) == 'adjPvalues'] = 'fdr'
  colnames(res)[colnames(res) == 'pvalues'] = 'pval'
  res <- res[,c('pval','fdr'), drop = F] %>% rownames_to_column('otu.id') 
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

metagenomeSeq.wrapper <- function(dat) {
  t1 = Sys.time()
  
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  mgs = newMRexperiment(counts=dat$counts)
  mgs = metagenomeSeq::cumNorm(mgs) # calculate the scaling factors, run cumNorm
  pd = pData(mgs)
  mod = model.matrix(~1 + dat$covariates, data = pd)
  if(length(table(dat$covariates))>2){
    #fitZig 
    settings = zigControl(maxit=100,verbose=TRUE)
    fit = fitZig(obj = mgs,mod=mod,
                 control=settings) #by default, the normalising factors for obj=MGS are included in the model
    res <- MRcoefs(fit, coef=colnames(mod)[2],by=colnames(mod)[2], number = dim(fData(mgs))[1], group=2)#
  }else{
    fit = fitZig(mgs, mod)
    res = MRfulltable(fit, number = nrow(assayData(mgs)$counts))
  }
  colnames(res)[colnames(res) == 'adjPvalues'] = 'fdr'
  colnames(res)[colnames(res) == 'pvalues'] = 'pval'
  res <- res[,c('pval','fdr'), drop = F] %>% rownames_to_column('otu.id') 
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

metagenomeSeq2.Wrench.wrapper <- function(dat) {
  t1 = Sys.time()
  
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  mgs = newMRexperiment(counts=dat$counts)
  normFactors(mgs) <- dat$normalizationFactors
  mod = model.matrix(~dat$covariates)
  fit = fitFeatureModel(mgs,mod) # MSeq recommned this over the MSeq1, which uses fitZig()
  res = MRfulltable(fit, number = nrow(assayData(mgs)$counts))
  colnames(res)[colnames(res) == 'adjPvalues'] = 'fdr'
  colnames(res)[colnames(res) == 'pvalues'] = 'pval'
  res <- res[,c('pval','fdr'), drop = F] %>% rownames_to_column('otu.id') 
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

wilcox.wrapper <- function(dat) {
  t1 = Sys.time()
  
  if (is.null(rownames(dat$prop))) {
    rownames(dat$prop) <- paste0('O', 1:nrow(dat$prop))
    otu.name <- rownames(dat$prop) 
  } else {
    otu.name <- rownames(dat$prop) 
  }
  prop <- dat$prop
  if(length(table(dat$covariates))>2){
    pval <- apply(prop, 1, function (y) cor.test(as.numeric(y), dat$covariates, method = 'spearman')$p.value)
  }else{
    cat('Covariate is 2-factor levels! \n')
    pval <- apply(prop, 1,  function (y) wilcox.test(y ~ as.factor(dat$covariates))$p.value)
  }
  
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res <- as.data.frame(t(rbind(pval, fdr)))%>% rownames_to_column('otu.id') 
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

ZicoSeq.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  grp.name <- names(dat$meta.dat)[1]
  zicoseq.obj <- ZicoSeq::ZicoSeq(meta.dat = dat$meta.dat, comm = dat$counts, grp.name = grp.name, adj.name = NULL,
                                  # Filtering criterion
                                  prev.filter = 0, min.prop = 0, abund.filter = 0, # since we have done preprocessing in the dat$counts!
                                  # Winsorization to replace outliers
                                  is.winsor = F, winsor.qt = 0.97,
                                  # Posterior sampling
                                  is.prior = T, prior.dist = c('BetaMix'), post.method = c('mean'), post.sample.no = 25, 
                                  # Link functions
                                  link.func = list(function (x) x^0.25, function (x) x^0.5, function (x) x^0.75), stats.combine.func = max,
                                  # Permutation
                                  perm.no = 99,  strata = NULL, 
                                  # Multiple stage normalization
                                  stage.no = 6, topK = NULL, stage.fdr = 0.75, stage.max.pct = 0.50,  
                                  # Tree-based FDR control and family-wise error rate control
                                  is.fwer = FALSE, is.tree.fdr = FALSE, tree = NULL, 
                                  verbose = TRUE, return.comm = FALSE, return.perm.F = FALSE)
  
  res <- cbind(as.data.frame(zicoseq.obj$p.raw),as.data.frame(zicoseq.obj$p.adj.fdr))%>% rownames_to_column('otu.id') 
  colnames(res) = c('otu.id','pval','fdr')
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

Rarefy.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  comm <- t(GUniFrac::Rarefy(t(dat$counts))$otu.tab.rff)
  if(length(table(dat$covariates))>2){
    pval <- apply(comm, 1, function (y) cor.test(as.numeric(y), dat$covariates, method = 'spearman')$p.value)
  }else{
    cat('Covariate is 2-factor levels! \n')
    pval <- apply(comm, 1,  function (y) wilcox.test(y ~ as.factor(dat$covariates))$p.value)
  }
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res <- as.data.frame(t(rbind(pval, fdr))) %>% rownames_to_column('otu.id')
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

DESeq.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  dat$counts = dat$counts + 1L # ADD 11/05/2020, for avoiding ERROR: every gene contains at least one zero, cannot compute log geometric means
  design <- as.formula(paste('~', names(dat$meta.dat[1])))
  dds <- DESeqDataSetFromMatrix(countData = dat$counts, colData = dat$meta.dat, design= design)
  dds <- DESeq2::DESeq(dds)
  res <- as.data.frame(as.matrix(results(dds)))%>% rownames_to_column('otu.id') %>% 
    mutate(pval = pvalue, fdr = p.adjust(pval, 'fdr'))%>% 
    dplyr::select(c('otu.id','pval','fdr'))
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

DESeq.gmpr.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  counts <- dat$counts[,names(dat$gmpr.size.factor)]
  meta <- dat$meta.dat[names(dat$gmpr.size.factor),, drop =F]
  design <- as.formula(paste('~', names(meta)))
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design)
  sizeFactors(dds) <- dat$gmpr.size.factor
  dds <- DESeq(dds)
  res <- as.data.frame(as.matrix(results(dds)))%>% rownames_to_column('otu.id') %>% 
    mutate(pval = pvalue, fdr = p.adjust(pval, 'fdr'))%>% 
    dplyr::select(c('otu.id','pval','fdr'))
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

DESeq.Wrench.wrapper <- function(dat) {
  t1 = Sys.time()
  # see link: https://bioconductor.org/packages/devel/bioc/vignettes/Wrench/inst/doc/vignette.html
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  design <- as.formula(paste('~', names(dat$meta.dat[1])))
  dds <- DESeqDataSetFromMatrix(countData =dat$counts, colData = dat$meta.dat, design= design)
  DESeq2::sizeFactors(dds) <- dat$normalizationFactors
  dds <- DESeq(dds)
  res <- as.data.frame(as.matrix(results(dds))) %>% rownames_to_column('otu.id') %>% mutate(pval = pvalue,fdr = p.adjust(pval, 'fdr')) %>%  dplyr::select(c('otu.id','pval','fdr'))
  sig <- sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

Aldex2.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  comm <- as.data.frame(dat$counts)
  
  if(length(table(dat$covariates))>2){
    covariate <- dat$covariates
    aldex.obj <- aldex.clr(comm, covariate)
    corr.test <- aldex.corr(aldex.obj, covariate)
    fdr <- corr.test[,'BH'] # wilcox method with BH adjustment
    names(fdr) <- rownames(corr.test)
  }else{
    covariate <- as.factor(dat$covariates)
    aldex.obj <- aldex(comm, covariate)
    fdr <- aldex.obj$wi.eBH # wilcox method with BH adjustment
    names(fdr) <- rownames(aldex.obj)
  }
  res <- aldex.obj %>% dplyr::select(c('wi.ep','wi.eBH'))%>% rownames_to_column('otu.id') 
  colnames(res)[2:3] = c('pval','fdr')
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

Omnibus.wrapper <- function(dat) { # grp.name can be categorical and numerical
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  covariate <- names(dat$meta.dat[1])
  Omnibus.obj <- ZISeq(dat$counts, dat$meta.dat, size.factor = dat$gmpr.size.factor,grp.name = covariate, method = 'omnibus',filter = F) # make sure the rownames(meta.dat) = colnames(otu.tab)
  res <- as.data.frame(Omnibus.obj$result) %>% 
    rownames_to_column('otu.id') %>%
    mutate(pval = p.value,fdr = p.adjust(pval,'fdr')) %>% 
    dplyr::select(c('otu.id','pval','fdr'))
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

ANCOMBC.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  dat$meta.dat$new <- dat$meta.dat$grp
  sam <- phyloseq::sample_data(dat$meta.dat)
  otu.tab <- otu_table(dat$counts, taxa_are_rows = T)
  phyloseq <- merge_phyloseq(otu.tab, sam)
  
  if(length(dat$covariates %>% unique())==2){
    out = ancombc(phyloseq = phyloseq, formula = names(dat$meta.dat[1]), p_adj_method = "fdr", lib_cut = 50, group = names(dat$meta.dat[1]))
  }else{
    out = ancombc(phyloseq = phyloseq, formula = names(dat$meta.dat[1]), p_adj_method = "fdr",lib_cut = 50, group = NULL, struc_zero = FALSE, neg_lb = FALSE)
  }
  
  res <- as.data.frame(out$res$p_val) %>% rownames_to_column('otu.id') %>% mutate(pval = grp,fdr = p.adjust(pval, 'fdr')) %>% dplyr::select(-c('grp','(Intercept)'))
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

RAIDA.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  df.raida <- as.data.frame(dat$meta.dat) %>% rownames_to_column('otu.id') %>% inner_join(as.data.frame(t(dat$counts))  %>% rownames_to_column('otu.id'))
  df.raida <- df.raida[order(-df.raida$grp),]
  rownames(df.raida) <- NULL
  n.lib <- c(as.numeric(table(df.raida$grp)[1]),as.numeric(table(df.raida$grp)[2]))
  c.data <- df.raida %>% dplyr::select(-'grp') %>% column_to_rownames('otu.id') %>% t(.) %>% as.data.frame(.)
  raida.obj <- raida(c.data, n.lib, show.ref.features = F, show.all.features = T)
  res <- raida.obj %>%
    rownames_to_column('otu.id') %>%
    mutate(fdr = p.adj, pval = p) %>%
    dplyr::select(otu.id, fdr, pval) 
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

dacomp.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  } # Tests of Association for Ordinal/Continuous Phenotypes
  dacomp.otu <- t(dat$counts) 
  dacomp.meta <- as.vector(dat$covariates)
  references = dacomp.select_references(X = dacomp.otu, median_SD_threshold = 0.6, verbose = T)
  
  if(length(table(dat$covariates))>2){
    dacomp.obj <- dacomp.test(X = dacomp.otu, y = dacomp.meta, ind_reference_taxa = references,
                              nr_perm = 1000, test = DACOMP.TEST.NAME.SPEARMAN, verbose = T,q = 0.05)
  }else{
    dacomp.obj <- dacomp.test(X = dacomp.otu, y = dacomp.meta, ind_reference_taxa = references,
                              nr_perm = 1000, test = DACOMP.TEST.NAME.WILCOXON, verbose = T,q = 0.05)
  }
  pval <- dacomp.obj$p.values.test
  names(pval) <- otu.name
  # Reference OTUs are NAs, shall be deleted
  ref.otu = colnames(dacomp.otu)[references$selected_references]
  pval <- pval[!(names(pval) %in% ref.otu)]
  res <- as.data.frame(pval) %>%
    rownames_to_column('otu.id') %>%
    dplyr::select(otu.id, pval) %>%
    mutate(fdr = p.adjust(pval, 'fdr')) 
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

ldm.wrapper <- function(dat) {#https://github.com/yijuanhu/LDM/
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  ldm.otu <<- t(dat$counts) %>% as.data.frame()
  ldm.obj <- ldm(formula=ldm.otu ~ grp, data = dat$meta.dat, fdr.nominal = 0.05) # otu: row sample, col taxa
  res <- as.data.frame(t(ldm.obj$p.otu.omni)) %>% # only select p-values for the OTU-specific tests based on the omnibus statistics
    rownames_to_column('otu.id') %>%
    dplyr::rename(pval = V1) %>%  mutate(fdr = p.adjust(pval, 'fdr'))
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

RioNorm.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  RioNorm.otu <- dat$counts %>% as.data.frame()
  RioNorm.class <- dat$covariates %>% as.numeric()
  # RioNorm.class <- RioNorm.class-1
  size_factor = hk_find(RioNorm.otu, min_avg_counts = 5)$size_factor
  log_sizefactor = log(size_factor)
  
  # Step 2: apply "overdisp_scoretest" function in the package for over-dispersion test
  # output will be two lists of OTUs, one with dispersed OTUs, the other with non-overdispersed OTUs
  scoretest = overdisp_scoretest(RioNorm.otu, RioNorm.class, log_sizefactor)
  ID_nondisp = scoretest$ID_nondisp
  ID_disp = scoretest$ID_disp
  
  # Step 3: apply "ZIP_test" function in the package to test differential abundance for non-overdispersed OTUs
  if(length(ID_nondisp) > 0){
    nondisp_OTU = RioNorm.otu[ID_nondisp,]
    nondisp_res = ZIP_test(nondisp_OTU, RioNorm.class, log_sizefactor)
  }else{
    nondisp_res = NULL
  }
  # Step 4: apply "ZINB_test" function in the package to test differential abundance for overdispersed OTUs
  if(length(ID_disp) > 0){
    disp_OTU = RioNorm.otu[ID_disp,]
    disp_res = ZINB_test(disp_OTU, RioNorm.class, log_sizefactor)
  }else{
    disp_res = NULL
  }
  # combine test results from ZIP and ZINB
  combined_res = apply(cbind(disp_res, nondisp_res),1,unlist) %>% as.data.frame()
  combined_res$padj = as.numeric(as.character(combined_res$padj)) # defualt is BH
  combined_res$pvalue = as.numeric(as.character(combined_res$pvalue)) 
  res <- combined_res %>% as.data.frame() %>% dplyr::rename(otu.id = id, fdr = padj)
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

glmquassi.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  comm <- dat$counts[,names(dat$gmpr.size.factor)]
  covariate <- dat$meta.dat[names(dat$gmpr.size.factor),, drop=F][,'grp']
  pvs <- list()
  tryCatch({
    for(i in 1:nrow(comm)){
      m1.op <- glm(comm[i,] ~ covariate, offset=log(dat$gmpr.size.factor),family=quasipoisson)
      coef.op <- coef(m1.op)
      pv.op <- wald.test(b = coef.op, Sigma = vcov(m1.op), Terms = 2)$result$chi2['P'] #coef(summary(m1.op))[2,4] 
      pvs[i] <- pv.op
    }}, error =function(e){cat(paste0(' ERROR : '),conditionMessage(e), "\n")})
  
  res <- melt(pvs) 
  res$L1 <- rownames(comm)
  colnames(res) <- c('pval','otu.id')
  res <- res[res$otu.id %in% otu.name,]
  res$fdr = p.adjust(res$pval, 'fdr') 
  res <- res[,c('otu.id','pval','fdr')]
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}

BBinomial.wrapper <- function(dat) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  sam <- phyloseq::sample_data(dat$meta.dat)
  otu.tab <- otu_table(dat$counts, taxa_are_rows = T)
  phyloseq <- merge_phyloseq(otu.tab, sam)
  # Scenario: Testing for differential abundance across Day, without controlling for anything else:
  corncob <- differentialTest(formula = ~ grp,
                              phi.formula = ~ 1,
                              formula_null = ~ 1,
                              phi.formula_null = ~ 1,
                              test = "Wald", boot = FALSE,
                              fdr_cutoff = 0.05,
                              data = phyloseq)
  res <- as.data.frame(cbind(corncob$p, corncob$p_fdr)) %>% rownames_to_column('otu.id')
  colnames(res)[2:3] <- c('pval','fdr')
  res <- res[res$otu.id %in% otu.name,]
  sig = sum(is.na(res$fdr)<=0.05)
  na <- sum(is.na(res$fdr))/nrow(res)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  return(list(sig = sig,na = na, time = time, res = res))
}


# RioNorm.wrapper(dat)$sig
# names <- c('wilcox','DESeq','DESeq.Wrench','DESeq.gmpr','metagenomeSeq2','metagenomeSeq2.Wrench',
#            'edgeR','edgeR.Wrench','edgeR.gmpr','ZicoSeq','Rarefy','Aldex2','Omnibus','dacomp',
#            'RAIDA','ldm','ANCOMBC','glmquassi','BBinomial')

