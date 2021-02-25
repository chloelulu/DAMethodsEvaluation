edgeR.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

edgeR.gmpr.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
  res <- res[,c('pval','fdr')] %>% rownames_to_column('otu.id')
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

edgeR.Wrench.wrapper <- function(dat, cutoff = 0.05, FDR = T) { 
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
  res <- res[,c('pval','fdr')] %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

metagenomeSeq2.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
  res <- res[,c('pval','fdr')] %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

metagenomeSeq.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
  res <- res[,c('pval','fdr')] %>% rownames_to_column('otu.id')
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

metagenomeSeq2.Wrench.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
  res <- res[,c('pval','fdr')] %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

wilcox.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
  res <- cbind(pval = pval, fdr = fdr) %>% as.data.frame() %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}

wilcox.gmpr.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if (is.null(rownames(dat$gmpr.counts))) {
    rownames(dat$gmpr.counts) <- paste0('O', 1:nrow(dat$gmpr.counts))
    otu.name <- rownames(dat$gmpr.counts) 
  } else {
    otu.name <- rownames(dat$gmpr.counts) 
  }
  prop <- t(t(dat$gmpr.counts)/colSums(dat$gmpr.counts))
  
  if(length(table(dat$covariates))>2){
    pval <- apply(prop, 1, function (y) cor.test(as.numeric(y), dat$gmpr.covariates, method = 'spearman')$p.value)
  }else{
    cat('Covariate is 2-factor levels! \n')
    pval <- apply(prop, 1,  function (y) wilcox.test(y ~ as.factor(dat$gmpr.covariates))$p.value)
  }
  
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res <- cbind(pval = pval, fdr = fdr) %>% as.data.frame() %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$gmpr.truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}

wilcox.Wrench.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if (is.null(rownames(dat$Wrench.nf))) {
    rownames(dat$Wrench.nf) <- paste0('O', 1:nrow(dat$Wrench.nf))
    otu.name <- rownames(dat$Wrench.nf) 
  } else {
    otu.name <- rownames(dat$Wrench.nf) 
  }
  prop <- t(t(dat$Wrench.nf)/colSums(dat$Wrench.nf))
  
  if(length(table(dat$covariates))>2){
    pval <- apply(prop, 1, function (y) cor.test(as.numeric(y), dat$gmpr.covariates, method = 'spearman')$p.value)
  }else{
    cat('Covariate is 2-factor levels! \n')
    pval <- apply(prop, 1,  function (y) wilcox.test(y ~ as.factor(dat$Wrench.covariates))$p.value)
  }
  
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res <- cbind(pval = pval, fdr = fdr) %>% as.data.frame() %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$Wrench.truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}

ttest.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  
  if (is.null(rownames(dat$prop))) {
    rownames(dat$prop) <- paste0('O', 1:nrow(dat$prop))
    otu.name <- rownames(dat$prop) 
  } else {
    otu.name <- rownames(dat$prop) 
  }
  prop <- dat$prop
  
  pval <- apply(prop, 1,  function (y) t.test(y ~ as.factor(dat$covariates))$p.value)

  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res <- cbind(pval = pval, fdr = fdr) %>% as.data.frame() %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}

ttest.gmpr.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  
  if (is.null(rownames(dat$gmpr.counts))) {
    rownames(dat$gmpr.counts) <- paste0('O', 1:nrow(dat$gmpr.counts))
    otu.name <- rownames(dat$gmpr.counts) 
  } else {
    otu.name <- rownames(dat$gmpr.counts) 
  }
  prop <- t(t(dat$gmpr.counts)/colSums(dat$gmpr.counts))
  
  pval <- apply(prop, 1,  function (y) t.test(y ~ as.factor(dat$gmpr.covariates))$p.value)

  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res <- cbind(pval = pval, fdr = fdr) %>% as.data.frame() %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$gmpr.truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}

ttest.Wrench.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  
  if (is.null(rownames(dat$Wrench.nf))) {
    rownames(dat$Wrench.nf) <- paste0('O', 1:nrow(dat$Wrench.nf))
    otu.name <- rownames(dat$Wrench.nf) 
  } else {
    otu.name <- rownames(dat$Wrench.nf) 
  }
  prop <- t(t(dat$Wrench.nf)/colSums(dat$Wrench.nf))
  pval <- apply(prop, 1,  function (y) t.test(y ~ as.factor(dat$Wrench.covariates))$p.value)
  
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res <- cbind(pval = pval, fdr = fdr) %>% as.data.frame() %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$Wrench.truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}

RioNorm.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}


ZicoSeq.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
    
    fdr <- zicoseq.obj$p.adj.fdr
    res <- as.data.frame(cbind(fdr = fdr)) %>% rownames_to_column('otu.id') 
    
    na <- sum(is.na(res$fdr))/nrow(res)
    sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
    time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
    
    if(FDR){
      res <- res %>% inner_join(dat$truth, by = 'otu.id')
      tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
      tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
      fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
      fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
      tpr <- ifelse(!tp, 0, tp/(tp+fn))
      fpr <- ifelse(!fp, 0, fp/(tn+fp))
      fdr <- ifelse(!fp, 0, fp/(tp+fp))
      F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
      MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
      return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
    }else{
      return(list(res = res, na = na, time = time, sig = sig))
    }
    
}

Rarefyttest.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  comm <- t(GUniFrac::Rarefy(t(dat$counts))$otu.tab.rff)
  if(length(table(dat$covariates))>2){
    cat('t-test does not support continuous variable! \n')
    # pval <- apply(comm, 1, function (y) cor.test(as.numeric(y), dat$covariates, method = 'spearman')$p.value)
  }else{
    cat('Covariate is 2-factor levels! \n')
    pval <- apply(comm, 1,  function (y) t.test(y ~ as.factor(dat$covariates))$p.value)
  }
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

Rarefy.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

DESeq.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
  res <- as.data.frame(as.matrix(results(dds)))
  pval <- res$pvalue
  names(pval) <- rownames(res)
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

DESeq.gmpr.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  
  design <- as.formula(paste('~', names(dat$meta.dat[1])))
  dds <- DESeqDataSetFromMatrix(countData = dat$counts, colData = dat$meta.dat, design = design)
  sizeFactors(dds) <- dat$gmpr.size.factor
  dds <- DESeq(dds)
  res <- as.data.frame(as.matrix(results(dds)))
  pval <- res$pvalue
  names(pval) <- rownames(res)
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

DESeq.Wrench.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
  res <- as.data.frame(as.matrix(results(dds)))
  pval <- res$pvalue
  names(pval) <- rownames(res)
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

Aldex2.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
    ## for package version 1.22.0
    # pval <- corr.test$spearman.ep
    # fdr <- corr.test[,'spearman.eBH'] 
    # names(fdr) <- names(pval)<- rownames(corr.test)
    ## for package version 1.18.0
    res <- as.data.frame(corr.test[,c('p','BH')])%>% rownames_to_column('otu.id') 
  }else{
    covariate <- as.factor(dat$covariates)
    aldex.obj <- aldex(comm, covariate)
    pval <- aldex.obj$wi.ep
    fdr <- aldex.obj$wi.eBH # wilcox method with BH adjustment
    names(fdr) <- names(pval) <- rownames(aldex.obj)
    res <- as.data.frame(cbind(pval, fdr)) %>% rownames_to_column('otu.id') 
  }
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

Aldex2we.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
    fdr <- corr.test[,'spearman.eBH'] 
    names(fdr) <- rownames(corr.test)
  }else{
    covariate <- as.factor(dat$covariates)
    aldex.obj <- aldex(comm, covariate)
    pval <- aldex.obj$we.ep
    fdr <- aldex.obj$we.eBH # welch method with BH adjustment
    names(fdr) <- names(pval) <- rownames(aldex.obj)
  }
  res <- as.data.frame(cbind(pval, fdr)) %>% rownames_to_column('otu.id') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

Omnibus.wrapper <- function(dat, cutoff = 0.05, FDR = T) { # grp.name can be categorical and numerical
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  covariate <- names(dat$meta.dat[1])
  # For match meta.dat and gmpr.size.factor
  counts <- dat$counts[,names(dat$gmpr.size.factor)]
  meta <- dat$meta[names(dat$gmpr.size.factor),,drop=F]
  Omnibus.obj <- ZISeq(counts,meta, size.factor = dat$gmpr.size.factor,grp.name = covariate, method = 'omnibus', filter = F, winsor =T)
  #Omnibus.obj <- ZISeq(dat$counts, dat$meta.dat, size.factor = dat$gmpr.size.factor,
  #                     grp.name = covariate, method = 'omnibus') # make sure the rownames(meta.dat) = colnames(otu.tab)
  #Omnibus.obj <- ZISeq(dat$counts, dat$meta.dat, size.factor =NULL, grp.name = covariate, method = 'omnibus',winsor = T,filter = F)
  res <- as.data.frame(Omnibus.obj$result) %>% 
    rownames_to_column('otu.id') %>%
    mutate(pval = p.value) %>% 
    dplyr::select(pval, otu.id) %>% mutate(fdr = p.adjust(pval,'fdr'))   
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

mbzinb.wrapper <- function(dat, cutoff = 0.05, FDR = T) { # grp.name can be categorical and numerical
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  
  mbzinb.data <- mbzinb.dataset(dat$gmpr.comm, dat$gmpr.meta.dat)
  mbzinb.data <- mbzinb::norm.factors(mbzinb.data,  norm.factors = dat$gmpr.size.factor)
  mbzinb.obj <- mbzinb.test(mbzinb.data, group = 'grp',filter.if.unfiltered =F)
  res <- mbzinb.results(mbzinb.obj,nreturn = nrow(dat$counts)) %>% rownames_to_column('otu.id') %>% mutate(pval = PValue, fdr = Padj) %>% dplyr::select(otu.id, pval, fdr)
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}

ANCOMBC.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
    out = ancombc(phyloseq = phyloseq, formula = names(dat$meta.dat[1]), p_adj_method = "fdr", lib_cut = 0, zero_cut = 1,group = names(dat$meta.dat[1]), struc_zero =T, conserve = T)
    #out = ancombc(phyloseq = phyloseq, formula = names(dat$meta.dat)[1], p_adj_method = "fdr", lib_cut = 800, group = names(dat$meta.dat[1]))
  }else{
    out = ancombc(phyloseq = phyloseq, formula = names(dat$meta.dat)[1], p_adj_method = "fdr",lib_cut = 0, zero_cut = 1, group = NULL, struc_zero = FALSE, neg_lb = FALSE,conserve = T)
  }
  res <- as.data.frame(as.matrix(cbind(pval = out$res$p_val$grp,fdr = out$res$q_val$grp))) # Caution: the code adjusted since the cluster makes errors
  colnames(res) = c('pval','fdr')
  res$otu.id <- rownames(out$res$p_val)
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}

RAIDA.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  df.raida <- as.data.frame(dat$meta.dat) %>%
    rownames_to_column('otu.id') %>%
    inner_join(as.data.frame(t(dat$counts))  %>% rownames_to_column('otu.id'))
  df.raida <- df.raida[order(-df.raida$grp),]
  rownames(df.raida) <- NULL
  n.lib <- c(as.numeric(table(df.raida$grp)[1]),as.numeric(table(df.raida$grp)[2]))
  c.data <- df.raida %>% dplyr::select(-'grp') %>% column_to_rownames('otu.id') %>% t(.) %>% as.data.frame(.)
  
  # c.data <- as.data.frame(dat$counts)
  # n.lib <- c(as.numeric(table(dat$meta.dat)[1]),as.numeric(table(dat$meta.dat)[2]))
  # raida.obj <- raida(c.data, n.lib, show.ref.features = F, show.all.features = T)
  
  raida.obj <- raida(c.data, n.lib, show.ref.features = T, show.all.features = T)
  res <- raida.obj$result %>%
    rownames_to_column('otu.id') %>%
    mutate(pval = p, fdr = p.adj) %>%
    dplyr::select(otu.id, pval, fdr)
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  # add for references features are regarded as non-DAs, thus it is not wrong theoratically
  res <- rbind(res,as.data.frame(cbind(otu.id = raida.obj$reference.features, pval = rep(1,length(raida.obj$reference.features)),
                            fdr = rep(1,length(raida.obj$reference.features)))))
  res$pval <- as.numeric(as.character(res$pval))
  res$fdr <- as.numeric(as.character(res$fdr))
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

dacomp.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
                              nr_perm = 1000, test = DACOMP.TEST.NAME.SPEARMAN, verbose = T,q = cutoff,
                              disable_DSFDR=T)
  }else{
    dacomp.obj <- dacomp.test(X = dacomp.otu, y = dacomp.meta, ind_reference_taxa = references,
                              nr_perm = 1000, test = DACOMP.TEST.NAME.WILCOXON, verbose = T,q = cutoff,
                              disable_DSFDR=T)
  }
  pval <- dacomp.obj$p.values.test
  names(pval) <- otu.name
  qval <- dacomp.obj$p.values.test.adjusted
  names(qval) <- otu.name
  ref.otu = colnames(dacomp.otu)[references$selected_references]
  # replace reference otu with 1, while NA is still NA
  pval[names(pval) %in% ref.otu] = 1
  qval[names(qval) %in% ref.otu] = 1
  res <- as.data.frame(cbind(pval, fdr=qval)) %>% rownames_to_column('otu.id')
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

ldm.wrapper <- function(dat, cutoff = 0.05, FDR = T) {#https://github.com/yijuanhu/LDM/
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  ldm.otu <<- t(dat$counts) %>% as.data.frame()
  ldm.obj <- ldm(formula=ldm.otu ~ grp, data = dat$meta.dat, fdr.nominal = cutoff) # otu: row sample, col taxa
  res <- as.data.frame(t(ldm.obj$p.otu.omni)) %>% # only select p-values for the OTU-specific tests based on the omnibus statistics
    rownames_to_column('otu.id') %>%
    dplyr::rename(pval = V1) %>% mutate(fdr = p.adjust(pval, 'fdr'))
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}

RioNorm.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
  
}

glmquassi.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  comm <- dat$counts[,names(dat$gmpr.size.factor)]
  covariate <- dat$meta.dat[names(dat$gmpr.size.factor),]
  pvs <- NULL
  tryCatch({
    for(i in 1:nrow(comm)){
    m1.op <- glm(comm[i,] ~ covariate, offset=log(dat$gmpr.size.factor),family=quasipoisson)
    pv.op <- wald.test(b = coef(m1.op), Sigma = vcov(m1.op), Terms = 2)$result$chi2['P'] #coef(summary(m1.op))[2,4] 
    pvs <- c(pvs,pv.op)
  }}, error =function(e){cat(paste0(' ERROR : '),conditionMessage(e), "\n")})
  
  res <- as.data.frame(pvs)
  res$otu.id <- rownames(comm)
  colnames(res) <- c('pval','otu.id')
  res <- res[res$otu.id %in% otu.name,]
  res$fdr = p.adjust(res$pval, 'fdr') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
    }
}

BBinomial.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
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
                   fdr_cutoff = cutoff,
                   data = phyloseq)
  res <- as.data.frame(cbind(corncob$p,corncob$p_fdr)) %>% rownames_to_column('otu.id')
  colnames(res)[2:3] <- c('pval','fdr')
  res <- res[res$otu.id %in% otu.name,]
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}

## ADD 12/31/2020
eBayW.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  
  comm <- as.data.frame(t(dat$counts))
  covariate <- as.factor(dat$covariates)
  levels(covariate) <- c(0,1) # package require the case = 1, control = 0. If use other value, will cause fail
  
  eBay.res <- eBay(otu.data=comm, group=covariate, cutf=cutoff, test.methods='wilcoxon', adj.m='BH')
  res <- as.data.frame(cbind(pval = rep(1,length(eBay.res$final.p)),fdr = eBay.res$final.p)) %>% rownames_to_column('otu.id')
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}

eBayt.wrapper <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  
  comm <- as.data.frame(t(dat$counts))
  covariate <- as.factor(dat$covariates)
  levels(covariate) <- c(0,1) # package require the case = 1, control = 0. If use other value, will cause fail
  
  eBay.res <- eBay(otu.data=comm, group=covariate, cutf=cutoff, test.methods='t', adj.m='BH')
  res <- as.data.frame(cbind(pval = rep(1,length(eBay.res$final.p)),fdr = eBay.res$final.p)) %>% rownames_to_column('otu.id')
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}


## add 12/28/2020
CLRBC.wrapper <- function(dat, cutoff = 0.05, FDR = T){
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  Y <- dat$counts
  n <- nrow(dat$meta.dat)
  Z <- cbind(as.factor(dat$meta.dat$grp), rep(1 , n)) 
  CLR_BC <- daa.fun(Y, Z)
  res <- as.data.frame(CLR_BC) %>% mutate(otu.id =rownames(Y))
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}


# load('/Users/m216453/Documents/Mayo_Research/SemiSimulation/0925/sim/new/sim/iter1/D3abundantloglinearSub_L1nOTU_L5binaryDL3balancednonenSam_L1nonemedium.Rdata')
## LinDA with BBmix imputation
LINDA.wrapper <- function(dat, cutoff = 0.05, FDR = T){
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  Y <- dat$counts
  Z <- cbind(as.numeric(dat$meta.dat$grp), rep(1 , nrow(dat$meta.dat))) 
  CLR_BC <- daa.mix(Y, Z)
  res <- as.data.frame(CLR_BC) %>% mutate(otu.id =rownames(Y))
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}



## temp store Jun's new method
daa.fun <- function(Y, Z, impute = T, test = T) {
  ## library(modeest)
  ## Y: The count matrix 
  ## Z: The annotation matrix with n row samples and d columns. Make sure the first column 
  ##    is the variable of interest and there should be a constant (all ones) column
  ## return: p-values
  
  if (test) {
    corr.pval <- cor.test(Z[, 1], colSums(Y), alternative = "two.sided", method = 'spearman')$p.value
    if (!is.na(corr.pval) & corr.pval <= 0.1) impute <- TRUE 
    else impute <- FALSE
  } 
  
  if (sum(Y == 0) == 0) {
    YY <- Y
  } else if (impute) {
    m <- nrow(Y)
    N <- colSums(Y)
    N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
    ind <- which(Y > 0)
    N.mat[ind] <- 0
    tmp <- N[max.col(N.mat)]
    YY <- Y + N.mat / tmp 
  } else {
    YY <- Y + 0.5
  }
  
  
  n <- ncol(Y)
  logY <- log(YY)
  W <- t(logY) - colMeans(logY)
  
  ZtZ.inv <- solve(t(Z) %*% Z)
  para.est <- ZtZ.inv %*% t(Z) %*% W
  W.est <- Z %*% para.est
  r2 <- (W - W.est) ^ 2
  sigma2.est <- colSums(r2) / (n - ncol(Z))
  
  z <- para.est[1, ]
  mode <- mlv(sqrt(n) * z, method = "meanshift", kernel = 'gaussian') / sqrt(n)
  
  tmp <- sqrt(sigma2.est * ZtZ.inv[1, 1])
  T.hat <- (z - mode) / tmp
  pval <- as.vector(2 * pt(-abs(T.hat), n - ncol(Z)))
  fdr <- p.adjust(pval, method = 'fdr')
  res <- cbind(pval, fdr)
  return(res)
}



daa.mix <- function(Y, Z, impute = FALSE, test = T) {
  ## library(modeest)
  ## Y: The count matrix 
  ## Z: The annotation matrix with n row samples and d columns. Make sure the first column 
  ##    is the variable of interest and there should be a constant (all ones) column
  ## return: p-values
  
  if (test) {
    corr.pval <- cor.test(Z[, 1], colSums(Y), alternative = "two.sided", method = 'spearman')$p.value
    if (!is.na(corr.pval) & corr.pval <= 0.1) impute <- TRUE 
    else impute <- FALSE
  } 
  
  if (sum(Y == 0) == 0) {
    YY <- Y
  } else if (impute) {
    depth <- colSums(Y)
    sample.no <- ncol(Y)
    Y.p <- apply(Y, 1, function (x) {
      
      err1 <- try(res <- bbmix.fit.MM(x, depth))
      
      # Handle error
      if (class(err1) != 'try-error') {
        
        prop1.1 <- rbeta(sample.no, shape1 = x + res$shape1.1, shape2 = res$shape1.2 + depth - x)
        prop1.2 <- rbeta(sample.no, shape1 = x + res$shape2.1, shape2 = res$shape2.2 + depth - x)
        prop <- ifelse(runif(sample.no) <= res$q1, prop1.1, prop1.2)
        
      } else {
        prop <- x / depth
        v <- var(prop)
        m <- mean(prop)
        
        a1 <- ((1 - m) / v - 1 / m) * m ^ 2
        a2 <- a1 * (1 / m - 1)
        post.sample.no = 1
        
        if (is.na(a1) | a1 < 0) {
          # uniform prior
          prop <- rbeta(sample.no * post.sample.no, shape1 = x + 1, shape2 = otu.no + depth - x)
        } else {
          prop <- rbeta(sample.no * post.sample.no, shape1 = x + a1, shape2 = a2 + depth - x)
        }	
      }
      return(prop)
    })
    rowSums(Y.p)
    Y.p[Y.p < 1e-6] <- 1e-6
    dim(Y.p)
    YY <- t(Y.p / rowSums(Y.p))
  } else {
    YY <- Y + 0.5
  }
  n <- ncol(Y)
  logY <- log(YY)
  W <- t(logY) - colMeans(logY)
  
  ZtZ.inv <- solve(t(Z) %*% Z)
  para.est <- ZtZ.inv %*% t(Z) %*% W
  W.est <- Z %*% para.est
  r2 <- (W - W.est) ^ 2
  sigma2.est <- colSums(r2) / (n - ncol(Z))
  
  z <- para.est[1, ]
  mode <- mlv(sqrt(n) * z, method = "meanshift", kernel = 'gaussian') / sqrt(n)
  
  tmp <- sqrt(sigma2.est * ZtZ.inv[1, 1])
  T.hat <- (z - mode) / tmp
  pval <- as.vector(2 * pt(-abs(T.hat), n - ncol(Z)))
  fdr <- p.adjust(pval, method = 'fdr')
  res <- cbind(pval, fdr)
  
  return(res)
}


bbmix.fit.MM <- function (ct, dep,  nIter = 10, winsor.qt = 1.0) {
  
  if (mean(ct == 0) >= 0.95 | sum(ct != 0) < 5) {
    stop('The number of nonzero is too small to fit the model! Consider removing it from testing!\n')
  }
  
  # Initialization
  prop0 <- ct / dep
  qt <- quantile(prop0, winsor.qt)
  prop0[prop0 >= qt] <- qt
  
  
  # This intialization may not be appropriate. You may consider testing other strategy.
  var1 <-  var(prop0[prop0 <= median(prop0)])
  mean1 <- mean(prop0) / 2
  
  var2 <- var(prop0[prop0 > median(prop0)])
  mean2 <- 3 * mean(prop0) / 2
  
  pi <- 0.5
  for (i in 1 : nIter) {
    shape1.1 <- ((1 - mean1) / var1 - 1 / mean1) * mean1 ^ 2
    shape1.2 <- shape1.1 * (1 / mean1 - 1)
    
    shape2.1 <- ((1 - mean2) / var2 - 1 / mean2) * mean2 ^ 2
    shape2.2 <- shape2.1 * (1 / mean2 - 1)
    
    m1 <- shape1.1 / (shape1.1 + shape1.2)
    s1 <- shape1.1 + shape1.2
    
    m2 <- shape2.1 / (shape2.1 + shape2.2)
    s2 <- shape2.1 + shape2.2
    
    f1 <- dbetabinom(ct, dep, m1, s1)
    f2 <- dbetabinom(ct, dep, m2, s2)
    
    q1 <-  pi * f1 /  (pi * f1 + (1 - pi) * f2)
    q2 <-  1 - q1
    
    pi <- mean(q1)
    
    # Rough estimation
    mean1 <- sum(prop0 * q1) / sum(q1)
    var1 <- sum((prop0 - mean1)^2 * q1) / sum(q1)
    
    mean2 <- sum(prop0 * q2) / sum(q2)
    var2 <- sum((prop0 - mean2)^2 * q2) / sum(q2)
  }
  
  return(list(shape1.1 = shape1.1, shape1.2 = shape1.2, shape2.1 = shape2.1, shape2.2 = shape2.2, pi = pi, q1 = q1))
  
}

glmquassi.wrapper0 <- function(dat, cutoff = 0.05, FDR = T) {
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  comm <- dat$counts[,names(dat$gmpr.size.factor)]
  covariate <- dat$meta.dat[names(dat$gmpr.size.factor),]
  pvs <- fcs <- NULL
  tryCatch({
    for(i in 1:nrow(comm)){
      m1.op <- glm(comm[i,] ~ covariate, offset=log(dat$gmpr.size.factor),family=quasipoisson)
      pv.op <- wald.test(b = coef(m1.op), Sigma = vcov(m1.op), Terms = 2)$result$chi2['P'] #coef(summary(m1.op))[2,4] 
      pvs <- c(pvs,pv.op)
      coef.op <- coef(m1.op)		
      fcs <- c(fcs,coef.op[2])
    }}, error =function(e){cat(paste0(' ERROR : '),conditionMessage(e), "\n")})
  
  res <- as.data.frame(cbind(pvs,fcs))
  res$otu.id <- rownames(comm)
  colnames(res) <- c('pval','fc','otu.id')
  res <- res[res$otu.id %in% otu.name,]
  res$fdr = p.adjust(res$pval, 'fdr') 
  
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}

Omnibus.wrapper0 <- function(dat, cutoff = 0.05, FDR = T) { # grp.name can be categorical and numerical
  t1 = Sys.time()
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  covariate <- names(dat$meta.dat[1])
  # For match meta.dat and gmpr.size.factor
  counts <- dat$counts[,names(dat$gmpr.size.factor)]
  meta <- dat$meta[names(dat$gmpr.size.factor),,drop=F]
  Omnibus.obj <- ZISeq(counts,meta, size.factor = dat$gmpr.size.factor,grp.name = covariate, method = 'omnibus', filter = F, winsor =T)
  #Omnibus.obj <- ZISeq(dat$counts, dat$meta.dat, size.factor = dat$gmpr.size.factor,
  #                     grp.name = covariate, method = 'omnibus') # make sure the rownames(meta.dat) = colnames(otu.tab)
  #Omnibus.obj <- ZISeq(dat$counts, dat$meta.dat, size.factor =NULL, grp.name = covariate, method = 'omnibus',winsor = T,filter = F)
  res <- as.data.frame(Omnibus.obj$result) %>% 
    rownames_to_column('otu.id') %>%
    mutate(pval = p.value) %>% 
    dplyr::select(pval,fcs = abund.LFC.grpbefore.est, otu.id) %>% mutate(fdr = p.adjust(pval,'fdr'))   
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}
ANCOMBC.wrapper0 <- function(dat, cutoff = 0.05, FDR = T) {
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
    out = ancombc(phyloseq = phyloseq, formula = names(dat$meta.dat[1]), p_adj_method = "fdr", lib_cut = 0, zero_cut = 1,group = names(dat$meta.dat[1]), struc_zero =T, conserve = T)
    #out = ancombc(phyloseq = phyloseq, formula = names(dat$meta.dat)[1], p_adj_method = "fdr", lib_cut = 800, group = names(dat$meta.dat[1]))
  }else{
    out = ancombc(phyloseq = phyloseq, formula = names(dat$meta.dat)[1], p_adj_method = "fdr",lib_cut = 0, zero_cut = 1, group = NULL, struc_zero = FALSE, neg_lb = FALSE,conserve = T)
  }
  res <- as.data.frame(as.matrix(cbind(pval = out$res$p_val$grp,fdr = out$res$q_val$grp, fcs= out$res$W))) # Caution: the code adjusted since the cluster makes errors
  colnames(res) = c('pval','fdr','fcs')
  res$otu.id <- rownames(out$res$p_val)
  na <- sum(is.na(res$fdr))/nrow(res)
  sig <- sum(res$fdr[!is.na(res$fdr)] <= cutoff)
  time <- difftime(Sys.time(),t1, units = 'secs') %>% as.numeric()
  
  if(FDR){
    res <- res %>% inner_join(dat$truth, by = 'otu.id')
    tp <- sum(res$fdr <= cutoff & res$diff.otu.ind ==TRUE, na.rm = T)
    tn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    fn <- sum(res$fdr > cutoff &  res$diff.otu.ind ==TRUE, na.rm = T)
    fp <- sum(res$fdr <= cutoff &  res$diff.otu.ind ==FALSE, na.rm = T)
    tpr <- ifelse(!tp, 0, tp/(tp+fn))
    fpr <- ifelse(!fp, 0, fp/(tn+fp))
    fdr <- ifelse(!fp, 0, fp/(tp+fp))
    F1 <- ifelse(!tp, 0 ,(2*tp)/(2*tp+fp+fn))
    MCC <- ifelse(!((tp*tn)-(fp*fn)),0,((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return(list(tpr = tpr,fdr = fdr, fpr = fpr, fp = fp, F1 = F1, MCC =MCC,res = res, na = na, time = time, sig = sig))
  }else{
    return(list(res = res, na = na, time = time, sig = sig))
  }
}
