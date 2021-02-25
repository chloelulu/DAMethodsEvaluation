edgeR.gmpr<- function(dat) {
  if (is.null(rownames(dat$gmpr.counts))) {
    rownames(dat$gmpr.counts) <- paste0('O', 1:nrow(dat$gmpr.counts))
    otu.name <- rownames(dat$gmpr.counts) 
  } else {
    otu.name <- rownames(dat$gmpr.counts) 
  }
  d <- DGEList(counts=dat$gmpr.counts, group=dat$gmpr.covariates$grp)
  d <- edgeR::calcNormFactors(d)
  d <- estimateCommonDisp(d)
  d <- estimateTagwiseDisp(d)
  fit <- exactTest(d)
  res <- topTags(fit, n=nrow(dat$gmpr.counts))
  res <- as.data.frame(res)
  res <- res[otu.name,]
  res <- cbind(res, pval = res[, 'PValue'], fdr = res[, 'FDR'])
  res <- as.data.frame(res[,c('pval','fdr')]) %>% rownames_to_column('otu.id')
  res[is.na(res$fdr),'fdr'] = 1
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}

edgeR <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  d=DGEList(counts=dat$counts, group=dat$covariates)
  d=edgeR::calcNormFactors(d)
  d=estimateCommonDisp(d)
  d=estimateTagwiseDisp(d)
  fit=exactTest(d)
  res=topTags(fit, n=nrow(dat$counts))
  res=as.data.frame(res)
  res <- res[otu.name,]
  res <- cbind(res, pval = res[, 'PValue'], fdr = res[, 'FDR'])
  res <- as.data.frame(res[,c('pval','fdr')]) %>% rownames_to_column('otu.id')   
  res[is.na(res$fdr),'fdr'] = 1
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
  
}


edgeR.Wrench <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  d=edgeR::DGEList( counts=dat$counts,
                    group = as.matrix(dat$covariates),
                    norm.factors=dat$compositionalFactors)
  d=edgeR::calcNormFactors(d)
  d=estimateCommonDisp(d)
  d=estimateTagwiseDisp(d)
  fit=exactTest(d)
  res=topTags(fit, n=nrow(dat$counts))
  res=as.data.frame(res)
  res <- res[otu.name,]
  res <- cbind(res, pval = res[, 'PValue'], fdr = res[, 'FDR'])
  res <- as.data.frame(res[,c('pval','fdr')]) %>% rownames_to_column('otu.id')
  res[is.na(res$fdr),'fdr'] = 1
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}


metagenomeSeq2.gmpr <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  mgs = newMRexperiment(counts=dat$gmpr.counts)
  mgs = metagenomeSeq::cumNorm(mgs)
  normFactors(mgs) <- dat$gmpr.size.factor
  mod = model.matrix(~dat$gmpr.covariates$grp)
  settings = zigControl(verbose = F)
  fit = fitFeatureModel(mgs,mod)
  res = MRfulltable(fit,number = nrow(assayData(mgs)$counts))
  res <- res[otu.name, ]
  pval = res$pvalues
  # pval[is.na(pval)] = 1 # change NA into 1 ~ non-significant
  fdr = res$adjPvalues # MRfulltable produce fdr adjusted pvalues
  fdr[is.na(fdr)] = 1
  res <- cbind(res, pval = pval, fdr = fdr)
  res <- as.data.frame(res[,c('pval','fdr')]) %>% rownames_to_column('otu.id')
  res[is.na(res$fdr),'fdr'] = 1
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}


metagenomeSeq2.Wrench <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  mgs = newMRexperiment(counts=dat$counts)
  mgs = metagenomeSeq::cumNorm(mgs)
  normFactors(mgs) <- dat$normalizationFactors
  mod = model.matrix(~dat$covariates)
  settings = zigControl(verbose = F)
  fit = fitFeatureModel(mgs,mod)
  res = MRfulltable(fit,number = nrow(assayData(mgs)$counts))
  res <- res[otu.name, ]
  pval = res$pvalues
  # pval[is.na(pval)] = 1 # change NA into 1 ~ non-significant
  fdr = res$adjPvalues # MRfulltable produce fdr adjusted pvalues
  fdr[is.na(fdr)] = 1
  res <- cbind(res, pval = pval, fdr = fdr)
  res <- as.data.frame(res[,c('pval','fdr')]) %>% rownames_to_column('otu.id')
  res[is.na(res$fdr),'fdr'] = 1
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}



metagenomeSeq2 <- function(dat) {
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
  fit = fitFeatureModel(mgs,mod) # MSeq recommned this over the MSeq1, which uses fitZig()
  res = MRfulltable(fit,number = nrow(assayData(mgs)$counts))
  res <- res[otu.name, ]
  pval = res$pvalues
  # pval[is.na(pval)] = 1 # change NA into 1 ~ non-significant
  fdr = res$adjPvalues # MRfulltable produce fdr adjusted pvalues(adjustMethod = "fdr")
  fdr[is.na(fdr)] = 1
  res <- cbind(res, pval = pval, fdr = fdr)
  res <- as.data.frame(res[,c('pval','fdr')]) %>% rownames_to_column('otu.id')
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}




wilcox <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  prop <- t(t(dat$counts) / colSums(dat$counts))
  covariate <- dat$covariates
  pval <- apply(prop, 1,  function (y) wilcox.test(y ~ covariate)$p.value)
  pval <- pval[otu.name]
  # pval[is.na(pval)] = 1 # change NA into 1 ~ non-significant
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  fdr[is.na(fdr)]=1
  res <- cbind(pval = pval, fdr = fdr) %>% as.data.frame() %>% rownames_to_column('otu.id')
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}


ZicoSeq <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  grp.name <- names(dat$meta.dat)[1]
  comm = as.data.frame(dat$counts)
  meta.dat = as.data.frame(dat$meta.dat)
  zicoseq.obj <- ZicoSeq::ZicoSeq(meta.dat = meta.dat, comm = comm, grp.name = grp.name, adj.name = NULL,
                         is.winsor = T, winsor.qt = 0.97,
                         # Posterior sampling
                         is.prior = T, 
                         prior.dist = c('BetaMix'), post.method = c('sample'), post.sample.no = 25, 
                         # Link functions
                         link.func = list(function (x) x^0.25, function (x) x^0.5, function (x) x^0.75), stats.combine.func = max,
                         # Permutation
                         perm.no = 99,  strata = NULL, 
                         # Multiple stage normalization
                         stage.no = 6, topK = NULL, stage.fdr = 0.75, stage.max.pct = 0.50,  
                         # Tree-based FDR control and family-wise error rate control
                         is.fwer = FALSE, is.tree.fdr = FALSE, tree = NULL, 
                         verbose = TRUE, return.comm = FALSE, return.perm.F = FALSE)
  pval <- zicoseq.obj$p.raw
  # pval[is.na(pval)] = 1
  fdr <- zicoseq.obj$p.adj.fdr
  fdr[is.na(fdr)] = 1 # change NA into 1 ~ non-significant
  res <- as.data.frame(cbind(pval = pval,fdr = fdr)) %>% rownames_to_column('otu.id')
  res[is.na(res$fdr),'fdr'] = 1
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}


Rarefy <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  comm.r <- t(GUniFrac::Rarefy(t(dat$counts))$otu.tab.rff)
  pval <- apply(comm.r, 1,  function (y) wilcox.test(y ~ dat$meta.dat$grp)$p.value)
  # pval[is.na(pval)] = 1 # change NA into 1 ~ non-significant
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  fdr[is.na(fdr)] = 1
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id')
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
  
}

DESeq <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  
  design <- as.formula(paste('~', names(dat$meta.dat[1])))
  dds <- DESeqDataSetFromMatrix(countData = dat$counts, colData = dat$meta.dat, design= design)
  dds <- DESeq2::DESeq(dds)
  res <- as.data.frame(as.matrix(results(dds)))
  pval <- res$pvalue
  names(pval) <- rownames(res)
  # pval[is.na(pval)] = 1
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  fdr[is.na(fdr)]=1
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id')
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
  
}

DESeq.gmpr <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  
  design <- as.formula(paste('~', names(dat$gmpr.covariates[1])))
  dds <- DESeqDataSetFromMatrix(countData = dat$gmpr.counts, colData = dat$gmpr.covariates, design= design)
  sizeFactors(dds) <- dat$gmpr.size.factor
  
  dds <-  DESeq2::DESeq(dds)
  res <- as.data.frame(as.matrix(results(dds)))
  pval <- res$pvalue
  names(pval) <- rownames(res)
  # pval[is.na(pval)] = 1
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  fdr[is.na(fdr)]=1
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id')
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}


DESeq.Wrench <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  design <- as.formula(paste('~', names(dat$meta.dat[1])))
  dds <- DESeqDataSetFromMatrix(countData =dat$counts, colData = dat$meta.dat, design= design)
  DESeq2::sizeFactors(dds) <- dat$normalizationFactors
  dds <-  DESeq2::DESeq(dds)
  res <- as.data.frame(as.matrix(results(dds)))
  pval <- res$pvalue
  names(pval) <- rownames(res)
  # pval[is.na(pval)] = 1
  fdr = p.adjust(pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  fdr[is.na(fdr)]=1
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id')
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
  
}




Aldex2 <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  comm <- as.data.frame(dat$counts)
  covariate <- as.character(dat$covariates)
  #aldex.obj <- aldex(comm, covariate)
  aldex.obj <- aldex(comm, covariate, test="t", effect=FALSE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
  pval <- aldex.obj$wi.ep
  names(pval) <- rownames(aldex.obj)
  # pval[is.na(pval)] = 1
  fdr = p.adjust(pval, 'fdr')
  fdr[is.na(fdr)]=1
  res <- as.data.frame(cbind(pval = pval, fdr = fdr)) %>% rownames_to_column('otu.id')
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}



Omnibus <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  covariate <- names(dat$meta.dat[1])
  Omnibus.obj <- ZISeq(dat$counts, dat$meta.dat, size.factor = NULL, winsor = TRUE, winsor.qt = 0.97,
                       grp.name = covariate, adj.name = NULL, method = 'omnibus',filter = F) # make sure the rownames(meta.dat) = colnames(otu.tab)
  
  res <- as.data.frame(Omnibus.obj$result) %>% 
    rownames_to_column('otu.id') %>%
    mutate(pval = p.value) %>% 
    dplyr::select(pval,otu.id)
  # res$pval[is.na(res$pval)] <- 1
  res <- res %>%
    mutate(fdr = p.adjust(pval,'fdr'))
  res[is.na(res$fdr),'fdr'] = 1
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}


ANCOM2 <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  covariate <- names(dat$meta.dat[1])
  meta_data <- dat$meta.dat %>% mutate(Sample.ID = rownames(.))
  feature_table = dat$counts; sample_var = "Sample.ID"; group_var = covariate
  # out_cut = 0; zero_cut = 0;
  prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, neg_lb = FALSE)
  
  feature_table <- prepro$feature_table # Preprocessed feature table
  meta_data <- prepro$meta_data # Preprocessed metadata
  struc_zero <- prepro$structure_zeros # Structural zero info
  
  ANCOM2.obj <- ANCOM(feature_table, meta_data, struc_zero, main_var = covariate, p_adj_method = "fdr",
                      alpha = 0.05, adj_formula = NULL, rand_formula = NULL)
  res <- as.data.frame(ANCOM2.obj$out) %>% mutate(otu.id = taxa_id) %>% dplyr::select(-taxa_id) %>% dplyr::rename(taxa = otu.id)
  
  sig <- sum(res$detected_0.7)
  return(list(sig = sig,res = res))
}


RAIDA <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts) 
  } else {
    otu.name <- rownames(dat$counts) 
  }
  df.raida <- as.data.frame(dat$meta.dat) %>% rownames_to_column('otu.id') %>% 
    inner_join(as.data.frame(t(dat$counts))  %>% rownames_to_column('otu.id'))
  df.raida <- df.raida[order(-df.raida$grp),] 
  rownames(df.raida) <- NULL
  n.lib <- c(as.numeric(table(df.raida$grp)[1]),as.numeric(table(df.raida$grp)[2]))
  df.raida <- df.raida %>% dplyr::select(-'grp') %>% column_to_rownames('otu.id') %>% t(.) %>% as.data.frame(.)
  raida.obj <- raida(df.raida, n.lib, show.ref.features=T)
  
  res <- raida.obj$result %>% 
    rownames_to_column('otu.id') %>% 
    mutate(pval = p) %>% 
    dplyr::select(otu.id, pval) %>% 
    mutate(fdr = p.adjust(pval, 'fdr')) %>% mutate(ref.raida = 0)
  res[is.na(res$fdr),'fdr'] = 1
  ref = as.data.frame(raida.obj$reference.features) %>% mutate(pval = 1, fdr = 1) %>% mutate(ref.raida = 1)
  colnames(ref)[1] = 'otu.id'
  res = rbind(res,ref)
  res = res %>% arrange(factor(otu.id, levels =rownames(df.raida)[rownames(df.raida) %in% res$otu.id]))
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}



# https://github.com/barakbri/dacomp
dacomp <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  dacomp.otu <- t(dat$counts)
  dacomp.meta <- as.vector(dat$covariates)
  references = dacomp.select_references(X = dacomp.otu, median_SD_threshold = 0.6, verbose = T)
  dacomp.obj <- dacomp.test(X = dacomp.otu, y = dacomp.meta, ind_reference_taxa = references,
                            nr_perm = 1000, test = DACOMP.TEST.NAME.WILCOXON, verbose = T,q = 0.05)
  pval <- dacomp.obj$p.values.test
  names(pval) <- otu.name
  
  res <- as.data.frame(pval) %>%
    rownames_to_column('otu.id') %>%
    dplyr::select(otu.id, pval) %>%
    mutate(fdr = p.adjust(pval, 'fdr')) %>% 
    mutate(ref.dacomp = 0)
  res[res$otu.id %in% colnames(dacomp.otu)[references$selected_references], 'ref.dacomp'] = 1
  res[is.na(res$fdr) & res$ref.dacomp ==1,'fdr'] = 1
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig, res = res))
}



ldm <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  ldm.otu <<- t(dat$counts) %>% as.data.frame()
  ldm.obj <- LDM::ldm(formula=ldm.otu ~ grp, data = dat$meta.dat)
  res <- as.data.frame(t(ldm.obj$p.otu.omni)) %>%
    rownames_to_column('otu.id') %>%
    mutate(pval = V1) %>%
    dplyr::select(otu.id, pval) %>% mutate(fdr = p.adjust(pval, method = 'fdr'))
  res[is.na(res$fdr),'fdr'] = 1
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}


RioNorm <- function(dat) {
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
  combined_res$padj = as.numeric(as.character(combined_res$padj)) 
  combined_res$pvalue = as.numeric(as.character(combined_res$pvalue)) 
  
  res <- combined_res %>% as.data.frame() %>% dplyr::rename(otu.id = id, fdr = padj) 
  res[is.na(res$fdr),'fdr'] = 1
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}



glmquassi <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  comm <- dat$counts
  covariate <- dat$covariates
  
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
  # res[is.na(res$pval),'pval'] = 1
  res$fdr = p.adjust(res$pval, 'fdr') # MRfulltable produce fdr adjusted pvalues
  res[is.na(res$fdr),'fdr'] = 1
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
  }


BBinomial <- function(dat) {
  if (is.null(rownames(dat$counts))) {
    rownames(dat$counts) <- paste0('O', 1:nrow(dat$counts))
    otu.name <- rownames(dat$counts)
  } else {
    otu.name <- rownames(dat$counts)
  }
  
  sam <- sample_data(dat$meta.dat)
  otu.tab <- otu_table(dat$counts, taxa_are_rows = T)
  phyloseq <- merge_phyloseq(otu.tab, sam)
  corncob <- differentialTest(formula = ~ grp,
                              phi.formula = ~ grp,
                              formula_null = ~ 1,
                              phi.formula_null = ~ grp,
                              test = "Wald", boot = FALSE,
                              data = phyloseq,
                              fdr_cutoff = 0.05)
  res <- as.data.frame(cbind(pval = corncob$p, fdr = corncob$p_fdr)) %>% rownames_to_column('otu.id')
  res <- res[res$otu.id %in% otu.name,]
  res[is.na(res$fdr),'fdr'] = 1
  sig <- sum(res$fdr <= 0.05)
  return(list(sig = sig,res = res))
}




