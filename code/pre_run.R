func <- function(part, paras) {
  set.seed(part)
  model.paras <- paras$model.paras
  nSams <- paras$nSams
  nOTUs <- paras$nOTUs
  diff.otu.pcts <- paras$diff.otu.pcts
  diff.otu.directs <- paras$diff.otu.directs
  diff.otu.modes <- paras$diff.otu.modes
  covariate.types <- paras$covariate.types
  covariate.eff.means <- paras$covariate.eff.means
  confounder.types <- paras$confounder.types
  depth.mus <- paras$depth.mus
  depth.conf.factors <- paras$depth.conf.factors
  models <- paras$models
  nSubs <- paras$nSubs
  include.top.otu <- paras$include.top.otu
  
  resdir <- paras$resdir
  prefix <- paras$prefix
  resdir <- gsub('/$', '', resdir)
  
  sum.list <- list()
  cat(date(), '\n')
  dir.create(paste0('../iter',part))
  
  for (diff.otu.pct in diff.otu.pcts) {
    for (diff.otu.mode in diff.otu.modes) {
      for (confounder.type in confounder.types) {
        for (depth.conf.factor in depth.conf.factors){
          for (covariate.type in covariate.types) {
            for (nSam in nSams) {
              for (nOTU in nOTUs) {
                for (depth.mu in depth.mus) {
                  for (model in models){
                    for (nSub in nSubs){
                      for (covariate.eff.mean in covariate.eff.means){
                        for (diff.otu.direct in diff.otu.directs){
                          Sim.obj <- SimulateSeq(otu.tab,
                                                 nOTU = nOTU, nSam = nSam, diff.otu.pct = diff.otu.pct, diff.otu.direct = diff.otu.direct, diff.otu.mode = diff.otu.mode,
                                                 covariate.type = covariate.type, covariate.eff.mean = covariate.eff.mean, confounder.type = confounder.type, depth.mu =depth.mu, depth.conf.factor = depth.conf.factor,
                                                 model = model, nSub = nSub,include.top.otu = include.top.otu)
                          cat('1.Simulation finshed! \n')
                          #-- original simulated otu.tab
                          otu.tab.sim <- Sim.obj$otu.tab.sim;dim(otu.tab.sim)
                          covariates <- Sim.obj$X
                          truth <- as.data.frame(cbind(Sim.obj$diff.otu.ind, Sim.obj$otu.names)) %>% dplyr::rename(diff.otu.ind = V1, otu.id = V2)
                          meta.dat <- as.data.frame(Sim.obj$X) %>% dplyr::rename(grp = V1)
                          rownames(meta.dat) <- colnames(otu.tab.sim)
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
                          
                          gmpr.covariates <- covariates[!is.na(size.factor),]
                          gmpr.truth <- truth[truth$otu.id %in% rownames(gmpr.counts),]
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
                            Wrench.truth <- truth[truth$otu.id %in% rownames(Wrench.nf),]
                            Wrench.meta.dat <- meta.dat[rownames(meta.dat) %in% colnames(Wrench.nf),, drop=FALSE]
                            cat('3.Wrench finshed! \n')
                          }
                          
                          ## raw otu table prevelance filtration
                          prop <- t(t(otu.tab.sim) / colSums(otu.tab.sim))
                          prop <- prop[rowSums(prop!=0) > 0.1 * ncol(prop), , drop=FALSE]
                          otu.tab.sim <- otu.tab.sim[rownames(prop), , drop=FALSE]
                          ## raw otu tableminimal abundance filtration
                          prop <- prop[rowMaxs(prop) > 0.002, , drop=FALSE]
                          otu.tab.sim <- otu.tab.sim[rownames(prop), , drop=FALSE] # used for analysis
                          cat('4.Preprocessing finshed! \n')
                          
                          
                          dat <- list(counts = otu.tab.sim, prop = prop, truth = truth, covariates = covariates, meta.dat = meta.dat,
                                      Wrench.nf = Wrench.nf, compositionalFactors= compositionalFactors, normalizationFactors = normalizationFactors, Wrench.covariates = Wrench.covariates,Wrench.truth= Wrench.truth,Wrench.meta.dat = Wrench.meta.dat,
                                      gmpr.comm=gmpr.comm, gmpr.counts = gmpr.counts, gmpr.size.factor = gmpr.size.factor, gmpr.covariates=gmpr.covariates, gmpr.truth = gmpr.truth, gmpr.meta.dat = gmpr.meta.dat)
                          save(dat,file=paste0('../iter',part,'/',depth.mu, diff.otu.mode, model, nSub, nOTU, covariate.type, depth.conf.factor, diff.otu.direct, confounder.type, nSam, covariate.eff.mean, diff.otu.pct,'.Rdata'))
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}




clsapply <- function(dat, func, df=NULL, tempdir="~/project/tempZ",
                     queque="7-day", timespan="540", mem="8G", req.pack='stats',
                     tidy=F) {
  # Apply to the row of dat by default (dim=2 for column)
  # Create temp dir
  if (!file.exists(tempdir)) {
    dir.create(tempdir)
  }
  setwd(tempdir)
  
  save(dat, func, df, otu.tab, dirmult.paras, file=("dat.RData"))
  # methods <- df$methods
  nJob <- length(dat)
  # Create R script
  txt <- paste('
               args=(commandArgs(TRUE))
               if(length(args)==0){
               print("No arguments supplied.")
               ##supply default values
               } else{
               part = args[1]
               part = as.integer(part)
                              }',
               paste(paste0('\nrequire(',req.pack, ')'), collapse="\n"),
               '\n date()
      load("dat.RData") 
      
                        if (is.list(dat)) {
                        item <- dat[[part]]
                        } 
                        if (is.vector(dat)) { 
                        item <- dat[part]
                        } 

                        source("../DM1.R")
                        source("../GMPR.R")
      pkg <- c("phyloseq","aod","reshape","MASS","corncob","readr","DESeq2", "ALDEx2", "metagenomeSeq", "edgeR", "GUniFrac", "grDevices", "dirmult", "exactRankTests","nlme", "dplyr", "magrittr", "tidyr", "protoclust", "ggplot2", "compositions","rmutil","tibble","reticulate","dacomp","LDM","Wrench","RioNorm2")
      lapply(pkg, require, character.only = TRUE)

                        res0 <- func(item, df)
                        date()
                        ')
  writeLines(txt, "dat_script.R")
  
  # Submit job
  cat("Submit jobs ...")
  prefix <- paste(c("J", sample(LETTERS, 4, repl=TRUE)), collapse="")
  for (part in 1:nJob) {
    rfile <- "dat_script.R"
    rout <- paste0(part, ".Rout")
    sh <- paste(
      "qsub",
      paste0("-N ", prefix, part),
      "-j y",
      "-cwd",
      "-q", queque,
      "-m abe",
      paste0("-l h_vmem=", mem),
      "-V",
      paste("-o ",  part, ".out", sep=""),
      paste("-e ",  part, ".err", sep=""),
      "-b y",
      paste("\'R CMD BATCH --no-save --no-restore \"--args ", part, "\" ",
            rfile, " ", rout, "\'", sep="")
    )
    print(sh)
    system(sh)
  }
}


