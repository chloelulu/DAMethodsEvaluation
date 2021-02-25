diagnose <- function(Y, Z, cutoff = 0.05){
  options(warn=-1)
  # Y: otu table with row otu, column sample
  # Z : covariates vector
  YY <- Y + 0.5
  prop <- t(t(YY)/colSums(YY))
  
  top5 <- tail(names(sort(rowMeans(prop)))) # top5 otu
  abundant <- tail(names(sort(rowMeans(prop)))) # top 25% abundant otu
  rare <- tail(names(sort(rowMeans(prop))))

  clrY <- t(clr(prop))
  ZZ <- scale(Z)
  
  fit <- lm(clrY ~ ZZ)
  
  est <- coef(fit)[2,]
  summary <- summary(fit)
  
  p <- NULL
  for(i in 1:length(summary)){
    p <- c(p,coef(summary[[i]])[2,4])
  }
  
  names(p) <- names(est)
  q <- p.adjust(p, method = 'fdr')
  idx <- which(q <= cutoff)
  
  pos.ind <- which(q<=cutoff & est>0)
  neg.ind <- which(q<=cutoff & est<0)
  
  pos.est <- est[pos.ind]
  neg.est <- est[neg.ind]
  
  ## check depth confounding
  dep <- colSums(Y)
  diff.seq.p <- cor.test(dep, as.numeric(Z), alternative = "two.sided", method = 'spearman')$p.value
  if(diff.seq.p <= 0.05){
    cat("Below inncluding rough diagnoses and recommendations for differential analysis method chosen: \n(1) Signficant sequencing depth confounding found! Please do rarefraction at first! \n")
  }else{
    cat("Below including diagnoses and recommendations for differential analysis method chosen: \n(1) No sequencing depth confounding found!\n")
  }

  ## calculate compositionality
  if(length(idx) > 0){
    prop0 <- rowMeans(prop[idx,,drop =F])
    compostionality <- abs(sum(exp(pos.est) * prop0[names(pos.ind)]) - sum(exp(abs(neg.est)) * prop0[names(neg.ind)]))
    
    ## check compisitonal include top 5 abundant taxa
    if(compostionality > 0 &  sum(names(idx) %in% top5) > 0){
      cat('(2) Your data is compositional! Please choose method under setting: compositional + top5 abundant taxa! \n')
    }else{
      cat('(2) Your data is compositional! Please choose method under setting: compositional! \n')
    }
  }else{
    compostionality <- 0
    cat('(2) Your data is not compositional! Please choose method under settings without "compostional" \n')
  }
  
  ## check signal density(how many taxa differentiate)
  if(length(idx) <= nrow(prop)*0.05){ cat('(3) Please chooes method under setting: singal density = Low ! \n')}
  if(length(idx) <= nrow(prop)*0.1 & length(idx) > nrow(prop)*0.05){ cat('(3) Please choose method under setting: singal density = Medium! \n') }
  if(length(idx) > nrow(prop)*0.1){ cat('(3) Please chooes setting: High singal density! \n')}
  
  ## check differential taxa mode
  if(sum(names(idx) %in% rare) > 0.8 *length(idx)){
    cat('(4) Please choose setting: differential taxa are rare taxa!\n')
  }else{
    cat('(4) Please choose setting: differential taxa are abundant taxa !\n')
    }
}


library(compositions)
require(GUniFrac)
data(throat.otu.tab)
data(throat.meta)
Y <- t(throat.otu.tab)
Z <- throat.meta$Age

diagnose(Y, Z)




