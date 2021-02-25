################################################################################################
######################  random sample n sample and compare with the rest  ########################
################################################################################################
pkg = c('dplyr','GUniFrac','vegan','ape','microbiome','tibble','reshape2','RColorBrewer','plyr','ggpubr','MASS','dirmult')
sapply(pkg, require, character = TRUE)
setwd('~/Documents/Mayo_Research/Semiparametric/ReArrange1204/')
source('code/Func.R')
QC <- function(otu.tab, prev = 0.1, minp = 0.002){
  prop = t(t(otu.tab)/colSums(otu.tab))
  prop <- prop[rowSums(prop!=0) > prev * ncol(prop), , drop=FALSE]
  otu.tab <- otu.tab[rownames(prop), , drop=FALSE]
  
  prop <- prop[rowMaxs(prop) > minp, , drop=FALSE]
  otu.tab <- otu.tab[rownames(prop), , drop=FALSE]
  
  idx.f <- which(colSums(otu.tab)>1000)
  # idx.f <- apply(otu.tab, 1, function(x) sum(x > 2)) > round(ncol(otu.tab) * 0.05)
  otu.tab <- otu.tab[,idx.f]
  return(otu.tab = otu.tab)
}


# load('~/Documents/Mayo_Research/DataLib/VaginalIntroitus_V35.RData') # OTU dataset
# load('~/Documents/Mayo_Research/DataLib/VaginalIntroitus_V35_dirmult.RData') # For saving computational time, already prepared dirmult.paras
# name = 'Vaginal'
# otu.tab = VaginalIntroitus_V35

load('~/Documents/Mayo_Research/DataLib/Stool_V35.RData') # OTU dataset
load('~/Documents/Mayo_Research/DataLib/Stool_V35_dirmult.RData') # For saving computational time, already prepared dirmult.paras
name = 'Stool'
otu.tab = Stool_V35
otu.tab <- QC(otu.tab = otu.tab)

load('~/Documents/Mayo_Research/DataLib/UrogenitalTract_V35_QC.RData') # OTU dataset
load('~/Documents/Mayo_Research/DataLib/UrogenitalTract_V35_dirmult.RData') # For saving computational time, already prepared dirmult.paras
name = 'UrogenitalTract'


idx = colnames(otu.tab)
half1.idx = sample(idx, 0.5 * length(idx))
half2.idx = idx[!(idx %in% half1.idx)]
otu.tab1 = otu.tab[,half1.idx]
otu.nonsample = otu.tab[,half2.idx]
dim(otu.tab1);dim(otu.nonsample)
otu.tab = otu.tab1


nOTU = 'nOTU_L5'
nSam = 'nSam_L1'
diff.otu.pct= 'low'
diff.otu.mode = 'abundant'
diff.otu.direct = 'balanced'
covariate.type = 'binary'
covariate.eff.mean = 'none'
confounder.type = 'none'
depth.mu = 'D3'
depth.conf.factor = 'none'
include.top.otu = FALSE
model = 'loglinear'
nSub = 'Sub_L1'
model.paras = NULL
par(mfrow = c(2,2))
Sim.obj <- SimulateSeq(otu.tab = otu.tab,
                       nOTU = nOTU, diff.otu.pct = diff.otu.pct, diff.otu.direct = diff.otu.direct, diff.otu.mode = diff.otu.mode,
                       covariate.type = covariate.type, covariate.eff.mean = covariate.eff.mean, confounder.type = confounder.type, depth.mu =depth.mu, depth.conf.factor = depth.conf.factor,
                       model = model, nSub = nSub,include.top.otu = include.top.otu)
otu.tab.dir <- Sim.obj$otu.tab.sim
dim(otu.tab.dir);dim(otu.nonsample)


# start compare
rary.test = Rarefy(t(otu.nonsample), depth = 1000)$otu.tab.rff %>% t()
rary.dir = Rarefy(t(otu.tab.dir), depth = 1000)$otu.tab.rff %>% t()
dim(rary.test);dim(rary.dir)

rary.test.r = t(t(rary.test)/colSums(rary.test))
rary.dir.r = t(t(rary.dir)/colSums(rary.dir))

# OTU prevalence
pre.test = as.data.frame(prevalence(rary.test, count=FALSE)) %>% rownames_to_column('OTU')
colnames(pre.test)[2] = 'raw'
pre.dir = as.data.frame(prevalence(rary.dir, count=FALSE)) %>% rownames_to_column('OTU')
colnames(pre.dir)[2] = 'dirchlet'
dim(pre.test);dim(pre.dir)



qq = inner_join(pre.dir, pre.test)
head(qq)
cols = c(Semiparametric = '#ff8849', `Non-trained half`= '#69be28')

head(qq)
qqplot3 <- function(qq, title = title, xdensity = '',digits = 2, sqrt= F){
  colnames(qq)[2:3] = c('Non-trained half','Semiparametric')
  fmt_dcimals <- function(decimals=0){
    # return a function responpsible for formatting the 
    # axis labels with a given number of decimals 
    function(x) as.character(round(x,decimals))
  }
  
  qq_24 = as.data.frame(stats::qqplot(qq[,'Non-trained half'], qq[,'Semiparametric'], plot.it=FALSE)) %>% dplyr::rename(Semiparametric = y, `Non-trained half` = x)
  qq24 = ggplot(qq_24) + 
    geom_abline(intercept = 0, color = 'grey10', cex = 1.5, alpha = 0.7) +
    geom_point(aes(x=`Non-trained half`, y= Semiparametric), color = '#ff8849') +
    theme_bw()+
    scale_x_continuous(labels =fmt_dcimals(digits)) +
    scale_y_continuous(labels =fmt_dcimals(digits)) +
    theme(axis.text.x = element_text(color="black", size = 18),
          plot.margin = unit(c(1,1,1,1), "cm"),
          axis.text.y = element_text(color="black", size = 18),
          axis.title = element_text(color="black", size = 18),
          legend.position = 'none')+
    labs(x='Non-Trained Half',y='Semiparametric')

  prev = qq_24 %>% melt()
  if(sqrt){
    prev = prev %>% mutate(value = sqrt(value))
  }
  
  mu <- ddply(prev, "variable", summarise, grp.mean=mean(value))
  p2 = ggplot(prev, aes(x=value, fill=variable)) +
    geom_density(alpha=0.5)+
    scale_fill_manual(values = cols) +
    theme_minimal() + theme_classic() +
    theme(axis.text = element_text(color="black", size = 18,vjust = 0.5),
          legend.text= element_text(color="black", size = 16),
          plot.margin = unit(c(1,1,1,1), "cm"),
          
          legend.title = element_text(color="black", size = 18),
          axis.title = element_text(color="black", size = 18),
          legend.position = 'right')+
    labs(x = xdensity, y = 'Density', fill = '')+
    scale_y_continuous(trans = sqrt_trans(),
                       breaks = trans_breaks("sqrt", function(x) x^2),
                       labels = trans_format("sqrt", math_format(.x^2)))
  # prev1230 = ggarrange(p1, p2, nrow = 2)
  plt = ggarrange(p2, qq24, nrow = 1, common.legend = TRUE)
  return(plt)
}
qqplot3(qq, title = 'OTU prevelance',xdensity = 'OTU prevelance')
ggsave(paste0('result/overfitting/',name,'_OTU prevelance_OverfittingTest.pdf'), plot = last_plot(), width = 8, height = 4, dpi = 60, units = 'in')


#  mean abundance
mn2 = as.data.frame(rowMeans(rary.test.r)) %>% rownames_to_column('OTU')
colnames(mn2)[2] ='raw'
mn4 = as.data.frame(rowMeans(rary.dir.r)) %>% rownames_to_column('OTU')
colnames(mn4)[2] = 'dirchlet'
mn_abund = inner_join(mn2, mn4)
qqplot3(mn_abund, title = 'OTU mean abundance',xdensity = 'sqrt(Abundance)', sqrt = T)
ggsave(paste0('result/overfitting/',name,'OTU mean abundance_OverfittingTest.pdf'), plot = last_plot(), width = 8, height = 4, dpi = 60, units = 'in')


# OTU abundance variation (sd / mu or arcsin(sqrt))
tran2 = asin(sqrt(rary.test.r))
tran4 = asin(sqrt(rary.dir.r))
var2 = apply(tran2, 1, function(x) var(x))
var4 = apply(tran4, 1, function(x) var(x))
var = as.data.frame(var2) %>% rownames_to_column('OTU') %>%
  inner_join(as.data.frame(var4) %>% rownames_to_column('OTU')) %>% 
  dplyr::rename(raw = var2,dirchlet = var4) #BBmix = var1,
qqplot3(var, title = 'OTU abundance variation',xdensity = 'sqrt(Variance of abundance)', digits = 3, sqrt = T)
ggsave(paste0('result/overfitting/',name,'OTU abundance variation_OverfittingTest.pdf'), plot = last_plot(), width = 8, height = 4, dpi = 60, units = 'in')

# sample distribution
div2 = vegan::diversity(rary.test, index = "shannon", MARGIN = 1, base = exp(1)) %>% as.matrix()%>% as.data.frame() %>% rownames_to_column('sampleid')%>% dplyr::rename(raw = V1)
div4 = vegan::diversity(rary.dir, index = "shannon", MARGIN = 1, base = exp(1)) %>% as.matrix()%>% as.data.frame() %>% rownames_to_column('sampleid')%>% dplyr::rename(dirchlet = V1)

div0 = inner_join(div2, div4);head(div0)
qqplot3(div0, title = 'Shannon Diversity',xdensity = 'Shannon Diversity')
ggsave(paste0('result/overfitting/',name,'alphadiversity_OverfittingTest.pdf'), plot = last_plot(), width = 8, height = 4, dpi = 60, units = 'in')


# PCoA compare sample dist
a2 = otu.nonsample
colnames(a2) = paste0('raw.', colnames(a2))
a4 = otu.tab.dir
colnames(a4) = paste0('Semiparametric.', colnames(a4))
a4 = a4[rownames(a2),] # make sure the OTU.id are the same

a = cbind(a2,a4)
rary = Rarefy(t(a), depth = 1000)$otu.tab.rff %>% t()
m2 = colnames(a2) %>% as.data.frame() %>% mutate(grp ='Non-Trained half')
m4 = colnames(a4) %>% as.data.frame() %>% mutate(grp ='Semiparametric')
m = rbind(m2, m4)
colnames(m)[1] = 'SampleID'
metadata <- sample_data(m %>% mutate(sampleid = SampleID) %>% column_to_rownames('sampleid'))
otutable <- otu_table(rary, taxa_are_rows = T)
rary <- merge_phyloseq(otutable,metadata) 
pslog <- transform_sample_counts(rary, function(x) log(1 + x))
out.pcoa.log <- ordinate(pslog,  method = "PCoA", distance = "bray")
evals <- out.pcoa.log$values[,1]
eng <- as.data.frame(out.pcoa.log$vectors)[,c(1:2)] %>% rownames_to_column('SampleID') %>% inner_join(as.data.frame(as.matrix(sample_data(rary))) )
pc1 <- round((out.pcoa.log$values$Eigenvalues[1])/sum(out.pcoa.log$values$Eigenvalues),2)
pc2 <- round((out.pcoa.log$values$Eigenvalues[2])/sum(out.pcoa.log$values$Eigenvalues),2)
col = c('Non-Trained half' = brewer.pal(8,'Dark2')[1], 'Semiparametric' = brewer.pal(8,'Dark2')[3] )

pcoa = ggplot(eng,aes(x = Axis.1, y = Axis.2,shape = grp, fill = grp)) +
  stat_ellipse(geom = "polygon", alpha = 0.3,type = "norm")+ 
  geom_point(size = 2)+
  scale_color_manual(values = c('black','black','black','black'))+
  scale_shape_manual(values = c(24, 23,21,11))+
  scale_fill_manual(values = col)+
  labs(fill = 'Dataset', shape = 'Dataset') +
  xlab(paste0('PC1, ',(paste0(pc1*100,'%')))) + ylab(paste0('PC2, ',paste0(pc2*100, '%'))) + 
  theme_bw() +
  theme(text = element_text(size = 16, color = "black"))
pcoa

ggsave(paste0('result/overfitting/',name,'PCoA_OverfittingTest.pdf'), plot = last_plot(), width = 8, height = 6, dpi = 60, units = 'in')


h2 = rary.test.r; h4 = rary.dir.r
dim(h2);dim(h4)
r2 = sqrt(h2);r4 = sqrt(h4)

TOPn = T
if(TOPn){
  top = names(sort(rowMeans(r2), decreasing = T))[1:50]
  r2 = r2[top,];r4 = r4[top,]
}

cbind = cbind(r2, r4)
col.scheme <- brewer.pal(9, 'YlOrRd')
max(r2);max(r4)
breaks <- c(0, seq(min(cbind[cbind != 0]), 0.2, len=6),seq(0.3,max(cbind),length = 3))

save_pheatmap_png <- function(x, filename, width=1800, height=1700, res = 200) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

heat2 = pheatmap::pheatmap(r2,legend = TRUE, cluster_rows = F, 
                           cluster_cols = F,clustering_method = 'ward.D2',color = col.scheme,
                           breaks = breaks,
                           annotation_legend = F,show_colnames = F,show_rownames = F,fontsize_row =6,fontsize_col =6,
                           main = paste0('Non-trained half'))
save_pheatmap_png(heat2, paste0('result/overfitting/',name,'_heat_Nontrained_half.png'))

heat4 = pheatmap::pheatmap(r4,legend = TRUE, cluster_rows = F, 
                           cluster_cols = F,clustering_method = 'ward.D2',color = col.scheme,
                           breaks = breaks,
                           annotation_legend = F,show_colnames = F,show_rownames = F,fontsize_row =6,fontsize_col =6,
                           main = paste0('Semiparametric'))
save_pheatmap_png(heat4, paste0('result/overfitting/',name,'_heat_SemiSimulation.png'))








