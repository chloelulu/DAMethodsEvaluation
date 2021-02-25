setwd('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/data/shuffle_data/v22/')
pkg =  c('patchwork','scales','stringr','ggpubr','microbiome',"eBay","modeest","ANCOMBC","aod","phyloseq","reshape","corncob","MASS","readr","DESeq2", "ALDEx2", "metagenomeSeq", "edgeR", "GUniFrac", "grDevices", "dirmult", "exactRankTests","nlme", "dplyr", "magrittr", "tidyr", "protoclust", "ggplot2", "compositions","rmutil","tibble","reticulate","dacomp","LDM","Wrench")
sapply(pkg, require, character = TRUE)
source('~/Documents/Mayo_project/Code/DailyCode.R')
files = list.files(getwd(),pattern = '.Rdata')
names = NULL
file = files[1]
for(file in files){
  a = rev(strsplit(file, "_")[[1]])[-c(1:2)]
  name = paste(rev(a), collapse = '_')
  # cat(file,':',name,'\n')
  names = c(names, name)
}
shuffles = unique(names)
# shuffles = shuffles[-grep("^$",shuffles)]
file = shuffles[1]
FDR <- list()
for(file in shuffles){
  fdrs = list()
  for(i in 1:50){
    try({
      load(paste0(file,'_',i,'_res.Rdata'))
      cat(file,'\n')
      fdr <- res.sum  %>% column_to_rownames('taxa') %>% dplyr::select(contains('fdr.')) #%>% replace(is.na(.), 1)
      fdr[,names(which(sapply(fdr, class)!='numeric'))] <- as.numeric(as.character(fdr[,names(which(sapply(fdr, class)!='numeric'))]))
      fdr[is.na(fdr)] <- 1
      fdr <- colSums(fdr <= 0.05, na.rm = T) %>% as.data.frame()
      colnames(fdr) = file
      fdr <- fdr %>% rownames_to_column('methods')
      colnames(fdr)[2] = i
      fdrs[[i]] = fdr
    }
    )
  }
  fdrss = reduce_list(fdrs)
  a = melt(fdrss, by = c('methods')) %>% dplyr::rename(iter = variable) # each method each iter findings
  # a = apply(fdrss %>% column_to_rownames('methods'),1, function(x) mean(x)) %>% as.data.frame()  %>% rownames_to_column('methods')# mean of each file in 50 iters
  colnames(a)[3] = file
  FDR[[file]] <- a
}
# # save(FDR, file = '/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/data/shuffle_data/result/v3_shuffle.Rdata')

for(i in 1:length(FDR)){
  cat(dim(FDR[[i]]),'\n')
}

load('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/data/shuffle_data/result/v1_shuffle.Rdata')
length(FDR)
FDR.df1 = reduce_list(data = FDR)
FDR.df1[is.na(FDR.df1)] = 0
load('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/data/shuffle_data/result/v3_shuffle.Rdata')
length(FDR)
FDR.df3 = reduce_list(data = FDR)
FDR.df3[is.na(FDR.df3)] = 0

load('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/data/shuffle_data/result/v2_shuffle.Rdata')
length(FDR)
FDR.df2 = reduce_list(FDR)
FDR.df2[is.na(FDR.df2)] = 0
head(FDR.df2)
# FDR.df = FDR.df2
# FDR.df = full_join(as.data.frame(FDR.df1), as.data.frame(FDR.df2))

FDR.df$methods = gsub('fdr.','',FDR.df$methods)
res1 = melt(FDR.df) %>% na.omit()
delete = c('DESeq2','Wrench+DESeq2','edgeR','Wrench+edgeR','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon','RioNorm2','LINDA2','LinDA')

res1$methods = as.character(res1$methods)
res1$methods[res1$methods =='Wilcox'] = 'TSS+Wilcoxon'
res1$methods[res1$methods =='ttest'] = 'TSS+t-test'
res1$methods[res1$methods =='Rarefy'] = 'Rarefy+Wilcoxon'
res1$methods = gsub('ANCOMBC','ANCOM-BC',res1$methods)
res1$methods = gsub('glmquassi','GLM(quasipoisson)',res1$methods)
res1$methods = gsub('DESeq2.gmpr','GMPR+DESeq2',res1$methods)
res1$methods = gsub('DESeq2.Wrench','Wrench+DESeq2',res1$methods)
res1$methods = gsub('edgeR.gmpr','GMPR+edgeR',res1$methods)
res1$methods = gsub('edgeR.Wrench','Wrench+edgeR',res1$methods)
res1$methods = gsub('MSeq2.Wrench','Wrench+metagenomeSeq',res1$methods)
res1$methods = gsub('MSeq2','metagenomeSeq',res1$methods)
res1$methods = gsub('eBayW','eBay(Wilcoxon)',res1$methods)
res1$methods = gsub('eBayt','eBay(t-test)',res1$methods)
res1$methods = gsub('BBinomial','Beta-binomial',res1$methods)
res1$methods = gsub('Aldex2we','Aldex2(t-test)',res1$methods)
res1$methods = gsub('^Aldex2$','Aldex2(Wilcoxon)',res1$methods)
res1$methods = gsub('^CLRBC$','LinDA',res1$methods)
res1$methods = gsub('^Rarefyttest$','Rarefy+t-test',res1$methods)
res1$methods = gsub('^ttest$','TSS+t-test',res1$methods)

df1 = res1 %>% filter(!(methods %in% delete))
unique(df1$methods)

df1$methods = as.character(df1$methods)
df1 = df1[order(df1$methods),]
uniq = unique(df1$methods)
for(i in 1:length(unique(df1$methods))){
  df1[(df1$methods ==uniq[i]),'label'] = letters[i]
}
df1$legend = paste0(df1$label,':',df1$methods)
md = aggregate( value ~ methods,df1, function(x) median(x))
md = md[order(md$value),]
ord<- within(df1, methods <- factor(methods, levels=md$methods))
head(df1)
size = 6
unique(ord1$methods)


cols1 = cols

s = unique(df1$legend)
for(i in 1:length(cols1)){
  if(names(cols1)[i] %in% gsub('.*:','',s)){
    names(cols1)[i] = s[which((names(cols1)[i] == gsub('.*:','',s)))]
  }
}

p1 = ggplot(ord, aes(x =reorder(label, value, FUN = median), y = value, fill= legend)) +
  theme_bw() +
  stat_boxplot(geom = "errorbar", width = 0.2,lwd=0.1) +
  geom_boxplot(width=0.7,outlier.size = 0.2, outlier.colour = 'grey30',lwd=0.1) +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1) +
  ylab('False discoveries') +
  xlab('') +
  labs(fill = '')+
  theme(axis.text.x = element_text(color="black", size =26),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color="black", size = 26),
        axis.title = element_text(color="black", size = 30),
        strip.text = element_text(size = 30),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        plot.title = element_text(size = 40,hjust = -0.2),
        plot.margin = rep(grid::unit(0.75,"in"),4),
        legend.position = 'none',
        legend.title=element_text(size=30),
        legend.text = element_text(size=30))+
  ggtitle("(a)") +
  scale_y_continuous(trans = sqrt_trans(),
                     breaks = trans_breaks("sqrt", function(x) x^2),
                     labels = trans_format("sqrt", math_format(.x^2))) +
  geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5)
p1
#ggsave('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/plot/Shuffle_Discoveries.pdf', width = 9, height = 7, dpi = 100)


head(ord)
ord1 = ord
ord1$value[ord1$value >0] = 1
head(ord1)
ord2 = aggregate(value ~ legend +methods+label +variable, ord1, function(x) sum(x)/100)
ord2 = data_summary(ord2, formula = value ~ legend+label+methods)
# ord2 = aggregate(value ~ methods + variable, ord1, function(x) mean(x))
head(ord2)
p2 = ggplot(ord2, aes(x = reorder(label, value), y = value, fill= legend)) +
  theme_bw() +
  geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7) +
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.4, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1) +
  scale_y_continuous(breaks = c(0,0.005,0.01))+
  ylab('FDIs') +
  xlab('') +
  labs(fill = '')+
  theme(axis.text.x = element_text(color="black", size =26),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color="black", size = 26),
        axis.title = element_text(color="black", size = 30),
        strip.text = element_text(size = 30),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        legend.position = 'none',
        legend.title=element_text(size=30),
        legend.text = element_text(size=30),
        plot.title = element_text(size=22))
p1 + p2
# ggsave('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/plot/Shuffle_DiscoveriesFDIs.pdf', width = 15, height = 8, dpi = 100)
# ggsave('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/Shuffle_DiscoveriesFDIs.pdf', width = 15, height = 8, dpi = 100)

## Shuffled Depth confounding data, more findings more errors
setwd('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/realdata/depth/')
# files = list.files(getwd(), pattern = 'Rdata$')
# names = NULL
# for(file in files){
#   a = rev(strsplit(file, "_")[[1]])[-c(1:2)]
#   name = paste(rev(a), collapse = '_')
#   # cat(file,':',name,'\n')
#   names = c(names, name)
# }
# shuffles = unique(names)
# 
# FDR <- list()
# for(file in shuffles){
#   fdrs = list()
#   for(i in 1:50){
#     try({
#       load(paste0(file,'_',i,'_res.Rdata'))
#       cat(file,'\n')
#       fdr <- res.sum  %>% column_to_rownames('taxa') %>% dplyr::select(contains('fdr.')) #%>% replace(is.na(.), 1)
#       fdr[,names(which(sapply(fdr, class)!='numeric'))] <- as.numeric(as.character(fdr[,names(which(sapply(fdr, class)!='numeric'))]))
#       fdr[is.na(fdr)] <- 1
#       fdr <- colSums(fdr <= 0.05, na.rm = T) %>% as.data.frame()
#       colnames(fdr) = file
#       fdr <- fdr %>% rownames_to_column('methods')
#       colnames(fdr)[2] = i
#       fdrs[[i]] = fdr
#     }
#     )
#   }
#   fdrss = reduce_list(fdrs)
#   a = melt(fdrss, by = c('methods')) %>% dplyr::rename(iter = variable)
#   # a = apply(fdrss %>% column_to_rownames('methods'),1, function(x) mean(x)) %>% as.data.frame() # mean of each file in 50 iters
#   colnames(a)[3] = file
#   FDR[[file]] <- a
# }
# save(FDR, file = '/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/depthconfounding.Rdata')
getwd()


load('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/depthconfounding.Rdata')
FDR.df = reduce_list(FDR)


FDR.df$methods = gsub('fdr.','',FDR.df$methods)
head(FDR.df)
res1 = melt(FDR.df) %>% na.omit()
head(res1)
delete = c('DESeq2','Wrench+DESeq2','edgeR','Wrench+edgeR','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon','RioNorm2','LINDA2','LinDA')

res1$methods = as.character(res1$methods)
res1$methods[res1$methods =='Wilcox'] = 'TSS+Wilcoxon'
res1$methods[res1$methods =='ttest'] = 'TSS+t-test'
res1$methods[res1$methods =='Rarefy'] = 'Rarefy+Wilcoxon'
res1$methods = gsub('ANCOMBC','ANCOM-BC',res1$methods)
res1$methods = gsub('glmquassi','GLM(quasipoisson)',res1$methods)
res1$methods = gsub('DESeq2.gmpr','GMPR+DESeq2',res1$methods)
res1$methods = gsub('DESeq2.Wrench','Wrench+DESeq2',res1$methods)
res1$methods = gsub('edgeR.gmpr','GMPR+edgeR',res1$methods)
res1$methods = gsub('edgeR.Wrench','Wrench+edgeR',res1$methods)
res1$methods = gsub('MSeq2.Wrench','Wrench+metagenomeSeq',res1$methods)
res1$methods = gsub('MSeq2','metagenomeSeq',res1$methods)
res1$methods = gsub('eBayW','eBay(Wilcoxon)',res1$methods)
res1$methods = gsub('eBayt','eBay(t-test)',res1$methods)
res1$methods = gsub('BBinomial','Beta-binomial',res1$methods)
res1$methods = gsub('Aldex2we','Aldex2(t-test)',res1$methods)
res1$methods = gsub('^Aldex2$','Aldex2(Wilcoxon)',res1$methods)
res1$methods = gsub('^CLRBC$','LinDA',res1$methods)
res1$methods = gsub('^Rarefyttest$','Rarefy+t-test',res1$methods)
res1$methods = gsub('^ttest$','TSS+t-test',res1$methods)

df1 = res1 %>% filter(!(methods %in% delete))
unique(df1$methods)

df1$methods = as.character(df1$methods)
df1 = df1[order(df1$methods),]
uniq = unique(df1$methods)
for(i in 1:length(unique(df1$methods))){
  df1[(df1$methods ==uniq[i]),'label'] = letters[i]
}
df1$legend = paste0(df1$label,':',df1$methods)
md = aggregate( value ~ methods,df1, function(x) median(x))
md = md[order(md$value),]
ord<- within(df1, methods <- factor(methods, levels=md$methods))
head(df1)
size = 6
unique(ord1$methods)


d1 = ggplot(ord, aes(x =reorder(label, value, FUN = median), y = value, fill= legend)) +
  theme_bw() +
  stat_boxplot(geom = "errorbar", width = 0.2,lwd=0.1) +
  geom_boxplot(width=0.7,outlier.size = 0.2, outlier.colour = 'grey30',lwd=0.1) +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1) +
  ylab('False discoveries') +
  xlab('') +
  labs(fill = '')+
  theme(axis.text.x = element_text(color="black", size =26),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color="black", size = 26),
        axis.title = element_text(color="black", size = 30),
        strip.text = element_text(size = 30),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        plot.title = element_text(size = 40,hjust = -0.2),
        plot.margin = rep(grid::unit(0.75,"in"),4),
        legend.position = 'none',
        legend.title=element_text(size=30),
        legend.text = element_text(size=30))+
  ggtitle("(c)") +
  scale_y_continuous(trans = sqrt_trans(),
                     breaks = trans_breaks("sqrt", function(x) x^2),
                     labels = trans_format("sqrt", math_format(.x^2))) +
  geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed', size =  0.5)


head(ord)


head(ord)
ord1 = ord
ord1$value[ord1$value >0] = 1
head(ord1)
ord2 = aggregate(value ~ legend +methods+label +variable, ord1, function(x) sum(x)/100)
ord2 = data_summary(ord2, formula = value ~ legend+label+methods)
# ord2 = aggregate(value ~ methods + variable, ord1, function(x) mean(x))
head(ord2)
d2 = ggplot(ord2, aes(x = reorder(label, value), y = value, fill=legend)) +
  theme_bw() +
  geom_bar(position = position_dodge2(width = 0.7, preserve = "single"), stat="identity", width = 0.7) +
  geom_errorbar(aes(ymax=ymax, ymin=ymin, linetype = NULL),  width=0.4, size = 0.2,position=position_dodge2(.7, preserve = "single")) +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1) +
  ylab('FDIs') +
  xlab('') +
  labs(fill = '')+
  theme(axis.text.x = element_text(color="black", size =26),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color="black", size = 26),
        axis.title = element_text(color="black", size = 30),
        strip.text = element_text(size = 30),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        legend.position = 'none',
        legend.title=element_text(size=30),
        legend.text = element_text(size=30),
        plot.title = element_text(size=22))
d1 +d2

# ggsave('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/plot/Shuffle_Depth_Discoveries.pdf', width =15, height = 7, dpi = 100)
# ggsave('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/Shuffle_Depth_Discoveries.pdf', width =15, height = 8, dpi = 100)

## majority votes
pct = 0.2
load(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/realdata/FDRs_realdata_',pct,'.Rdata'))
load(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/realdata/TPRs_realdata_',pct,'.Rdata'))
func <- function(FDRs){
  FDR1 <- as.data.frame(do.call(rbind, FDRs)) %>% melt(.)
  colnames(FDR1)[1]= 'methods'
  FDR1$methods = as.character(FDR1$methods)
  FDR1$methods[FDR1$methods =='Wilcox'] = 'TSS+Wilcoxon'
  FDR1$methods[FDR1$methods =='Rarefy'] = 'Rarefy+Wilcoxon'
  FDR1$methods = gsub('ANCOMBC','ANCOM-BC',FDR1$methods)
  FDR1$methods = gsub('glmquassi','GLM(quasipoisson)',FDR1$methods)
  FDR1$methods = gsub('DESeq2.gmpr','GMPR+DESeq2',FDR1$methods)
  FDR1$methods = gsub('DESeq2.Wrench','Wrench+DESeq2',FDR1$methods)
  FDR1$methods = gsub('edgeR.gmpr','GMPR+edgeR',FDR1$methods)
  FDR1$methods = gsub('edgeR.Wrench','Wrench+edgeR',FDR1$methods)
  FDR1$methods = gsub('MSeq2.Wrench','Wrench+metagenomeSeq',FDR1$methods)
  FDR1$methods = gsub('MSeq2','metagenomeSeq',FDR1$methods)
  FDR1$methods = gsub('eBayW','eBay(Wilcoxon)',FDR1$methods)
  FDR1$methods = gsub('eBayt','eBay(t-test)',FDR1$methods)
  FDR1$methods = gsub('BBinomial','Beta-binomial',FDR1$methods)
  FDR1$methods = gsub('Aldex2we','Aldex2(t-test)',FDR1$methods)
  FDR1$methods = gsub('^Aldex2$','Aldex2(Wilcoxon)',FDR1$methods)
  FDR1$methods = gsub('^Rarefyttest$','Rarefy+t-test',FDR1$methods)
  FDR1$methods = gsub('^ttest$','TSS+t-test',FDR1$methods)
  FDR1$methods = gsub('^Omnibus$','mbzinb',FDR1$methods)
  cols <- c('ZicoSeq'=brewer.pal(9,'Paired')[5], 'Omnibus' = brewer.pal(9,'PuBuGn')[7], 'mbzinb' = brewer.pal(9,'PuBuGn')[8], 
            'GLM(quasipoisson)' = brewer.pal(9,'PuRd')[3],'LinDA'=brewer.pal(9,'Set1')[8],'LINDA2'=brewer.pal(9,'Set1')[7],
            'RAIDA'=brewer.pal(11,'BrBG')[7],'RioNorm2' = brewer.pal(11,'BrBG')[8],
            'eBay(Wilcoxon)' = brewer.pal(8,'RdPu')[4],'eBay(t-test)' = brewer.pal(8,'RdPu')[6], 'ANCOM-BC' = brewer.pal(8,'Set2')[7],
            'LDM'=brewer.pal(9,'Set1')[1], 'DACOMP'=brewer.pal(8,'Dark2')[6],'Aldex2(Wilcoxon)'=brewer.pal(8,'YlGn')[3],'Aldex2(t-test)'=brewer.pal(8,'YlGn')[2],
            'Beta-binomial'=brewer.pal(11,'BrBG')[2],'DESeq2'=brewer.pal(9,'Greys')[7], 'Wrench+DESeq2'= brewer.pal(9,'Greys')[3], 'GMPR+DESeq2'=brewer.pal(9,'Greys')[5],
            'TSS+Wilcoxon'=brewer.pal(9,'Purples')[7],'Rarefy+Wilcoxon'= brewer.pal(9,'Set1')[5], 'GMPR+Wilcoxon'= brewer.pal(9,'Purples')[5], 'Wrench+Wilcoxon'= brewer.pal(9,'Purples')[3], 
            'TSS+Spearman'=brewer.pal(9,'Purples')[7],'Rarefy+Spearman'=brewer.pal(9,'Set1')[5], 'GMPR+Spearman'= brewer.pal(9,'Purples')[5], 
            'TSS+t-test'=brewer.pal(9,'GnBu')[7],'Rarefy+t-test'= brewer.pal(9,'GnBu')[5], 'GMPR+t-test'= brewer.pal(9,'GnBu')[5], 'Wrench+t-test'= brewer.pal(9,'GnBu')[3], 
            'edgeR'=brewer.pal(9,'Greens')[7], 'Wrench+edgeR'=brewer.pal(9,'Greens')[3], 'GMPR+edgeR'=brewer.pal(9,'Greens')[5], 
            'metagenomeSeq'=brewer.pal(9,'Blues')[7], 'Wrench+metagenomeSeq'=brewer.pal(9,'Blues')[5])
  
  mn = aggregate(value~methods, FDR1 , function(x) median(x[!is.na(x)]))
  mn = mn[order(mn$value),]
  FDR1 <- within(FDR1, methods <- factor(methods, levels=mn$methods))
  delete = c('DESeq2','Wrench+DESeq2','edgeR','Wrench+edgeR','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon','RioNorm2','LINDA2','LinDA')
  
  FDR1 <- FDR1 %>% filter(!(methods %in% delete))
}
data1 <- func(FDRs = FDRs)
data2 <- func(FDRs = TPRs)
head(data1)
data1$methods = as.character(data1$methods)
data1 = data1[order(data1$methods),]
uniq = unique(data1$methods)
for(i in 1:length(unique(data1$methods))){
  data1[(data1$methods ==uniq[i]),'label'] = letters[i]
}
data1$legend = paste0(data1$label,':',data1$methods)
head(data1)
md = aggregate( value ~ label,data1, function(x) median(x))
md = md[order(md$value),]
data1$label = as.factor(data1$label)
data1 <- within(data1, label <- factor(label, levels=md$label))
head(data1)
unique(data1$legend)
v1 = ggplot(data1,aes(x = label, y = value, fill =legend)) +
  stat_boxplot(geom = "errorbar", width = 0.5,lwd=0.2) +
  geom_boxplot(aes(fill = legend),outlier.size = 0.7, outlier.colour = 'grey', width = 0.9) +
  theme_bw() +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1) +
  theme(axis.text.x = element_text(color="black", size =26),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color="black", size = 26),
        axis.title = element_text(color="black", size = 30),
        strip.text = element_text(size = 30),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        plot.title = element_text(size = 40,hjust = -0.2),
        plot.margin = rep(grid::unit(0.75,"in"),4),
        legend.position = 'none',
        legend.title=element_text(size=30),
        legend.text = element_text(size=30))+
  ggtitle("(b)") +
  ylab('FDR') + xlab('') + guides(fill = F) +
  geom_hline(yintercept = 0.05, colour = brewer.pal(8,'Set1')[3], linetype = 'dashed', size = 0.5) +
  geom_hline(yintercept = 0.1, colour = brewer.pal(8,'Set2')[6], linetype = 'dashed', size = 0.5)  +
  geom_hline(yintercept = 0.2, colour = brewer.pal(8,'Set1')[1], linetype = 'dashed', size = 0.5) 
v1


head(data2)
(aggregate(value ~ methods, data2, function(x) median(x)) %>% filter(value <=0.7 &value >0.5))[,1]
(aggregate(value ~ methods, data2, function(x) mean(x)) %>% filter(value >0.6))[,1]
data2$methods = as.character(data2$methods)
data2 = data2[order(data2$methods),]
uniq = unique(data2$methods)
for(i in 1:length(unique(data2$methods))){
  data2[(data2$methods ==uniq[i]),'label'] = letters[i]
}
data2$legend = paste0(data2$label,':',data2$methods)
head(data2)
md = aggregate( value ~ label,data2, function(x) median(x))
md = md[order(md$value),]
data2$label = as.factor(data2$label)
data2 <- within(data2, label <- factor(label, levels=md$label))
head(data2)
v2 = ggplot(data2,aes(x =label, y = value, fill = legend)) +
  stat_boxplot(geom = "errorbar", width = 0.5,lwd=0.2) +
  geom_boxplot(aes(fill = legend),outlier.size =0.7, outlier.colour = 'grey', width = 0.9) +
  theme_bw() +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1) +
  theme(axis.text.x = element_text(color="black", size =26),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color="black", size = 26),
        axis.title = element_text(color="black", size = 30),
        strip.text = element_text(size = 30),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        legend.position = 'none',
        legend.title=element_text(size=30),
        legend.text = element_text(size=30),
        plot.title = element_text(size=22))+
  ylab('Power') + xlab('') + guides(fill = F) 

v3 = ggplot(data2,aes(x =legend, y = value)) +
  stat_boxplot(geom = "errorbar", width = 0.5,lwd=0.2) +
  geom_boxplot(aes(fill = legend),outlier.size =0.7, outlier.colour = 'grey', width = 0.9) +
  theme_bw() +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1) +
  theme(axis.text.x = element_blank(),#element_text(color="black", size =26, angle = 90, hjust = 1, vjust = 0.25),
        axis.ticks = element_blank(),
        axis.text.y = element_text(color="black", size = 26),
        axis.title = element_text(color="black", size = 30),
        strip.text = element_text(size = 30),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        panel.border = element_rect(colour = "black",size= 1),
        legend.title=element_text(size=30),
        legend.text = element_text(size=30),
        plot.title = element_text(size=22))+
  labs(x = '',y = 'Power', fill = '')
leg <- get_legend(v3)
legend = as_ggplot(leg)
ggsave(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/truthfinding_composites_legend.pdf'), width = 22, height = 5, dpi = 100)

v1 +v2
# ggsave(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/FDRs_realdata_',pct,'.pdf'), width = 12, height = 8, dpi = 100)
((p1 + p2)/(v1 +v2)/(d1 +d2))|legend
ggsave(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/truthfinding_composites.pdf'), width = 22, height = 15, dpi = 100)



