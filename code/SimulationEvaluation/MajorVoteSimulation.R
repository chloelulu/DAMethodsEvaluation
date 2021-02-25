setwd("/Users/m216453/Documents/Mayo_Research/SemiSimulation/0925/sim/new/sim/allC1/")
source('~/Documents/Mayo_project/code/DailyCode.R')
files = list.files(getwd(), pattern = 'allC1_summarymatrix*')
files
load(files[j])
names = gsub(gsub('.Rdata','',gsub('.*\\-','',files[j])),'',names(res_seqs))
names = names[grep('binary',names)] # caution, only extract binary

file = files[1]
methods = NULL
for (file in files) {
  name = gsub('allC1_summarymatrix1-','',file)
  name = gsub('.Rdata','',name)
  name = gsub('.*\\-','',name)
  name
  methods = c(methods, name)
}
methods = methods[-grep('edgeR|DESeq2|RioNorm2|Wilcox.gmpr|Wilcox.Wrench|ttest.gmpr|ttest.Wrench|Omnibus',methods)]
methods = unique(methods)

# each filename in different method  has 1 voted FDR/TPR, cal mean of 100 iterations
# i = iters; methods = methods; names[1] means names
total = list()
for(pct in seq(0.2,0.8,0.2)){
  dfs = list()
  for(k in 1:length(names)){
    FDR.dfs = list()
    for(i in 1:100){
      FDR  = list()
      for(method in methods){
        
        load(paste0("allC1_summarymatrix",i,"-",method,".Rdata"))
        df  = res_seqs[[paste0(names[k],method)]] # one file name
        if(method %in% c('Aldex2','Aldex2we')){
          fdr = 'BH'
        }else{
          fdr = 'fdr'
        }
        
        if(is.null(df)){
          cat(method,'\n')
        }
        if(!(is.null(df))){
          df = df %>% dplyr::select(c('otu.id',fdr))
          colnames(df)[2] = method
          FDR[[method]] = df
        }
      }
      FDR.df = reduce_list(FDR)
      FDR.df[is.na(FDR.df)] = 1
      
      # extract the truth ids 
      load(paste0("allC1_summarymatrix",i,"-Wilcox.Rdata")) 
      truth = res_seqs[[paste0(names[k],method)]] %>% dplyr::select(c('otu.id', 'diff.otu.ind'))# res_seqs
      test1 = full_join(truth,FDR.df)
      # calculate voted/truth FDR/TPR of names[1]
      cutoff = 0.05
      test1$truth = rowSums(test1[,c(3:ncol(test1))] <= cutoff)/(ncol(test1)-3)
      
      vote.fdr = vote.tpr = truth.fdr = truth.tpr = NULL
      for(method in methods){
        try({
          tp <- sum(test1[,method] <= cutoff & test1$truth > pct, na.rm = T)
          tn <- sum(test1[,method] > cutoff &  test1$truth < pct, na.rm = T)
          fn <- sum(test1[,method] > cutoff &  test1$truth > pct, na.rm = T)
          fp <- sum(test1[,method] <= cutoff &  test1$truth < pct, na.rm = T)
          tp0 <- sum(test1[,method] <= cutoff & test1$diff.otu.ind == TRUE, na.rm = T)
          tn0 <- sum(test1[,method] > cutoff &  test1$diff.otu.ind == FALSE, na.rm = T)
          fn0 <- sum(test1[,method] > cutoff &  test1$diff.otu.ind == TRUE, na.rm = T)
          fp0 <- sum(test1[,method] <= cutoff &  test1$diff.otu.ind == FALSE, na.rm = T)
          tpr <- ifelse(!tp, 0, tp/(tp+fn))
          fdr <- ifelse(!fp, 0, fp/(tp+fp))
          
          tpr0 <- ifelse(!tp0, 0, tp0/(tp0+fn0))
          fdr0 <- ifelse(!fp0, 0, fp0/(tp0+fp0))
          names(tpr) =names(fdr) = names(tpr0) =names(fdr0) = method
          
          vote.fdr = c(vote.fdr,fdr)
          vote.tpr = c(vote.tpr,tpr)
          
          truth.fdr = c(truth.fdr,fdr0)
          truth.tpr = c(truth.tpr,tpr0)
          
          res = cbind(vote.fdr, vote.tpr, truth.fdr, truth.tpr)
          
        })
      }
      
      FDR.dfs[[i]] = res # all methods in all iterations
    }
    
    vote.fdrs = vote.tprs = truth.fdrs = truth.tprs = NULL
    for(i in 1:100){
      test1 = FDR.dfs[[i]]
      vote.fdrs = cbind(vote.fdrs,test1[,'vote.fdr'])
      vote.tprs = cbind(vote.tprs,test1[,'vote.tpr'])
      truth.fdrs = cbind(truth.fdrs,test1[,'truth.fdr'])
      truth.tprs = cbind(truth.tprs,test1[,'truth.tpr'])
    }
    df = cbind(rowMeans(vote.fdrs),rowMeans(vote.tprs),rowMeans(truth.fdrs),rowMeans(truth.tprs))
    colnames(df) = c('vote.fdr','vote.tpr','truth.fdr','truth.tpr')
    dfs[[names[k]]] = df
  }

  
  
  vote.fdrs = vote.tprs = truth.fdrs = truth.tprs = NULL
  for(m in 1:length(dfs)){
    test1 = dfs[[m]]
    vote.fdrs = cbind(vote.fdrs,test1[,'vote.fdr'])
    vote.tprs = cbind(vote.tprs,test1[,'vote.tpr'])
    truth.fdrs = cbind(truth.fdrs,test1[,'truth.fdr'])
    truth.tprs = cbind(truth.tprs,test1[,'truth.tpr'])
  }
  summary = cbind(rowMeans(vote.fdrs),rowMeans(vote.tprs),rowMeans(truth.fdrs),rowMeans(truth.tprs))
  colnames(summary) = c('vote.fdr','vote.tpr','truth.fdr','truth.tpr')
  save(summary, file = paste0('summaryOfVote',pct,'.Rdata'))
  total[pct] = summary
}


load(paste0('summaryOfVote0.2.Rdata'))
y = melt(summary) %>% separate(Var2,c('vote','grp'))
y$vote = gsub('vote','20% vote',y$vote)

load(paste0('summaryOfVote0.4.Rdata'))
y1 = melt(summary)%>% separate(Var2,c('vote','grp')) %>% filter(vote == 'vote')
y1$vote = gsub('vote','40% vote',y1$vote)


load(paste0('summaryOfVote0.6.Rdata'))
y2 = melt(summary)%>% separate(Var2,c('vote','grp')) %>% filter(vote == 'vote')
y2$vote = gsub('vote','60% vote',y2$vote)


load(paste0('summaryOfVote0.8.Rdata'))
y3 = melt(summary)%>% separate(Var2,c('vote','grp')) %>% filter(vote == 'vote')
y3$vote = gsub('vote','80% vote',y3$vote)

yy = rbind(y, y1) %>% rbind(y2) %>% rbind(y3) 
head(yy)
yy$grp = gsub('fdr','FDR',yy$grp)
yy$grp = gsub('tpr','Power',yy$grp)


yy$Var1 = as.character(yy$Var1)
yy$Var1[yy$Var1 =='Wilcox'] = 'TSS+Wilcoxon'
yy$Var1[yy$Var1 =='ttest'] = 'TSS+t-test'

yy$Var1[yy$Var1 =='Rarefy'] = 'Rarefy+Wilcoxon'
yy$Var1 = gsub('ANCOMBC','ANCOM-BC',yy$Var1)
yy$Var1 = gsub('glmquassi','GLM(quasipoisson)',yy$Var1)
yy$Var1 = gsub('DESeq2.gmpr','GMPR+DESeq2',yy$Var1)
yy$Var1 = gsub('DESeq2.Wrench','Wrench+DESeq2',yy$Var1)
yy$Var1 = gsub('edgeR.gmpr','GMPR+edgeR',yy$Var1)
yy$Var1 = gsub('edgeR.Wrench','Wrench+edgeR',yy$Var1)
yy$Var1 = gsub('MSeq2.Wrench','Wrench+metagenomeSeq',yy$Var1)
yy$Var1 = gsub('MSeq2','metagenomeSeq',yy$Var1)
yy$Var1 = gsub('eBayW','eBay(Wilcoxon)',yy$Var1)
yy$Var1 = gsub('eBayt','eBay(t-test)',yy$Var1)
yy$Var1 = gsub('BBinomial','Beta-binomial',yy$Var1)
yy$Var1 = gsub('Aldex2we','Aldex2(t-test)',yy$Var1)
yy$Var1 = gsub('^Aldex2$','Aldex2(Wilcoxon)',yy$Var1)
yy$Var1 = gsub('^CLRBC$','LinDA',yy$Var1)
yy$Var1 = gsub('^Rarefyttest$','Rarefy+t-test',yy$Var1)

head(yy)
yy$vote = gsub('truth','Truth',yy$vote) 
yy$vote <- factor(yy$vote, levels = c('Truth','20% vote','40% vote','60% vote','80% vote'))
delete = c('DESeq2','Wrench+DESeq2','edgeR','Wrench+edgeR','GMPR+t-test','Wrench+t-test','GMPR+Wilcoxon','Wrench+Wilcoxon','RioNorm2','LINDA2','LinDA')

yy = yy %>% filter(!(Var1 %in% delete))
ggplot(yy, aes(x = Var1, y = value, fill = Var1)) +
  geom_bar(stat="identity",position = position_dodge2(width = 0.9, preserve = "single"),  width = 0.9) +
  scale_fill_manual(values = cols) +
  facet_grid(vote ~grp, scales = 'free_y') +
  geom_smooth(aes(group = vote), color = 'red', se = F, size =2) +
  theme_bw()+
  labs(x = '', fill = '', y = '') +
  theme(axis.text.x = element_text(color="black", size =26, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(color="black", size = 26),
        axis.title = element_text(color="black", size = 30),
        strip.text = element_text(size = 30),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        legend.title=element_text(size=30),
        legend.text = element_text(size=30), legend.position = 'none',
        plot.title = element_text(size=22))
# ggsave('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/plot/vote.pdf', width = 14, height =18, dpi = 100)
ggsave(filename = '/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/VoteOfSimulation.pdf', width = 14, height = 18, dpi = 100)




