library(tibble); library(dplyr)

## We calculate compositionality all simulated dataset allD5 on cluster.
setwd('~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/calcomp/')
dir = "~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/calcomp"
files = list.files(dir)
DF = NULL
for(i in 1:length(files)){
  file = files[i]
  load(file)
  df$iter = paste0('iter',i)
  DF = rbind(DF, df)
}

colnames(DF)[1:2] = c('Effect size','Compositionality')
levels(DF$signal.denity)
DF$signal.denity = gsub('low','Low',DF$signal.denity)
DF$signal.denity = gsub('medium','Medium',DF$signal.denity) 
DF$signal.denity = gsub('high','High',DF$signal.denity)
DF$diff.otu.mode = gsub('abundant','Abundant',DF$diff.otu.mode)
DF$diff.otu.mode = gsub('rare','Rare',DF$diff.otu.mode)
DF$signal.denity = as.factor(DF$signal.denity)
levels(DF$signal.denity)
DF <- within(DF, signal.denity <- factor(signal.denity, levels=c('Low','Medium','High')))

DF$`Effect size` = as.factor(DF$`Effect size`)
levels(DF$`Effect size`) <- c('Weak','Moderate','Strong')

ggplot(DF, aes(x = diff.otu.mode, y = Compositionality, fill = `Effect size`)) +
  geom_boxplot(outlier.size = 0.2, width = 0.5) +
  facet_grid(.~signal.denity) +
  scale_fill_brewer(palette = 'Set2') +
  theme_bw() +
  theme(axis.text.x = element_text(color="black", size =20),
        axis.text.y = element_text(color="black", size = 20),
        axis.title = element_text(color="black", size = 20),
        strip.text = element_text(size = 20),
        strip.background = element_rect(fill="white",color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size= 1),
        legend.title=element_text(size=20),
        legend.text = element_text(size=20),
        plot.title = element_text(size=20)) +
  labs(fill = 'Effect size') +
  xlab('')
ggsave("~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/Compostionality.pdf", width = 10, height = 3, dpi = 100)

