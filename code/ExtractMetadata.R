setwd('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/data/')
files = list.files(pattern = 'Rdata$')
files
file = files[1]

df = matrix(data = NA, nrow = length(files), ncol = 4)
for(i in 1:length(files)){
  file = files[i]
  load(file)
  otu = nrow((phy@otu_table)@.Data)
  sample = ncol((phy@otu_table)@.Data)
  if(length(unique(phy@sam_data$HMP_BODY_SITE))==2){
    status = paste0(unique(phy@sam_data$HMP_BODY_SITE), collapse = ';')
  }else if(length(unique(phy@sam_data$grp)) ==2){
    status = paste0(unique(phy@sam_data$HMP_BODY_SUBSITE), collapse = ';')
  }
  df[i,1] = gsub('.Rdata','',file)
  df[i,2] = otu
  df[i,3] = sample
  df[i,4] = status
}


curated = as.data.frame(df)
colnames(curated) = c('dataset_name','taxa number','sample size','study_condition')


library(readxl)
methods_summary <- read_excel("~/Documents/Mayo_Research/SemiSimulation/methods_summary.xlsx", sheet = "Sheet3") %>% column_to_rownames('dataset_name')
head(methods_summary)
files = list.files(path = '/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/data/HMP/',pattern = 'Rdata$')
files.name = gsub('.Rdata','',files)

for(file in files){
  load(file)
  cat(file,'\n')
  name =gsub('.Rdata','',file)
  otu = nrow((phy@otu_table)@.Data)
  sample = ncol((phy@otu_table)@.Data)
  methods_summary[name,'sample size'] = sample
  methods_summary[name,'taxa number'] = otu
}
dim(methods_summary)

CRC.files = list.files(path = '/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/data/CRCdata/',pattern = 'Rdata$')
file = CRC.files[1]
for(file in CRC.files){
  load(paste0('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/data/CRCdata/',file))
  name = gsub('.Rdata','',file)
  otu = nrow((phy@otu_table)@.Data)
  sample = ncol((phy@otu_table)@.Data)
  
  methods_summary[name,'sample size'] = sample
  methods_summary[name,'taxa number'] = otu
}

methods_summary <- methods_summary %>% rownames_to_column('dataset_name')

CRC = list.files(path = '/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/data/CRCdata/',pattern = 'Rdata$')
CRC = gsub('.Rdata','',CRC)

methods_summary1 = methods_summary %>% filter(dataset_name %in% c(files.name,CRC))
# methods_summary1$`taxa number` = NA
# methods_summary1 = methods_summary1 %>% column_to_rownames('dataset_name')
# setwd('/Users/m216453/Documents/Mayo_Research/SemiSimulation/ReArrange1204/data/HMP/')
# for(file in files){
#   load(file)
#   name = gsub('.Rdata','',file)
#   otu = nrow((phy@otu_table)@.Data)
#   methods_summary1[name,'taxa number'] = otu
# }



methods_summary1 = methods_summary1 %>% rownames_to_column('dataset_name')
methods_summary1 = methods_summary1[,c("dataset_name","taxa number","sample size","study_condition","body_site","disease","age_category",
                                       "gender","country","sequencing_platform","PMID")]
head(curated)
head(methods_summary1)
m1 = methods_summary1[,c("dataset_name","taxa number","sample size","study_condition")]
m2 = methods_summary1 %>% dplyr::select(-c("taxa number","sample size","study_condition"))
data = rbind(m1, curated)
data = full_join(data, m2)

dim(data)
# data$dataset_name = gsub('^HanniganGD_2017_1$','HanniganGD_2017',data$dataset_name)
data$dataset_name = gsub('\\ ','_',data$dataset_name)

res.files = list.files(path = '~/Documents/Mayo_Research/SemiSimulation/ReArrange1204/result/realdata/result_preprocessing/',pattern = '_res.Rdata$')
res.files = gsub('_res.Rdata','',res.files)

sum(unique(data$dataset_name) %in% c(res.files,CRC))
data1 = data %>% filter(dataset_name %in% c(res.files,CRC)) # final datasets successfully analysis 
write.csv(data1, file = '/Users/m216453/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/result/realdataanalysis/DataSetInfo.csv', row.names = F)
