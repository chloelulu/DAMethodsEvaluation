# A semi-parametric simulation framework for evaluation of microbial differential abundance analysis

## Main function of simulation function
```
SimulateSeq(otu.tab,
            model = c('loglinear','heterogeneity','nonmonotone'), nSub = 0.3,
            nSam = 100, nOTU = 500, 
            diff.otu.pct = 0.1, diff.otu.direct = c('balanced','unbalanced'), diff.otu.mode = c('abundant','rare','mix'), include.top.otu = FALSE, k.top.otu = 5, 
            covariate.type = c("binary","continuous"), covariate.eff.mean = 0.5, grp.ratio = 1, covariate.eff.sd = 0, weight.abund.func = function (x) 0,
            confounder.type = c('none',"binary","continuous","both"), conf.cov.cor = 0.6, conf.diff.otu.pct = 0.05, conf.nondiff.otu.pct = 0.1, confounder.eff.mean = 0, confounder.eff.sd = 0,
            depth.mu = 10000, depth.theta = 5,  depth.conf.factor = c('none','DL1','DL2','DL3'))
```

Arguments
```
otu.tab
```
row: OTU, column: sample


## Example of simulation with semi-paramteric simulation framework
```
setwd('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/')
source('code/SimulationEvaluation/DM1.R')
source('code/SimulationEvaluation/wrappers.R')
load('data/SimulationEvaluation/Vaginal.RData')
Sim.obj <- SimulateSeq(otu.tab = otu.tab, nOTU = 200, nSam = 50, covariate.eff.mean = 0, 
                       model = 'loglinear', nSub = 0, depth.conf.factor = 'none', covariate.type = 'binary')
```


