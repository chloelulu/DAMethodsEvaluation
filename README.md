# A semi-parametric simulation framework for evaluation of microbial differential abundance analysis

## Main function of simulation function
```
SimulateSeq(otu.tab, model ='loglinear',
            nSam = 100, nOTU = 500, 
            diff.otu.pct = 0.1, diff.otu.direct = 'balanced', diff.otu.mode = 'abundant',
            include.top.otu = FALSE, k.top.otu = 5, 
            covariate.type = "binary", covariate.eff.mean = 0.5, 
            confounder.type = 'none', depth.mu = 10000, depth.conf.factor = 'none')
```

Arguments
`otu.tab` : row are taxa, column are samples
``


## Simulation example
```
source('code/SimulationEvaluation/Simulation.R')
load('data/SimulationEvaluation/Vaginal.RData')
Sim.obj <- SimulateSeq(otu.tab = otu.tab, nOTU = 200, nSam = 50, covariate.eff.mean = 0, 
                       model = 'loglinear', nSub = 0, depth.conf.factor = 'none', covariate.type = 'binary')
```


