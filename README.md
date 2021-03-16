## A semi-parametric simulation framework for evaluation of microbial differential abundance analysis

### Main function of simulation function
```
SimulateSeq(otu.tab, model ='loglinear',
            nSam = 100, nOTU = 500, 
            diff.otu.pct = 0.1, diff.otu.direct = 'balanced', diff.otu.mode = 'abundant',
            include.top.otu = FALSE, k.top.otu = 5, 
            covariate.type = "binary", covariate.eff.mean = 0.5, 
            confounder.type = 'none', depth.mu = 10000, depth.conf.factor = 'none')
```

### Arguments
- `otu.tab`:     row are taxa, column are samples.   
- `nOTU` :    a number, indicates the number of taxa to be simulated.  
- `nSam`:     a number, indicates the number of samples to be simulated.   
- `diff.otu.pct `     a numerical fraction between 0 and 1, indicates the percentage of differential taxa.   
`diff.otu.direct`     Options include 'balanced', 'unbalanced'.  
`diff.otu.mode`     Options include 'abundant', 'rare', 'both'. 'abundant' means differential taxa come from top 1/4 abundant taxa, 'rare' means differential taxa come from tail 1/4 abundant taxa, 'both' means half of the differential taxa come from top 1/4 abundant taxa, while the rest half come from tail 1/4 abundant taxa.  
`covariate.eff.mean`     effect size to be created between case group and control group.   
`confounder.type`     covariate confounder. Options include 'none',  'continuous', 'binary', 'both'.   
`depth.mu`     mu to generate the sequencing depth.  
`depth.conf.factor`     depth confounding strength. Options include 'none'.  


### Value 
a list with components:
-  ***otu.tab.sim***, a data.frame of simulated taxa
- ***X***, a vector of simulated covariate
- ***diff.otu.ind***, a vector of *TRUE/FALSE* indicates the truth of simulated differential taxa
- ***otu.names***, a vector of simulated taxa names

### Simulation example
```
source('code/SimulationEvaluation/Simulation.R')
load('data/SimulationEvaluation/Vaginal.RData')
Sim.obj <- SimulateSeq(otu.tab = otu.tab, nOTU = 200, nSam = 50, covariate.eff.mean = 0, 
                       model = 'loglinear', nSub = 0, depth.conf.factor = 'none', covariate.type = 'binary')
```



