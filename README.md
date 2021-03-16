## A semi-parametric simulation framework for evaluation of microbial differential abundance analysis

### Main function of simulation function
```
SimulateSeq(otu.tab, model ='loglinear',
            nSam = 100, nOTU = 500,
            diff.otu.pct = 0.1, diff.otu.direct = 'balanced', diff.otu.mode = 'abundant',
            include.top.otu = FALSE, k.top.otu = 5,
            covariate.type = "binary", covariate.eff.mean = 0,
            confounder.type = 'none', depth.mu = 10000, depth.conf.factor = 0)
```

### Arguments
- **`otu.tab`**:     row are taxa, column are samples.   
- **`nOTU`** :    a number, indicates the number of taxa to be simulated.  Default is 500.
- **`nSam`**:     a number, indicates the number of samples to be simulated.  Default is 100.
- **`model`**:    indicates the relationship modeled between covariate and taxa. Options *`loglinear`*, *`heterogeneity`*, *`nonmonotone`*.
- **`diff.otu.pct`**:     a numerical fraction between 0 and 1, indicates the percentage of differential taxa. Default is 0.1.   
- **`diff.otu.direct`**:     Options include *`balanced`*, *`unbalanced`*.  Default is *`balanced`*.
- **`diff.otu.mode`**:    Options include *`abundant`*, *`rare`*, *`both`*. *`abundant`* means differential taxa come from top 1/4 abundant taxa, *`rare`* means differential taxa come from tail 1/4 abundant taxa, *`both`* means half of the differential taxa come from top 1/4 abundant taxa, while the rest half come from tail 1/4 abundant taxa.  Default is *`abundant`*.
- **`include.top.otu`**: *`TRUE`* or *`FALSE`* will be used, respectively, indicates top n abundant taxa are differential taxa or not. Default is *`FALSE`*.
- **`k.top.otu`**: a number to indicates top n number of taxa are will be included in differential taxa. Default is 5.
- **`covariate.eff.mean`**:     effect size to be created between case group and control group. Default is 0.
- **`confounder.type`**:     covariate confounder. Options include *`none`*,  *`continuous`*, *`binary`*, *`both`*.   Default is *`none`*.
- **`depth.mu`**:     mu to generate the sequencing depth.  Default is 10000.
- **`depth.conf.factor`**:    a number indicates depth confounding strength. If 0 is given, no depth confounding will be generated. Default is 0.


### Value
a list with components:
- **`otu.tab.sim`**:    a data.frame of simulated taxa.
- **`covariate`**:    a vector of simulated covariate.
- **`diff.otu.ind`**:   a vector of *`TRUE`* or *`FALSE`* indicates the truth of simulated differential taxa, *`TRUE`* means this taxon is differential taxon, while *`FALSE`* means this taxon is not differential taxon.
- **`otu.names`**:    a vector of simulated taxa names.

### Simulation example
```
source('code/SimulationEvaluation/Simulation.R')
load('data/SimulationEvaluation/Vaginal.RData')
Sim.obj <- SimulateSeq(otu.tab = otu.tab, nOTU = 200, nSam = 50,
                      covariate.eff.mean = 0, model = 'loglinear', nSub = 0,
                      depth.conf.factor = 'none', covariate.type = 'binary')
```
