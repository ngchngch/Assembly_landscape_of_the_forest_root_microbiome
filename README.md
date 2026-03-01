# Assembly Landscape of the Forest Root Microbiome
## Keystone Evaluation Framework Based on "Assembly Landscape"
We propose a novel framework for exploring keystone species or taxa within complex microbiomes. 
The assembly rules of the microbiome are inferred using a landscape analogy from community ecology, specifically through energy landscape analysis (Suzuki et al., 2021).
By extending this approach, we quantify the topographic shifts of the inferred landscape in relation to the abundance of focal species (or taxa). 
In this context, we explicitly conceptualize keystone species (or taxa) as those whose increases or decreases in abundance reorganize the rules governing microbiome assembly.

# Workflows
## Validation using generalized Lotka-Volterra simulations
Initially, we evaluated the performance of the proposed framework using generalized Lotka-Volterra models.

### Example
#### 0. Load Libraries and Original Functions
```
library(ggplot2)
library(dplyr)
library(tidyverse)
library(doParallel)
library(igraph)
library(deSolve)
library(vegan)
library(cluster)
source("packages/functions_for_examples.R")
```
#### 1. Generating community dynamics comprasing 80 species

```
# Simulation
total_sample <- 200
rep =200
nrep_eval <- 3
n.core <- 8

time=100
Nsp=c(80)
connectance=0.2
A_self = -1
```
```
#read parameters of Model1
A <- readRDS("examples/datasets/gLV_simulation/Amatrix_seed0.2.rds")[[1]]
r <- readRDS("examples/datasets/gLV_simulation/rparameter_seed0.2_Nsp80.rds")[[1]]

# %% Run simulation

set = expand.grid(Nsp, connectance, 1:rep)|>
  set_names(c("N", "C", "rep"))|>
  as.matrix()

cluster = makeCluster(n.core)
registerDoParallel(cluster)  

dynamics = foreach(n=1:nrow(set),
                   .packages = c("igraph", "tidyverse", "deSolve", "doParallel"))%dopar%{
                     
                     s = paste(set[n,1:2], collapse="_")
                     N=set[n,1]
                     r_set <- r
                     
                     mat = gen_dyn(time=time, dt=1000,
                                   A=A, r=r_set, N0=runif(N), type="gLV")
                     
                     return( list(dynamics=mat) )
                     
                   }|> set_names(apply(set,1,paste, collapse="_"))

stopCluster(cluster)
```
#### 2. Sample Community Compositions from Simulated Dynamics
```
##create community matrix
non_pertabate = map(dynamics, ~.x[[1]])

sim_set <- matrix(unlist(strsplit(names(non_pertabate),"_")),ncol=3,byrow=TRUE)

mat0 <- NULL
for(i in 1:nrow(sim_set[sim_set[,2]==as.character(connectance),])){
  mat0 <- rbind(mat0,non_pertabate[[which(sim_set[,2]==as.character(connectance))[i]]][sample(1:time,total_sample/rep),])
}

d_bray <- vegdist(mat0,method="bray")

qt <- quantile(mat0, probs = c(0.1,0.9))
bin_th <- seq(qt[1],qt[2], length.out = 30)

cr <- c()
for(i in 1:length(bin_th)){
  cat("\r", i, "/", length(bin_th))
  bth <- bin_th[i]
  mat2 <- mat0
  mat2[mat2<=bth] <- 0
  mat2[mat2>0] <- 1
  d_jac <- vegdist(mat2,method="jaccard")
  cr[i] <- cor(d_bray, d_jac,method = "kendall")
}

max_cr <- which.max(cr)
best_bin_th <- bin_th[max_cr[order(max_cr)[1]]]

mat2 <- mat0
mat2[mat2<=best_bin_th] <- 0
mat2[mat2>0] <- 1
```
#### 3. Evaluate Species-Specific Influences on Community Dynamics
To make this process feasible on a local computer, we replicate this evaluation three times.
```
cluster = makeCluster(n.core)
registerDoParallel(cluster)  

kness = foreach(n=1:nrep_eval,
                .packages = c("igraph", "tidyverse", "deSolve", "doParallel"))%dopar%{
                  
                  diff <- NULL
                  com_wis <- NULL; com_wois <- NULL; b_com <- NULL
                  
                  for(i in 1:N){
                    init <- c(runif(N-1),0.5)[order(c(c(1:N)[-i],i))]
                    com_woi = gen_dyn(time=time, dt=1000, A=A[-i,-i], r=r[-i],
                                      N0=init[-i], type="gLV")
                    com_wi = gen_dyn(time, dt=1000, A, r,
                                     init, type="gLV")[,-i]
                    
                    abruptness_diffstart = calc_aburptness(before=com_woi,
                                                           after=com_wi,
                                                           ab_th = best_bin_th)
                    
                    names(abruptness_diffstart) <- paste0(names(abruptness_diffstart),"_diffstart")
                    
                    diff = rbind(diff, c(set[n,], id=i, 
                                         abruptness_diffstart[1:2]))
                    
                  }
                  
                  keystoneness = 
                    data.frame(diff,r=r) |>
                    left_join(
                      calc_centrality(A),
                      by="id"
                    )
                  # 
                  return( list(keystoneness=keystoneness) )
                  
                }

stopCluster(cluster)

keystonness = map_df(kness, ~.x[[1]])

mBC2 <- foreach(i = 1:Nsp, .combine="c")%do% {
  return(mean(keystonness[keystonness$id == i,
                          "BC_diffstart"],na.rm=TRUE))
}


mJC2 <- foreach(i = 1:Nsp, .combine="c")%do% {
  return(mean(keystonness[keystonness$id == i,
                          "Jaccard_diffstart"],na.rm=TRUE))
}

key_res <- data.frame(sp=1:Nsp,
                      mBC_diffstart=mBC2,
                      mJC_diffstart=mJC2
)
```
#### 4. Energy Landscape Analysis
We evaluated species-specific influence on inferred assembly landscape topography. This process is too computationally intensive for a local computer; thus, we present a version with a reduced number of iterations in the parameter fitting and randomization process.
```
seed <- 12
n_itr <- 16
set.seed(seed)

cluster = makeCluster(threads)
registerDoParallel(cluster)

mdp <- foreach(sp=1:80, .packages = c("rELA","doParallel"))%dopar%{
  bin0 <- as.matrix(mat2[,-sp])
  bin <- bin0[,colMeans(bin0)>0.1&colMeans(bin0)<0.9]
  res <- Quant_landshift(bin=bin, #binary community matrix (samples x species)
                         abundance_focal = mat0[,sp], #abundance of the focal species
                         qt_seq = c(0.5), # abundance quantiles to evaluate
                         itr_fit=n_itr, qth=10^-5,SS.itr=150000 #ELA parameters
                         )
  return(res) 
}

stopCluster(cluster)
```
Null model simulations were also performed using permutations of the focal species' abundance.
```
nrandamization <- 1000

stdDtop <- NULL
for(sp in 1:80){
  cat(sprintf("%s\n",sp))
  
  mdp_rand <- foreach(rand=1:nrandamization, .packages = c("rELA","doParallel"))%dopar%{
    bin0 <- as.matrix(mat2[,-sp])
    bin <- bin0[,colMeans(bin0)>0.1&colMeans(bin0)<0.9]
    res <- Quant_landshift(bin=bin, #binary community matrix (samples x species)
                           abundance_focal = mat0[sample(nrow(mat0)),sp], #abundance of the focal species
                           qt_seq = c(0.5), # abundance quantiles to evaluate
                           itr_fit=n_itr, qth=10^-5,SS.itr=150000 #ELA parameters
    )
    return(res) 
  }
  
  mdpr_all <- do.call(rbind,mdp_rand)
  
  z_dtopo <- (mdp[mdp$ra=="perc50","d_land"]-mean(mdpr_all[mdpr_all$ra=="perc50","d_land"]))/sd(mdpr_all[mdpr_all$ra=="perc50","d_land"])
  p_dtopo <- sum(mdpr_all[mdpr_all$ra=="perc50","d_land"]>=mdp[mdp$ra=="perc50","d_land"])/nrandamization
  
  stdDtop <- rbind(stdDtop,data.frame(sp=sp, z_dtopo=z_dtopo, p_dtopo=p_dtopo))
  cat("|\n")
}
```
#### 5. Comparing the Proposed Index and Community-Scale Influence of Each Species
```
infulences <- merge(stdDtop,key_res,by="sp")
cor.test(influences$z_dtopo,influences$mBC_diffstart,method="spearman")
plot(influences$z_dtopo,influences$mJC_diffstart)

cor.test(influences$z_dtopo,influences$mJC_diffstart,method="spearman")
plot(influences$z_dtopo,influences$mJC_diffstart)
```
### Scripts
**Analysis in the SuperComputer System**
working_directory_in_supercomputer/Script/05_01_gLV_simulation_Nsp80_260201.R
working_directory_in_supercomputer/Script/05_01_00_summarize_gLV_simulation_Nsp80_260201.R
working_directory_in_supercomputer/Script/05_02_ELA_withRA_multistep_eachseed_Nsp80_260201.R
working_directory_in_supercomputer/Script/05_03_randELA_withRA_4step_Nsp80_260201.R
working_directory_in_supercomputer/Script/05_04_summarize_ELA_withRA_4step_Nsp80_260202.R
working_directory_in_supercomputer/Script/05_05_Zconv_ELA_withRA_4step_Nsp80_260204.R
working_directory_in_supercomputer/Script/05_06_merge_result_eachseed_Nsp80_260131.R

## Bioinfomatics
We combined the root-tip fungal community datasets described in our previous study ([Noguchi and Toju *et al.*, 2024](https://doi.org/10.1002/ecm.1469)) with newly obtained prokaryotic data. The sequncing outputs of six Miseq runs were processed respectively (bioinfomatics pipelines were described in the corresponding "RunXX" directories and these outputs are in the directory "Base_data/Bioinfomatics/seqtab") and converted to a sample-OTU matrixusing the scripts in "Base_data/Bioinfomatics/Script".

We then applied coverage-based rarefaction to the 1,270 root fungal and prokaryotic community data-sets.

### Script
**Merge multiple sequence outputs**
Base_data/Bioinfomatics/Script/01_merge_and_decontm_2libararies.R
**Rarefaction**
Script_in_local_computer/01_LOO_covrarefy.R

## Energy Landscape Analysis

 In the family-level taxonomic composition matrix, relative read counts for each family were binarized using the threshold. To make the subsequent energy landscape analysis computationally feasible, we prioritized families by their contribution to overall community structure as measured by PerMANOVA (*R²*). Among candidate family sets ranked by *R²*, we selected the set whose binarized pattern best matched the abundance-based community structure. Energy landscape analysis ([Suzuki *et al.*, 2021](https://doi.org/10.1002/ecm.1469)) was then performed using this selected family set together with host plant genera (encoded as dummy variables) as explanatory variables.

### Example
#### 0. 



### Scripts
**Energy landscape analysis in the SuperComputer System**
working_directory_in_supercomputer/Script/02_06_ELA.R
working_directory_in_supercomputer/Script/02_07_assemblygraph_onlyBasin.R

**Some graphics**
These script run after the analyses about flow diagrams of energy landscape.
Script_in_local_computer/03_11_SSheatmap_fullELA_recolor_250501.R Script_in_local_computer/03_12_DG_fullELA_recolor.R

## Statistical Inference of Keystone Taxa


Starting from the original data matrix with OTUs annotated as the focal genus removed, we performed coverage-based rarefaction. Binarization used the same family set as in the energy landscape analysis described above. In parallel, we rarefied the full community matrix and applied a centered log-ratio (CLR) transformation to genus-level compositions. 


We then performed energy landscape analysis including host plant genera (dummy variables) and the CLR-transformed relative abundance of the focal genus as external variables. "Keystoneness" indices were computed by comparing energy landscapes inferred under two conditions: (1) without the focal genus and (2) with the focal genus fixed at representative abundances (25%, 50%, and 75% quantiles of its observed relative abundance), using community assembly simulations.

### Script
**Community assembly simulations**
working_directory_in_supercomputer/Script/03_01_ELA_withRA_4step.R

**Null model simulations**
working_directory_in_supercomputer/Script/03_02_randELA_withRA_4s_fixPS_1_3000.R
working_directory_in_supercomputer/Script/03_02_randELA_withRA_4s_fixPS_3001_6000.R
working_directory_in_supercomputer/Script/03_02_randELA_withRA_4s_fixPS_6001_9000.R
working_directory_in_supercomputer/Script/03_02_randELA_withRA_4s_fixPS_9001_10500.R
working_directory_in_supercomputer/Script/03_02_summarize_randELA_withRA_fixP_250306.R

**Evaluation variances of the keystoneness metrics**
working_directory_in_supercomputer/Script/03_03_ELA_withRA_4step_rep.R
working_directory_in_supercomputer/Script/03_04_landscape_change_repuroducibility.R

**Flow diagrams**
working_directory_in_supercomputer/Script/03_05_states_flow_diagram.R

**graphics**

Script_in_local_computer/03_06_Zhistgram_250321.R Script_in_local_computer/03_07_graphics_Zconv_landchanges_biplot_250312.R Script_in_local_computer/03_08_Zeven_abundance_occurence_250507.R Script_in_local_computer/03_08_Zland_abundance_occurence_250507.R Script_in_local_computer/03_10_02_graphics_Fullstates_flow_Spl_250813.R 


## Additional analyses
### Scripts

Script_in_local_computer/02_08_hostpreference_Family.R
Script_in_local_computer/04_rarefaction_barplot.R

## Repository Contents


- `Base_data/` — Raw datasets used in this study.
- `Output/` — Results produced on the local computer (some large folders excluded).
- `Output_supercomputer/` — Results produced on the supercomputer (some large folders excluded).
- `Script_in_local_computer/` — R scripts used for analyses on the local computer.
- `packages/` — R packages and custom source code used for local analyses.
- `working_directory_in_supercomputer/` — Working directory structure used on the supercomputer, containing the analysis scripts (`Script/`) and additional data prepared on the local computer (`Import_data/`, `color/`). Note: `Base_data/` and `packages/` referenced here are not duplicated in this repository.

