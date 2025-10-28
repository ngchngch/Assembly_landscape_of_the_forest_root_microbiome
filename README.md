# Assembly Landscape of the Forest Root Microbiome

# Workflows

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

