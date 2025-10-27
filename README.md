# Assembly Landscape of the Forest Root Microbiome

# Workflows

## Energy Landscape Analysis


We applied coverage-based rarefaction to the 1,270 root fungal and prokaryotic community data-sets (Fungal community; [Noguchi and Toju *et al.*, 2024](https://doi.org/10.1002/ecm.1469)). In the family-level taxonomic composition matrix, relative read counts for each family were binarized using the threshold described in Figure 6. To make the subsequent energy landscape analysis computationally feasible, we prioritized families by their contribution to overall community structure as measured by PerMANOVA (*R²*). Among candidate family sets ranked by *R²*, we selected the set whose binarized pattern best matched the abundance-based community structure. Energy landscape analysis was then performed using this selected family set together with host plant genera (encoded as dummy variables) as explanatory variables.

*Fig. 6| Input data for energy landscape analysis.*


## Statistical Inference of Keystone Taxa


Starting from the original data matrix with OTUs annotated as the focal genus removed, we performed coverage-based rarefaction. Binarization used the same family set as in the energy landscape analysis described above. In parallel, we rarefied the full community matrix and applied a centered log-ratio (CLR) transformation to genus-level compositions. 


We then performed energy landscape analysis including host plant genera (dummy variables) and the CLR-transformed relative abundance of the focal genus as external variables. "Keystoneness" indices were computed by comparing energy landscapes inferred under two conditions: (1) without the focal genus and (2) with the focal genus fixed at representative abundances (25%, 50%, and 75% quantiles of its observed relative abundance), using community assembly simulations.

*Fig. 7| Input data for the keystone exploration.*


## Repository Contents


- `Base_data/` — Raw datasets used in this study.
- `Output/` — Results produced on the local computer (some large folders excluded).
- `Output_supercomputer/` — Results produced on the supercomputer (some large folders excluded).
- `Script_in_local_computer/` — R scripts used for analyses on the local computer.
- `packages/` — R packages and custom source code used for local analyses.
- `working_directory_in_supercomputer/` — Working directory structure used on the supercomputer, containing the analysis scripts (`Script/`) and additional data prepared on the local computer (`Import_data/`, `color/`). Note: `Base_data/` and `packages/` referenced here are not duplicated in this repository.

