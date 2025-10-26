# Assembly_landscape_of_the_forest_root_microbiome
In this repository, we showed "assembly landscape" concept and algorisms of keystone inference based on it.

# Relationship between clasical stability landscape and empirically reconstructed assembly landscapes
The stability landscape concept illustrates the relationship between community structure and stability (left), whereas assembly landscapes, inferred from empirical data, describe the statistical relationship between community states and their probabilities of observation (right). Both convergence toward true attractors and prolonged transients around a “ghost attractor” on a conceptual stability landscape can result in clusters of observations in empirical datasets. In statistical analyses of assembly landscapes, assembly toward true attractors is represented by deep basins, while long transients around ghost attractors can be inferred as shallow basins. Assembly landscapes estimated by the statistical framework of energy landscape analysis (described below) are referred to as energy landscapes.

![Stability landscape concept and empirically reconstructed assembly landscapes.](figures/Fig1a_assemblylandscape_concept.png)

# empirical estimation using energy landscape analysis
Based on the Ising model in statistical physics, the probability of observing a given community state (membership) can be expressed as a linear combination of effects from implicit (unobserved) factors, explicit (observed) factors, and associations among species/taxa (Suzuki et al., 2021 [![DOI](https://img.shields.io/badge/DOI-10.1002%2Fecm.1469-blue.svg)](https://doi.org/10.1002/ecm.1469)). 

  The lower the "energy" index, the more probable the corresponding community state is. The landscape topography is assumed to change depending on background environmental conditions.

![Model of energy landscape analysis](figures/Fig1b_model.png)

Based on the maximum-entropy fitting of the model parameters, the "energy" of each community state is inferred. The topography of the "energy landscape" is assumed to change depending on background environmental conditions. Community assembly is assumed to proceed toward basin bottoms within the "assembly graphs," in which community states that differ by the presence or absence of a single species are connected by one step.

![Statistically inferred assembly landscapes](figures/Fig1c_inferred_landscapes.png)

# Keystone concept based on assembly landscape reorganization
Given that assembly landscapes represent community assembly rules, a “keystone microbes” can be defined as one whose changes in abundance drastically change the landscape architecture. The extent of the landscape reorganization is quantitatively evaluated in terms of changes in the composition and distribution of basin bottoms (Δtopography; left) and changes in the evenness of basin distributions (Δevenness; right) along the abundance gradient of a focal species/taxon.

![Large shift in assembly landscape](figures/Fig1de_landchange_concept.png)

The dependence of energy landscape topography on the abundance of a focal species/taxon is inferred (detailed workflow are described below). Changes in assembly landscape architecture are then evaluated by comparing the topography inferred in the absence of the focal species/taxon with that estimated under the assumption of an intermediate abundance level (the median abundance across samples in which the genus was present). Specifically, across 20,000 community assembly simulations on the inferred landscapes, mean Jaccard distance of corresponding destinations (basin bottoms) calculated (Δtopography). Likewise, changes in the evenness of basin distributions were calculated as the difference in Shannon entropy between the two landscapes under different abundance scenarios (Δevenness).

![Keystone indexes](figures/Fig1g_indexes.png)

# Workflows in this study
## Energy landscape analysis
We performed coverage-based rarefaction for the root samples. In a matrix of family-level taxonomic compositions, relative read counts of each family were binarized with a threshold described. To make the subsequent statistical analysis (energy landscape analysis) computationally feasible, we prioritized and selected families based on their contribution to overall community structure, evaluated by PerMANOVA (R²; Supplementary Tables S1 and S2). The family set showing the highest correspondence between its binarized pattern and the abundance-based community structure was selected among candidate sets determined by their R² values. Energy landscape analysis was then conducted using the family-set with host plant genera (encoded as dummy variables) as explanatory variables.

![Energylandscape estimation](figures/FigS1a_ela_prep.png)

## Statistical inference of keystone taxa
From the original data matrix excluding OTUs annotated as a focal genus, coverage-based rarefaction was conducted. Binarization was applied to the same set of families as in above energy landscape analysis. In parallel, the whole community matrix was rarefied, and the genus-level compositions were subjected to CLR transformation. Energy landscape analysis was performed using host plant genera (dummy variables) and CLR-transformed relative read counts (+1) of the focal genus as external factors. “Keystoneness” indices were calculated based on comparisons between the inferred energy landscapes under conditions without the focal genus and those where it was present at the representative abundances (25%, 50%, and 75% quantiles of its relative abundance) in community assembly simulations.

![statistical inference of keystone taxa](figures/FigS1b_explore_keytstone.png)
