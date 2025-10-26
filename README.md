# Assembly_landscape_of_the_forest_root_microbiome
In this repository, we showed "assembly landscape" concept and algorisms of keystone inference based on it.

# Relationship between clasical stability landscape and empirically reconstructed assembly landscapes
The stability landscape concept illustrates the relationship between community structure and stability (left), whereas assembly landscapes, inferred from empirical data, describe the statistical relationship between community states and their probabilities of observation (right). Both convergence toward true attractors and prolonged transients around a “ghost attractor” on a conceptual stability landscape can result in clusters of observations in empirical datasets. In statistical analyses of assembly landscapes, assembly toward true attractors is represented by deep basins, while long transients around ghost attractors can be inferred as shallow basins. Assembly landscapes estimated by the statistical framework of energy landscape analysis (described below) are referred to as energy landscapes.

![Stability landscape concept and empirically reconstructed assembly landscapes.](figures/Fig1a_assemblylandscape_concept.png)

# empirical estimation using energy landscape analysis
Based on the Ising model in statistical physics, the probability of observing a given community state (membership) can be expressed as a linear combination of effects from implicit (unobserved) factors, explicit (observed) factors, and associations among species/taxa (Suzuki et al., 2021 [![DOI](10.1002/ecm.1469)](10.1002/ecm.1469)). 
The lower the "energy" index, the more probable the corresponding community state is. The landscape topography is assumed to change depending on background environmental conditions.

![Model of energy landscape analysis](figures/Fig1b_model.png)

![Statistically inferred assembly landscapes](figures/Fig1c_inferred_landscapes.png)

# Keystone concept based on assembly landscape reorganization
![Large shift in assembly landscape](figures/Fig1de_landchange_concept.png)

![Keystone indexes](figures/Fig1g_indexes.png)

# Workflows in this study
![Energylandscape estimation](figures/FigS1a_ela_prep.png)

![statistical inference of keystone taxa](figures/FigS1b_explore_keytstone.png)
