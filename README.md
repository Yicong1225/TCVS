# TCVS
The TCVS package implements the Tree-guided Compositional Variable Selection method to identify outcome-associated components within high-dimensional microbial compositional data. 
It accommodates compositional covariates (e.g., relative abundances of taxa). This method enhances variable selection by incorporating auxiliary knockoff copies of microbiome features,
which are treated as noise, and utilizing hierarchical tree structures inherent to microbial taxa. By integrating these elements, TCVS refines the selection process to achieve more precise outcomes.
This approach is particularly effective in enhancing selection accuracy, proving advantageous when the signal strength surpasses minimal thresholds. This dual strategy of leveraging tree structure 
and augmenting with knockoff features enables TCVS to provide robust insights into the complex interactions within microbiome data.

To install the package, download the package from this site to a local drive and install and load the package in R:
```{r download and load}
install.packages("TCVS_1.0.tar.gz", repos=NULL)
library(TCVS)
```
