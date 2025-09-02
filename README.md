# CARLIS
 Covariate-assisted replicability analysis for GWAS

```R
## Install dependency packages if necessary
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")

## Install CARLIS
install.packages("devtools")
devtools::install_github("YanLi15/CARLIS")

## Load CARLIS
library(CARLIS)
```

