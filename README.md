# Microbial genomics extensions / extra functions in R


### Dependencies

* dplyr
* ggplot2
* stringr
* gplots
* RColorBrewer  
* ggtree
* phytools
* ape
* castor
* purrr
* glue
* tibble
* seqinr
* genoplotR
* vegan
* UpSetR
* cowplot
* Peptides


# Install dependencies

```R
install.packages(c("dplyr", "ggplot2", "stringr", "gplots", "RColorBrewer", "phytools", "ape", "castor", "purrr", "glue", "tibble", "seqinr","genoPlotR","vegan", "UpSetR", "cowplot", "Peptides"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ggtree"))
```


# Install micro.gen.extra

```R
install.packages("devtools")
library(devtools)
install_github("djw533/micro.gen.extra")
```
